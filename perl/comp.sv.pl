#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  comp.sv.pl - 

=head1 SYNOPSIS
  
  comp.sv.pl [-help] [-qry qry-genome] [-tgt tgt-genome]

  Options:
    -h (--help)   brief help message
    -q (--qry)    qry genome
    -t (--tgt)    tgt genome

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use File::Path qw/make_path remove_tree/;
use File::Basename;
use Tabix;
use Bio::DB::Fasta;
use Common;
use Location;
use Gtb;
use Gal;
use Bed;
use List::Util qw/min max sum/;

my $help_flag;
my ($qry, $tgt) = qw/HM004 HM101/;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "qry|q=s"  => \$qry,
  "tgt|t=s"  => \$tgt,
) or pod2usage(2);
pod2usage(1) if $help_flag;

my @qrys = qw/
  HM058 HM125 HM056 HM129 HM060
  HM095 HM185 HM034 HM004 HM050 
  HM023 HM010 HM022 HM324 HM340
/;
@qrys = qw/HM004/;

my $dir = "$ENV{'misc3'}/$qry\_$tgt/31_sv";
-d $dir || make_path($dir);
chdir $dir || die "cannot chdir to $dir\n";

my $tdir = "$ENV{'genome'}/$tgt";
my $qdir = "$ENV{'genome'}/$qry";
my $cdir = "$ENV{'misc3'}/$qry\_$tgt/23_blat";

my $fgt = "$tdir/51.gtb";
my $fgq = "$qdir/51.gtb";
my $tdb = Bio::DB::Fasta->new("$tdir/51.fas");
my $qdb = Bio::DB::Fasta->new("$qdir/51.fas");

my $tgene = Tabix->new(-data => "$tdir/51.tbl.gz");
my $qgene = Tabix->new(-data => "$qdir/51.tbl.gz");

cat_var($tdir, $qdir, $cdir, "01.stb", "01.tlc");
runCmd("sort.header.pl -f stb -i 01.stb -o 02.sort.stb");
runCmd("sort.header.pl -f tlc -i 01.tlc -o 02.sort.tlc");
refine_tlc("02.sort.tlc", "05.refine.tlc");
#sv_stb2bed("02.sort.stb", "08");
#sv_tlc2bed("05.refine.tlc", "08");
stb_ovlp_cds("02.sort.stb", "11.cds.tbl");

sub cat_var {
  my ($tdir, $qdir, $dir, $fos, $fot) = @_;
  my $fi = "$dir/31.9/idm";
  open(my $fhi, "<$fi") || die "cannot read $fi\n";
  open(my $fhs, ">$fos") || die "cannot write $fos\n";
  open(my $fht, ">$fot") || die "cannot write $fot\n";
  print $fhs join("\t", qw/tchr tbeg tend btbeg btend type 
    qchr qbeg qend bqbeg bqend/)."\n";
  print $fht join("\t", qw/tdchr tdbeg tdend tichr tibeg 
    qdchr qdbeg qdend qichr qibeg type/)."\n";

  my $gta = Tabix->new(-data => "$dir/31.5/gax.gz");
  my $gqa = Tabix->new(-data => "$dir/41.5/gax.gz");
  my $gtr = Tabix->new(-data => "$dir/31.9/gax.gz");
  my $gqr = Tabix->new(-data => "$dir/41.9/gax.gz");

  my $gapt = Tabix->new(-data => "$tdir/16.gap.bed.gz");
  my $gapq = Tabix->new(-data => "$qdir/16.gap.bed.gz");

  while( <$fhi> ) {
    chomp;
    next if /(^\#)|(^\s*$)/;
    my ($tid, $tbeg, $tend, $tsrd, $qid, $qbeg, $qend, $qsrd, $cid, $lev) 
      = split "\t";
    $lev == 1 || next;
    my $tlen = $tend - $tbeg - 1;
    my $qlen = $qend - $qbeg - 1;
   
    my ($flagt, $flagq) = (0, 0);
    if($tlen > 0) {
      my $ary = read_gap($gapt, $tid, $tbeg + 1, $tend - 1);
      $flagt = @$ary == 0 ? 1 : 0;
      $flagq = @$ary == 0 ? $flagq : 0;
    }
    if($qlen > 0) {
      my $ary = read_gap($gapq, $qid, $qbeg + 1, $qend - 1);
      $flagq = @$ary == 0 ? 1 : 0;
      $flagt = @$ary == 0 ? $flagt : 0;
    }
    if($flagt) {
      my $stats = idm_cat($tid, $tbeg+1, $tend-1, $lev, $gtr, $gta, $gqr);
      for (@$stats) {
        my ($tb, $te, $type, $qi, $qb, $qe) = @$_;
        if($type eq "del") {
          print $fhs join("\t", $tid, $tb, $te, $tbeg, $tend, 
            'DEL', $qid, $qbeg, $qbeg, $qbeg, $qend)."\n";
        } elsif($type eq "cnv") {
          print $fhs join("\t", $tid, $tb, $te, $tbeg, $tend, 
            'CNL', $qid, $qbeg, $qbeg, $qbeg, $qend)."\n";
        } else {
          $type eq "tlc" || die "unknown type: $type\n";
          print $fht join("\t", $tid, $tb, $te, '', '', 
            $qi, $qb, $qe, $qid, $qbeg, 'tTLC')."\n";
        }
      }
    }
    if($flagq) {
      my $stats = idm_cat($qid, $qbeg+1, $qend-1, $lev, $gqr, $gqa, $gtr);
      for (@$stats) {
        my ($qb, $qe, $type, $ti, $tb, $te) = @$_;
        if($type eq "del") {
          print $fhs join("\t", $tid, $tbeg, $tbeg, $tbeg, $tend, 
            'INS', $qid, $qb, $qe, $qbeg, $qend)."\n";
        } elsif($type eq "cnv") {
          print $fhs join("\t", $tid, $tbeg, $tbeg, $tbeg, $tend, 
            'CNG', $qid, $qb, $qe, $qbeg, $qend)."\n";
        } else {
          $type eq "tlc" || die "unknown type: $type\n";
          print $fht join("\t", $ti, $tb, $te, $tid, $tbeg, 
            $qid, $qb, $qe, '', '', 'qTLC')."\n";
        }
      }
    }
  }
  close $fhi;
  close $fhs;
  close $fht;
}
sub idm_cat {
  my ($id, $beg, $end, $lev, $gr, $ga, $gqr) = @_;
  my @stats;
  if($end - $beg + 1 < 30) {
    push @stats, [$beg, $end, 'del'];
    return \@stats;
  }
  
  my $hr = read_chain($id, $beg, $end, $gr);
  my $lev2 = min( map {$hr->{$_}->[7]} keys(%$hr) );
  my @cids = grep {$hr->{$_}->[7] == $lev2} keys(%$hr);
  my $loct = [];
  for my $cid (@cids) {
    my ($ti, $tb, $te, $qi, $qb, $qe, $len) = @{$hr->{$cid}};
    my $hx = read_chain($qi, $qb, $qe, $gqr);
    my $leno = check_ovlp($ti, $tb, $te, $hx);
    $leno / ($te - $tb + 1) >= 0.8 || next;
    push @stats, [$tb, $te, 'tlc', $qi, $qb, $qe];
    push @$loct, [$tb, $te];
  }
  my $len1 = locAryLen($loct);
  my $len2 = locAryLen(posMerge($loct));
  $len1 == $len2 || die "$id:$beg-$end error\n".Dumper($hr)."\n";
 
  my $locp = [];
  my $hp = read_chain($id, $beg, $end, $ga);
  for my $cid (keys(%$hp)) {
    my ($ti, $tb, $te, $qi, $qb, $qe, $len, $lev2) = @{$hp->{$cid}};
    push @$locp, [$tb, $te];
  }
  $locp = posMerge($locp);
  
#  $beg =~ /^1074/ && die "$id:$beg-$end\n".Dumper(read_gax($ga, $id, $beg, $end)); 
  my ($locd) = posDiff([[$beg, $end]], $locp);
  for (@$locd) {
    my ($tb, $te) = @$_;
    push @stats, [$tb, $te, 'del'];
  }
  my $locu = posSubtract($locp, $loct);
  for (@$locu) {
    my ($tb, $te) = @$_;
    push @stats, [$tb, $te, 'cnv'];
  }
  return \@stats;
}
sub refine_tlc {
  my ($fi, $fo) = @_;
  my $t = readTable(-in => $fi, -header => 1);
  stb_tlc2bed($t, "$fo.1.bed");
  runCmd("mergeBed -i $fo.1.bed -c 4 -o collapse > $fo.2.bed");
  
  open(my $fhb, "<$fo.2.bed") or die "cannot read $fo.2.bed\n";
  my @idxs_rm;
  while(<$fhb>) {
    chomp;
    my ($chr, $beg, $end, $idxstr) = split "\t";
    my @idxs = split(",", $idxstr);
    my $h;
    @idxs > 1 || next;
    for my $i (0..$#idxs-1) {
      my ($idx1, $idx2) = ($idxs[$i]-1, $idxs[$i+1]-1);
      my ($flag, $idxk, $idxr) = check_rec_tlc($t, $idx1, $idx2);
      $flag == 1 || next;
      !exists $h->{$idxk} || die "$idxk in >1 rec-tlc\n";
      $h->{$idxk} = $idxr;
      $t->setElm($idxk, "tichr", $t->elm($idxr, 'tichr'));
      $t->setElm($idxk, "tibeg", $t->elm($idxr, 'tibeg'));
      $t->setElm($idxk, "type", 'TLC');
      push @idxs_rm, $idxr;
    }
  }
  close $fhb;

  $t->delRows(\@idxs_rm);
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho $t->tsv(1);
  close $fho;

  runCmd("rm $fo.*.bed");
}
sub stb_tlc2bed {
  my ($ti, $fo) = @_;
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  for my $i (0..$ti->lastRow) {
    my ($tdchr, $tdbeg, $tdend) = $ti->row($i);
    print $fho join("\t", $tdchr, $tdbeg, $tdend, $i+1)."\n";
  }
  close $fho;
}
sub check_rec_tlc {
  my ($t, $idx1, $idx2) = @_;
  my ($type1, $type2) = map {$t->elm($_, "type")} ($idx1, $idx2);
  my ($idxk, $idxr) = ($type1 eq "qTLC" && $type2 eq "tTLC") ? 
    ($idx2, $idx1) : ($type1 eq "tTLC" && $type2 eq "qTLC") ? 
    ($idx1, $idx2) : (undef, undef);
  return 0 unless defined $idxk;
  my @ary1 = $t->row($idxk);
  my @ary2 = $t->row($idxr);
  my ($ti1, $tb1, $te1, $qi1, $qb1, $qe1) = @ary1[0..2,5..7];
  my ($ti2, $tb2, $te2, $qi2, $qb2, $qe2) = @ary2[0..2,5..7];
  my ($tl1, $ql1) = ($te1 - $tb1 + 1, $qe1 - $qb1 + 1);
  my ($tl2, $ql2) = ($te2 - $tb2 + 1, $qe2 - $qb2 + 1);
  return 0 if $ti1 ne $ti2 || $qi1 ne $qi2;
  my $to = min($te1, $te2) - max($tb1, $tb2) + 1;
  my $qo = min($qe1, $qe2) - max($qb1, $qb2) + 1;
  return 0 if $to/$tl1 < 0.9 || $to/$tl2 < 0.9 ||
    $qo/$ql1 < 0.9 || $qo/$ql2 < 0.9;
  return (1, $idxk, $idxr);
}
sub read_chain {
  my ($chr, $beg, $end, $gax) = @_;
  my $h = {};
  my $ary = read_gax($gax, $chr, $beg, $end);
  for (@$ary) {
    my ($cid, $tid, $tb, $te, $tsrd, $qid, $qb, $qe, $qsrd, $lev) = @$_;
    $h->{$cid} ||= [$tid, $tb, $te, $qid, $qb, $qe, 0, $lev];
    $h->{$cid}->[1] = min($tb, $h->{$cid}->[1]); 
    $h->{$cid}->[2] = max($te, $h->{$cid}->[2]); 
    $h->{$cid}->[4] = min($qb, $h->{$cid}->[4]); 
    $h->{$cid}->[5] = max($qe, $h->{$cid}->[5]); 
    $h->{$cid}->[6] += $te - $tb + 1;
  }
  return $h;
}
sub check_ovlp {
  my ($tid, $tbeg, $tend, $h) = @_;
  my $loc = [];
  for my $cid (keys(%$h)) {
    my ($qi, $qb, $qe, $ti, $tb, $te, $len) = @{$h->{$cid}};
    $ti eq $tid || next;
    push @$loc, [$tb, $te];
  }
  my ($ref, $olen) = posOvlp([[$tbeg, $tend]], $loc);
  return $olen;
}

sub sv_tlc2bed {
  my ($fi, $fo) = @_;
  open(my $fhg, ">$fo.tlc.gan.bed") or die "cannot write $fo\n";
  open(my $fhl, ">$fo.tlc.los.bed") or die "cannot write $fo\n";
  my $t = readTable(-in => $fi, -header => 1);
  for my $i (0..$t->lastRow) {
    my ($tdchr, $tdbeg, $tdend, $tichr, $tibeg, 
        $qdchr, $qdbeg, $qdend, $qichr, $qibeg, $type) = $t->row($i);
    print $fhg join("\t", $qdchr, $qdbeg-1, $qdend)."\n";
    print $fhl join("\t", $tdchr, $tdbeg-1, $tdend)."\n";
  }
  close $fhg;
  close $fhl;
}
sub sv_stb2bed {
  my ($fi, $fo) = @_;
  open(my $fhi, ">$fo.ins.bed") or die "cannot write $fo\n";
  open(my $fhg, ">$fo.gan.bed") or die "cannot write $fo\n";
  open(my $fhd, ">$fo.del.bed") or die "cannot write $fo\n";
  open(my $fhl, ">$fo.los.bed") or die "cannot write $fo\n";
  my $t = readTable(-in => $fi, -header => 1);
  for my $i (0..$t->lastRow) {
    my ($tchr, $tbeg, $tend, $btbeg, $btend, $type,
      $qchr, $qbeg, $qend, $bqbeg, $bqend) = $t->row($i);
    if($type eq "INS") {
      print $fhi join("\t", $qchr, $qbeg-1, $qend)."\n";
    } elsif($type eq "CNG") {
      print $fhg join("\t", $qchr, $qbeg-1, $qend)."\n";
    } elsif($type eq "DEL") {
      print $fhd join("\t", $tchr, $tbeg-1, $tend)."\n";
    } elsif($type eq "CNL") {
      print $fhl join("\t", $tchr, $tbeg-1, $tend)."\n";
    } else {
      die "unknonw type: $type\n";
    }
  }
  close $fhi;
  close $fhg;
  close $fhd;
  close $fhl;
}
sub read_cds {
  my ($con, $chr, $beg, $end) = @_;
  my $iter = $con->query($chr, $beg - 1, $end);
  my @ary;
  return \@ary if ! $iter->get();
  while (my $line = $con->read($iter)) {
    my @ps = split("\t", $line);
    my ($chr, $beg1, $end1, $srd, $id, $type, $cat) = @ps;
    $type eq "cds" || next;
    push @ary, [$chr, max($beg, $beg1), min($end, $end1), $srd, $id];
  }
  return \@ary;
}
sub stb_ovlp_cds {
  my ($fi, $fo) = @_;

  my ($hq, $ht);
  my $tq = readTable(-in => $fgq, -header => 1);
  my $tt = readTable(-in => $fgt, -header => 1);
  for my $i (0..$tq->lastRow) {
    my ($id, $clocS) = map {$tq->elm($i, $_)} qw/id cloc/;
    my $len = locAryLen(locStr2Ary($clocS));
    $hq->{$id} = $len;
  }
  for my $i (0..$tt->lastRow) {
    my ($id, $clocS) = map {$tt->elm($i, $_)} qw/id cloc/;
    my $len = locAryLen(locStr2Ary($clocS));
    $ht->{$id} = $len;
  }

  my $t = readTable(-in => $fi, -header => 1);
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho join("\t", $t->header, qw/svlen gid glen gslen/)."\n";
  my ($hg, $con, $chr, $beg, $end);
  for my $i (0..$t->lastRow) {
    my ($tchr, $tbeg, $tend, $btbeg, $btend, $type,
      $qchr, $qbeg, $qend, $bqbeg, $bqend) = $t->row($i);
    if($type eq "INS" || $type eq "CNG") {
      ($hg, $con, $chr, $beg, $end) = ($hq, $qgene, $qchr, $qbeg, $qend);
    } elsif($type eq "DEL" || $type eq "CNL") {
      ($hg, $con, $chr, $beg, $end) = ($ht, $tgene, $tchr, $tbeg, $tend);
    } else {
      die "unknown type: $type\n";
    }
    
    my $svlen = $end - $beg + 1;
    my $ary = read_cds($con, $chr, $beg, $end);
    next if(@$ary == 0);
    my $h;
    for (@$ary) {
      my ($chr, $beg, $end, $srd, $id) = @$_;
      $h->{$id} ||= 0;
      $h->{$id} += $end - $beg + 1;
    }
    my @notes;
    for my $id (keys(%$h)) {
      exists $hg->{$id} || die "no len for $id\n";
      my $glen = $hg->{$id};
      my $gslen = $h->{$id};
      print $fho join("\t", $t->row($i), $svlen, $id, $glen, $gslen)."\n";
    }
  }
}

__END__

