#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  idm.refine.pl -

=head1 SYNOPSIS
  
  idm.refine.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input (IDM) file
    -o (--out)    output file
    -t (--tgt)    tgt genome (31.5/gax.gz, 31.9/gax.gz)
    -q (--qry)    qry genome (41.5/gax.gz, 41.9/gax.gz)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use Common;
use Location;
use Gal;
use Bed;
use List::Util qw/min max sum/;

my ($fi, $fo) = ('') x 2;
my ($qry, $tgt) = ('HM004', 'HM101');
my ($fhi, $fho);
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
  "qry|q=s"  => \$qry,
  "tgt|t=s"  => \$tgt,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

sv_tbl2bed("$fo.ins.tbl", "$fo.ins.bed");
sv_tbl2bed("$fo.gan.tbl", "$fo.gan.bed");
sv_tbl2bed("$fo.del.tbl", "$fo.del.bed");
sv_tbl2bed("$fo.los.tbl", "$fo.los.bed");
sv_tlc2bed("$fo.tlc.tbl", "$fo.tlc.gan.bed", "$fo.tlc.los.bed");
sub sv_tbl2bed {
  my ($fi, $fo) = @_;
  open(my $fhi, "<$fi") or die "cannot read $fi\n";
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  while(<$fhi>) {
    chomp;
    /(^chr\s)|(^\s*$)/ && next;
    my ($chr, $beg, $end, $locs) = split "\t";
    my ($c, $b, $e) = parse_locstr($locs);
    print $fho join("\t", $c, $b-1, $e)."\n";
  }
  close $fhi;
  close $fho;
}
sub sv_tlc2bed {
  my ($fi, $fg, $fl) = @_;
  open(my $fhi, "<$fi") or die "cannot read $fi\n";
  open(my $fhg, ">$fg") or die "cannot write $fg\n";
  open(my $fhl, ">$fl") or die "cannot write $fl\n";
  while(<$fhi>) {
    chomp;
    /(^tid)|(^\s*$)/ && next;
    my ($tid, $tbeg, $tend, $qid, $qbeg, $qend) = split "\t";
    print $fhg join("\t", $qid, $qbeg-1, $qend)."\n";
    print $fhl join("\t", $tid, $tbeg-1, $tend)."\n";
  }
  close $fhi;
  close $fhg;
  close $fhl;
}
__END__
cat_var($fi, "$fo.1.ins.tbl", "$fo.1.del.tbl", "$fo.1.gan.tbl", "$fo.1.los.tbl", "$fo.1.tlc.tbl");
runCmd("sort.header.pl -i $fo.1.ins.tbl -o $fo.ins.tbl");
runCmd("sort.header.pl -i $fo.1.del.tbl -o $fo.del.tbl");
runCmd("sort.header.pl -i $fo.1.gan.tbl -o $fo.gan.tbl");
runCmd("sort.header.pl -i $fo.1.los.tbl -o $fo.los.tbl");
runCmd("sort.header.pl -i $fo.1.tlc.tbl -o $fo.2.tlc.tbl");
runCmd("awk 'BEGIN{OFS=\"\\t\"} {if(NR>1) {print \$1, \$2, \$3, NR-1}}' \\
  $fo.2.tlc.tbl > $fo.2.tlc.bed");
runCmd("mergeBed -i $fo.2.tlc.bed -c 4 -o collapse > $fo.3.tlc.bed");
refine_bal("$fo.2.tlc.tbl", "$fo.3.tlc.bed", "$fo.tlc.tbl");
runCmd("rm $fo.[1-3].*");

sub check_rec_tlc {
  my ($t, $idx1, $idx2) = @_;
  my ($type1, $type2) = map {$t->elm($_, "type")} ($idx1, $idx2);
  my ($idxk, $idxr) = ($type1 eq "q" && $type2 eq "t") ? ($idx2, $idx1) :
    ($type1 eq "t" && $type2 eq "q") ? ($idx1, $idx2) : (undef, undef);
  return 0 unless defined $idxk;
  my ($ti1, $tb1, $te1, $qi1, $qb1, $qe1) = $t->row($idxk);
  my ($ti2, $tb2, $te2, $qi2, $qb2, $qe2) = $t->row($idxr);
  my ($tl1, $ql1) = ($te1 - $tb1 + 1, $qe1 - $qb1 + 1);
  my ($tl2, $ql2) = ($te2 - $tb2 + 1, $qe2 - $qb2 + 1);
  return 0 if $ti1 ne $ti2 || $qi1 ne $qi2;
  my $to = min($te1, $te2) - max($tb1, $tb2) + 1;
  my $qo = min($qe1, $qe2) - max($qb1, $qb2) + 1;
  return 0 if $to/$tl1 < 0.9 || $to/$tl2 < 0.9 ||
    $qo/$ql1 < 0.9 || $qo/$ql2 < 0.9;
  return (1, $idxk, $idxr);
}

sub refine_bal {
  my ($fi, $fb, $fo) = @_;
  my $t = readTable(-in => $fi, -header => 1);

  open(my $fhb, "<$fb") or die "cannot read $fb\n";
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
      $t->setElm($idxk, "tdloc", $t->elm($idxr, "tdloc"));
      $t->setElm($idxk, "qiloc", $t->elm($idxr, "qiloc"));
      $t->setElm($idxk, "type", 'r');
      push @idxs_rm, $idxr;
    }
  }
  close $fhb;

  $t->delRows(\@idxs_rm);
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho $t->tsv(1);
  close $fho;
}
sub cat_var {
  my ($fi, $fs, $fd, $fg, $fl, $ft) = @_;
  open(my $fhi, $fi) || die "cannot read $fi\n";
  open(my $fhs, ">$fs") || die "cannot write $fs\n";
  print $fhs join("\t", qw/chr beg end loc tloc qloc/)."\n";
  open(my $fhd, ">$fd") || die "cannot write $fd\n";
  print $fhd join("\t", qw/chr beg end loc tloc qloc/)."\n";
  open(my $fhg, ">$fg") || die "cannot write $fg\n";
  print $fhg join("\t", qw/chr beg end loc tloc qloc/)."\n";
  open(my $fhl, ">$fl") || die "cannot write $fl\n";
  print $fhl join("\t", qw/chr beg end loc tloc qloc/)."\n";
  open(my $fht, ">$ft") || die "cannot write $ft\n";
  print $fht join("\t", qw/tid tbeg tend qid qbeg qend type tiloc tdloc
    qiloc qdloc/)."\n";

  my $dir = "$ENV{'misc3'}/$qry\_$tgt/23_blat";
  my $gta = Tabix->new(-data => "$dir/31.5/gax.gz");
  my $gqa = Tabix->new(-data => "$dir/41.5/gax.gz");
  my $gtr = Tabix->new(-data => "$dir/31.9/gax.gz");
  my $gqr = Tabix->new(-data => "$dir/41.9/gax.gz");

  my $gapt = Tabix->new(-data => "$ENV{'genome'}/$tgt/16.gap.bed.gz");
  my $gapq = Tabix->new(-data => "$ENV{'genome'}/$qry/16.gap.bed.gz");

  while( <$fhi> ) {
    chomp;
    next if /(^\#)|(^\s*$)/;
    my ($tid, $tbeg, $tend, $tsrd, $qid, $qbeg, $qend, $qsrd, $cid, $lev) 
      = split "\t";
    $lev == 1 || next;
    my $tlen = $tend - $tbeg - 1;
    my $qlen = $qend - $qbeg - 1;
    my $qstr = "$qid:$qbeg-$qend";
    my $tstr = "$tid:$tbeg-$tend";
   
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
          print $fhd join("\t", $tid, $tb, $te, "$tid:$tb-$te", $tstr, $qstr)."\n";
        } elsif($type eq "cnv") {
          print $fhl join("\t", $tid, $tb, $te, "$tid:$tb-$te", $tstr, $qstr)."\n";
        } else {
          $type eq "tlc" || die "unknow type: $type\n";
          print $fht join("\t", $tid, $tb, $te, $qi, $qb, $qe, 't', $tstr, '', '', $qstr)."\n";
        }
      }
    }
    if($flagq) {
      my $stats = idm_cat($qid, $qbeg+1, $qend-1, $lev, $gqr, $gqa, $gtr);
      for (@$stats) {
        my ($qb, $qe, $type, $ti, $tb, $te) = @$_;
        if($type eq "del") {
          print $fhs join("\t", $tid, $tbeg, $tbeg, "$qid:$qb-$qe", $tstr, $qstr)."\n";
        } elsif($type eq "cnv") {
          print $fhg join("\t", $tid, $tbeg, $tbeg, "$qid:$qb-$qe", $tstr, $qstr)."\n";
        } else {
          $type eq "tlc" || die "unknow type: $type\n";
          print $fht join("\t", $ti, $tb, $te, $qid, $qb, $qe, 'q', '', $tstr, $qstr, '')."\n";
        }
      }
    }
  }
  close $fhi;
  close $fhs;
  close $fhd;
  close $fhl;
  close $fhg;
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


__END__
