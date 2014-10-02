#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  comp.ortho.pl - find orthologs in pairwise comparison

=head1 SYNOPSIS
  
  comp.ortho.pl [-help] [-qry qry-genome] [-tgt tgt-genome]

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
use List::Util qw/min max sum/;

my ($qry, $tgt) = qw/HM034 HM101/;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "qry|q=s"  => \$qry,
  "tgt|t=s"  => \$tgt,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$qry || !$tgt;

my $qdir = "$ENV{'genome'}/$qry";
my $tdir = "$ENV{'genome'}/$tgt";
my $cdir = "$ENV{'misc3'}/$qry\_$tgt/23_blat";

my $gax = Tabix->new(-data => "$cdir/31.9/gax.gz");
my $snp = Tabix->new(-data => "$cdir/31.9/snp.gz");

my $fgq = "$qdir/51.gtb";
my $fgt = "$tdir/51.gtb";
my $qdb = Bio::DB::Fasta->new("$qdir/51.fas");
my $tdb = Bio::DB::Fasta->new("$tdir/51.fas");

my $qgene = Tabix->new(-data => "$qdir/51.tbl.gz");

my $fo = "$cdir/ortho.tbl";
open(my $fho, ">$fo") or die "cannot write $fo\n";
print $fho join("\t", qw/tid tlen qid qlen slen tseq qseq/)."\n";

my $tq = readTable(-in => $fgq, -header => 1);
my $hq;
for my $i (0..$tq->lastRow) {
  my ($qid, $qpar, $qchr, $qb, $qe, $qsrd, 
    $elocS, $ilocS, $clocS, $flocS, $tlocS, $phase) = $tq->row($i);
  my $qloc = locStr2Ary($clocS);
  $hq->{$qid} = [0, locAryLen($qloc)];
}

my $tt = readTable(-in => $fgt, -header => 1);
for my $i (0..$tt->lastRow) {
  my ($tid, $tpar, $tchr, $tb, $te, $tsrd, 
    $elocS, $ilocS, $clocS, $flocS, $tlocS, $phase, 
    $src, $conf, $cat1, $cat2, $cat3, $note) = $tt->row($i);
  my $tloc = locStr2Ary($clocS); 
  $tloc = $tsrd eq "+" ? [ map {[$tb+$_->[0]-1, $tb+$_->[1]-1]} @$tloc ]
    : [ map {[$te-$_->[1]+1, $te-$_->[0]+1]} @$tloc ];
  my $tlen = locAryLen($tloc);
  
  my $h;
  for (@$tloc) {
    my ($tbeg, $tend) = @$_;
    my $ary = read_gax($gax, $tchr, $tbeg, $tend);
    for (@$ary) {
      my ($cid, $tid1, $tb1, $te1, $tsrd1, $qid, $qb, $qe, $qsrd) = @$_;
      my $ary2 = read_cds($qgene, $qid, $qb, $qe);
      for (@$ary2) {
        my ($chr, $beg, $end, $srd, $gene) = @$_;
        $h->{$gene} ||= 0;
        $h->{$gene} += $end - $beg + 1;
      }
    }
  }
  if(!defined($h)) {
    print $fho join("\t", $tid, $tlen, ("") x 5)."\n";
    next;
  }
  my @qids = sort {$h->{$b} <=> $h->{$a}} keys(%$h);
  my $qid = $qids[0];
  my $qlen = $hq->{$qid}->[1];
  my $slen = $h->{$qid};
  if($slen / $tlen >= 0.2 && $slen / $qlen >= 0.2) {
    $hq->{$qid}->[0] = 1;
    my $tseq = $tdb->seq($tid);
    my $qseq = $qdb->seq($qid);
    $tseq =~ s/X$//i;
    $qseq =~ s/X$//i;
    print $fho join("\t", $tid, $tlen, $qid, $qlen, $slen, $tseq, $qseq)."\n";
  } else {
    print $fho join("\t", $tid, $tlen, ("") x 5)."\n";
  }
}

for my $qid (sort(keys(%$hq))) {
  $hq->{$qid}->[0] == 0 || next;
  my $qlen = $hq->{$qid}->[1];
  print $fho join("\t", ("") x 2, $qid, $qlen, ('') x 3)."\n";
}
close $fho;

sub needle_ident {
  my ($seq1, $seq2) = @_;
  runCmd("echo -e \">seq1\\n$seq1\" > seq1.fas", 0); 
  runCmd("echo -e \">seq2\\n$seq2\" > seq2.fas", 0); 
  runCmd("needle seq1.fas seq2.fas -gapopen 10 -gapextend 0.5 \\
    -aformat pair -outfile aln.needle", 0);
  open(my $fha, "<aln.needle") or die "cannot read aln.needle\n";
  my ($len, $idt, $sim, $gap, $sco);
  while(<$fha>) {
    chomp;
    $len = $1 if /^\# Length:\s*(\d+)/;
    $idt = $1 if /^\# Identity:\s*(\d+)/;
    $sim = $1 if /^\# Similarity:\s*(\d+)/;
    $gap = $1 if /^\# Gaps:\s*(\d+)/;
    $sco = $1 if /^\# Score:\s*([\d\.]+)/;
  }
  close $fha;
  runCmd("rm seq1.fas seq2.fas aln.needle", 0);
  my $ident = sprintf "%.03f", $idt / $len;
  return $ident;
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


__END__

