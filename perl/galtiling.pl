#!/usr/bin/perl -w
use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Tabix;
use Common;
use Gal;
use Location;
use Data::Dumper;
use List::Util qw/min max sum/;

my ($fi, $fo) = ('') x 2;
my ($chr, $beg, $end) = ('') x 3;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "chr|c=s"   => \$chr,
  "beg|b=i"   => \$beg,
  "end|e=i"   => \$end,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$chr || !$beg || !$end;

my $QRY = "HM056";
my $dir = "/home/youngn/zhoup/Data/misc3/$QRY\_HM101/23_blat";

my $fgal = "$dir/23.gal";
my $fgax = "$dir/23.gax.gz";
my $gax = Tabix->new(-data=>$fgax);

open(my $fh, "<$fgal") or die "cannot read $fgal\n";
my $hs;
while(<$fh>) {
  chomp;
  next if /^(\#)|(id)/;
  my @ps = split "\t";
  my ($id, $tId, $tBeg, $tEnd, $tSrd, $tSize, 
    $qId, $qBeg, $qEnd, $qSrd, $qSize,
    $ali, $mat, $mis, $qN, $tN, $ident, $score, $tLocS, $qLocS) = @ps;
  ! exists $hs->{$id} || die "$id twice\n";
  $hs->{$id} = [$qId, int($score)];
}

my $ary = read_gax($gax, $chr, $beg, $end, "t");

my @locs;
my @scores;
my ($hc, $hi);
my $idx = 0;
for (@$ary) {
  my ($id, $tid, $tb, $te, $tsrd, $qid, $qb, $qe, $qsrd) = @$_;
  my $score = $hs->{$id}->[1];
  
  if(!exists $hc->{$id}) {
    $hc->{$id} = $idx;
    $hi->{$idx} = $id;
    $idx ++;
  }
  print join("\t", $tb, $te, $qid, $score)."\n"; 
  push @locs, [$tb, $te, $hc->{$id}];
  push @scores, $score;
}
my $ref = tiling(\@locs, \@scores, 2);
for (@$ref) {
  my ($beg, $end, $idx, $idxs) = @$_;
  my $id = $hi->{$idx};
  print join("\t", $beg, $end, @{$hs->{$id}}, $end-$beg+1)."\n";
}
