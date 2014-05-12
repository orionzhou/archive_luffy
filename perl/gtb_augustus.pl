#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gtb_augustus.pl - fix Gtb files output by Augustus

=head1 SYNOPSIS
  
  gtb_augustus.pl [-help] [-in input-file] [-out output-file]

  Options:
      -h (--help)   brief help message
      -i (--in)     input Gtb
      -o (--out)    output Gtb

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;
use Location;
use Gtb;

my ($fi, $fo) = ('') x 2;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

my ($fhi, $fho);
open ($fhi, $fi) || die "Can't open file $fi: $!\n";

if ($fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

my $t = readTable(-inh=>$fhi, -header=>1);
close $fhi;

for my $i (0..$t->nofRow-1) {
  my ($id, $par, $chr, $beg, $end, $srd, 
    $locES, $locIS, $locCS, $loc5S, $loc3S, $phase, 
    $src, $conf, $cat1, $cat2, $cat3, $note) = $t->row($i);
  next if $cat2 ne "mRNA";
  die "no CDS-loc for $id\n" unless $locCS;
  
  my ($rb, $re) = (1, $end - $beg + 1);
  my $locC = locStr2Ary($locCS);
  $locC = [ sort {$a->[0] <=> $b->[0]} @$locC ];
  if($locC->[0]->[0] > $rb) {
    my $los = $locC->[0]->[0] - $rb;
    $beg += $los if $srd eq "+";
    $end -= $los if $srd eq "-";
    ($rb, $re) = (1, $end - $beg + 1);
    $locC = [ map {[$_->[0] - $los, $_->[1] - $los]} @$locC ];
    my ($locI) = posDiff([[$rb, $re]], $locC);

    $t->setElm($i, "beg", $beg);
    $t->setElm($i, "end", $end);
    $t->setElm($i, "locE", locAry2Str($locC));
    $t->setElm($i, "locC", locAry2Str($locC));
    $t->setElm($i, "locI", locAry2Str($locI));
  }
  if($locC->[-1]->[1] < $re) {
    my $ros = $re - $locC->[-1]->[1];
    $end -= $ros if $srd eq "+";
    $beg += $ros if $srd eq "-";
    ($rb, $re) = (1, $end - $beg + 1);
    my ($locI) = posDiff([[$rb, $re]], $locC);
    
    $t->setElm($i, "beg", $beg);
    $t->setElm($i, "end", $end);
    $t->setElm($i, "locI", locAry2Str($locI));
  }
}
print $fho $t->tsv(1);
close $fho;



__END__
