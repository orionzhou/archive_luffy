#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gal2chain.pl - convert GAL file to Chain format

=head1 SYNOPSIS
  
  gal2chain.pl [-help] [-in input-file] [-out output-file]

  Options:
    -help   brief help message
    -in     input file
    -out    output file

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
use Gal;
use List::Util qw/min max sum/;

my ($fi, $fo) = ('', '');
my ($fhi, $fho);
my $help_flag;

#----------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;

if ($fi eq '' || $fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "Can't open file $fi: $!\n";
}

if ($fo eq '' || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

while( <$fhi> ) {
  chomp;
  next if /(^id)|(^\#)|(^\s*$)/;
  my $ps = [ split "\t" ];
  next unless @$ps == 20;
  my ($id, $tId, $tBeg, $tEnd, $tSrd, $tSize, 
    $qId, $qBeg, $qEnd, $qSrd, $qSize,
    $ali, $mat, $mis, $qN, $tN, $ident, $score, $tLocS, $qLocS) = @$ps;
  $tSrd eq "+" || die "$id: tSrd -\n";

  my ($qloc, $tloc) = (locStr2Ary($qlocS), locStr2Ary($tlocS));
  @$qloc == @$tloc || die "unequal pieces\n";
  my $nBlock = @$qloc;
  
  my ($ctb, $cte) = ($tb - 1, $te);
  my ($cqb, $cqe) = $qSrd eq "+" ? ($qb - 1, $qe) : 
    ($qsize-$qe+1 - 1, $qsize-$qb+1);
  print $fho join(" ", "chain", $score, $tId, $tSize, $tSrd, $ctb, $cte,
    $qId, $qSize, $qSrd, $cqb, $cqe, $id)."\n";
  if($nBlock > 1) {
    for my $i (0..$nBlock-2) {
      my $len = $tloc->[$i]->[1] - $tloc->[$i]->[0] + 1;
      my $dt = $tloc->[$i+1]->[0] - $tloc->[$i]->[1] - 1;
      my $dq = $qloc->[$i+1]->[0] - $qloc->[$i]->[1] - 1;
      print $fho join("\t", $len, $dt, $dq)."\n";
    }
  }
  my $tlocl = $tloc->[$nBlock-1];
  my $lenl = $tlocl->[1] - $tlocl->[0] + 1;
  print $fho $lenl."\n\n";
}
close $fhi;
close $fho;

__END__
