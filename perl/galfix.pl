#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  galfix.pl - Check and fix a GAL file

=head1 SYNOPSIS
  
  galfix.pl [-help] [-in input-file] [-out output-file]

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

my ($fi, $fo) = ('') x 2;
my ($fq, $ft) = ('') x 2; 
my ($fhi, $fho);
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
  "qry|q=s"  => \$fq,
  "tgt|t=s"  => \$ft,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

if ($fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "Can't open file $fi: $!\n";
}

if ($fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

print $fho join("\t", @HEAD_GAL)."\n";

my $cnt = 0;
while( <$fhi> ) {
  chomp;
  next if /(^id)|(^\#)|(^\s*$)/;
  my $ps = [ split "\t" ];
  next unless @$ps == 20;
  my ($id, $tId, $tBeg, $tEnd, $tSrd, $tSize, 
    $qId, $qBeg, $qEnd, $qSrd, $qSize,
    $ali, $mat, $mis, $qN, $tN, $ident, $score, $tLocS, $qLocS) = @$ps;
  $tSrd eq "+" || die "$id: tSrd -\n";
  
  my ($qLoc, $tLoc) = (locStr2Ary($qLocS), locStr2Ary($tLocS));
  my @qlens = map {$_->[1]-$_->[0]+1} @$qLoc;
  my $ref = tiling($qLoc, \@qlens, 2);
  my (@rqloc, @rtloc);
  for (@$ref) {
    my ($rqb, $rqe, $idx) = @$_;

    my ($qb, $qe) = @{$qLoc->[$idx]};
    my ($tb, $te) = @{$tLoc->[$idx]};
    my $rtb = $rqb - $qb + $tb;
    my $rte = $rqe - $qb + $tb;

    if(@rqloc == 0 || $rtb > $rtloc[-1]->[1]) { 
        push @rqloc, [$rqb, $rqe];
        push @rtloc, [$rtb, $rte];
    }
  }
  my ($nqLocS, $ntLocS) = (locAry2Str(\@rqloc), locAry2Str(\@rtloc));
  if($nqLocS ne $qLocS) {
    @$ps[18,19] = ($ntLocS, $nqLocS);
    $cnt ++;
  }
  print $fho join("\t", @$ps)."\n";
}
print STDERR "$cnt rows fixed\n";
close $fhi;
close $fho;


__END__
