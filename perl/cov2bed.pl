#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  cov2bed.pl - covert a COVerage file to Bed format

=head1 SYNOPSIS
  
  cov2bed.pl [-help] [-in input-file] [-out output-file]

  Options:
      -help   brief help message
      -in     input
      -out    output

=cut
  
#### END of POD documentation.
#-----------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;

#----------------------------------- MAIN -----------------------------------#
my ($fi, $fo) = ('') x 2;
my ($fhi, $fho);
my $help_flag;
GetOptions(
    "help|h"  => \$help_flag,
    "in|i=s"  => \$fi,
    "out|o=s" => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;

if ($fi eq '' || $fi eq "stdin" || $fi eq "-") {
    $fhi = \*STDIN;
} else {
    open ($fhi, "<$fi") || die "cannot read $fi: $!\n";
}

if ($fo eq '' || $fo eq "stdout" || $fo eq "-") {
    $fho = \*STDOUT;
} else {
    open ($fho, ">$fo") || die "cannot write $fo: $!\n";
}

while(<$fhi>) {
    chomp;
    my @ps = split "\t";
    next unless @ps == 3;
    my ($chr, $pos, $cov) = @ps;
    print $fho join("\t", $chr, $pos-1, $pos, $cov)."\n";
}
close $fhi;
close $fho;

exit 0;
