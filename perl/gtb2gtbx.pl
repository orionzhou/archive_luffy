#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gtb2gtbx.pl - convert a Gtb file to Gtbx format

=head1 SYNOPSIS
  
  gtb2gtbx.pl [-help] [-ref refseq-fasta] [-in input-file] [-out output-file]

  Options:
      -help   brief help message
      -in     input file
      -out    output file
      -ref    reference sequence file

=head1 DESCRIPTION

  This program converts an input Gtb file to an output Gtbx file

=cut
  
#### END of POD documentation.
#-----------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;
use Gtb;

my ($fi, $fo, $fr) = ('') x 3;
my $help_flag;

#----------------------------------- MAIN -----------------------------------#
GetOptions(
    "help|h"   => \$help_flag,
    "in|i=s"   => \$fi,
    "out|o=s"  => \$fo,
    "ref|r=s"  => \$fr,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo || !$fr;

my ($fhi, $fho);
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

gtb2gtbx($fhi, $fho, $fr);



__END__
