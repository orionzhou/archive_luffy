#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  toGal.pl - convert an input file to GAL format

=head1 SYNOPSIS
  
  toGal.pl [-help] [-in input-file] [-fmt input-format] [-out output-file]

  Options:
      -help   brief help message
      -in     input file
      -out    output file
      -fmt    input format (default: PSL)
      -qry    query-seq file 
      -tgt    target-seq file

=head1 DESCRIPTION

  This program converts an input file of specificed format to an output GAL file

=cut
  
#### END of POD documentation.
#-----------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;
use Gal;

my ($fi, $fo, $fmt) = ('', '', 'psl');
my ($fhi, $fho);
my ($fq, $ft) = ('') x 2; 
my $help_flag;

#----------------------------------- MAIN -----------------------------------#
GetOptions(
    "help|h"   => \$help_flag,
    "in|i=s"   => \$fi,
    "out|o=s"  => \$fo,
    "fmt|f=s"  => \$fmt,
    "qry|q=s"  => \$fq,
    "tgt|t=s"  => \$ft,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

if ($fi eq "stdin") {
    $fhi = \*STDIN;
} else {
    open ($fhi, $fi) || die "Can't open file $fi: $!\n";
}

if ($fo eq "stdout") {
    $fho = \*STDOUT;
} else {
    open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

if($fmt =~ /^psl$/i) {
    pod2usage(2) if !$fq || !$ft;
    psl2Gal($fhi, $fho, $fq, $ft);
} elsif($fmt =~ /^chain$/i) {
    pod2usage(2) if !$fq || !$ft;
    chain2Gal($fhi, $fho, $fq, $ft);
} elsif($fmt =~ /^net$/i) {
    net2Gal($fhi, $fho);
} elsif($fmt =~ /^blast$/i) {
    blast2Gal($fhi, $fho);
} else {
    die "unknown input format: $fmt\n";
}


__END__
