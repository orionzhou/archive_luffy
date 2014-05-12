#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gal2bb.pl - convert Gal to BigBed file

=head1 SYNOPSIS
  
  gal2bb.pl [-help] [-i input-file] [-o output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (Gal format)
    -o (--out)    output file (BigBed format)
    -s (--size)   chrom-size file for target genome

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin";
use Common;
use Data::Dumper;
use File::Path qw/make_path remove_tree/;
use File::Basename;
use List::Util qw/min max sum/;

my ($fi, $fs, $fo) = ('', '', '');
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
  "size|s=s" => \$fs,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fs;

my $base = basename($fi, ".gal");
$fo ||= "$base.bb";

runCmd("gal2psl.pl -i $base.gal -o $base.psl");
runCmd("pslToBed $base.psl stdout | \\
  fixbedbygal.pl -i - -g $base.gal -o $base.bed");
runCmd("bedSort $base.bed $base.bed");
runCmd("bedToBigBed -tab $base.bed $fs $fo");
runCmd("rm $base.bed $base.psl");



