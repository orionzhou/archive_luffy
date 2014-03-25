#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  idxgal.pl - create index files for a GAL file

=head1 SYNOPSIS
  
  idxgal.pl [-help] [-i input-gal] [-s genome-size-file]

  Options:
      -h (--help)   brief help message
      -i (--in)     input file (Gal format)
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

my ($fi, $fs) = ('', '');
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "size|s=s" => \$fs,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fs;

my $base = basename($fi, ".gal");
runCmd("galexpand.pl -i $fi -o $base.gax");
runCmd("sort -k1,1 -k2,2n -k3,3n $base.gax -o $base.gax");
runCmd("bgzip -c $base.gax > $base.gax.gz");
runCmd("tabix -s 1 -b 2 -e 3 $base.gax.gz");

runCmd("gal2psl.pl -i $base.gal -o $base.psl");
runCmd("pslToBed $base.psl $base.bed");
runCmd("bedSort $base.bed $base.bed");
runCmd("bedToBigBed -tab $base.bed $fs $base.bb");
runCmd("rm $base.bed $base.psl");



