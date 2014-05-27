#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  idxgal.pl - create index files for a GAL file

=head1 SYNOPSIS
  
  idxgal.pl [-help] [-i input-file] [-o output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (Gal format)
    -o (--out)    output prefix

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

my ($fi, $fo) = ('', '');
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

runCmd("gal2gax.pl -i $fi -o $fo");
runCmd("sort -k1,1 -k2,2n -k3,3n $fo -o $fo");
runCmd("bgzip -c $fo > $fo.gz");
runCmd("tabix -s 1 -b 2 -e 3 $fo.gz");

