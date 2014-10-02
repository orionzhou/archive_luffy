#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  sort.header.pl - sort a tabular file with header by coordinate

=head1 SYNOPSIS
  
  sort.header.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input (Tabular) file
    -o (--out)    output file

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

my ($fi, $fo) = ('') x 2;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "in|i=s"  => \$fi,
  "out|o=s" => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

runCmd("(head -n 1 $fi && tail -n +2 $fi | sort -k1,1 -k2,2n -k3,3n) > $fo");

exit 0;
