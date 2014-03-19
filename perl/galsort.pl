#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  galsort.pl - sort a GAL file by tId + tBeg + tEnd fields

=head1 SYNOPSIS
  
  galsort.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file
    -o (--out)    output file (optional, default: --input)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;

my ($fi, $fo) = ('', '');
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi;

-s $fi || die "$fi not exist!\n";

if(!$fo || !(-s $fo)) {
  $fo = $fi;
}

runCmd("sort -k2,2 -k3,3n -k4,4n $fi -o $fo", 1);

__END__
