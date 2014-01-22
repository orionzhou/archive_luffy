#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gtb2bigbed.pl - convert a Gtb file to BIGBED format

=head1 SYNOPSIS
  
  gtb2bigbed.pl [-help] [-in input-file] [-out output-file]

  Options:
      -help   brief help message
      -in     input file
      -out    output file
      -size   chromosome size file

=head1 DESCRIPTION

  This program converts an input Gtb file to an output BIGBED file

=cut
  
#### END of POD documentation.
#-----------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;

my ($fi, $fo, $fs) = ('') x 3;
my $help_flag;

#----------------------------------- MAIN -----------------------------------#
GetOptions(
    "help|h"   => \$help_flag,
    "in|i=s"   => \$fi,
    "out|o=s"  => \$fo,
    "size|s=s" => \$fs
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo || !$fs;

my $ft = "$fo.tmp";
runCmd("gtb2bed.pl -i $fi -o $ft", 1);
runCmd("bedSort $ft $ft", 1);
runCmd("bedToBigBed -tab $ft $fs $fo", 1);
runCmd("rm $ft", 1);



__END__
