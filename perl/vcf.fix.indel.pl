#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  vcf.fix.indel.pl - LeftAlignAndTrim InDels

=head1 SYNOPSIS
  
  vcf.fix.indel.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input VCF file
    -o (--out)    output file
    -r (--ref)    ref-fasta file (default: $genome/HM101/11_genome.fasta)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use Common;
use Gatk;
use List::Util qw/min max sum/;

my ($fi, $fo) = ('') x 2;
my $fr = "\$genome/HM101/11_genome.fasta";
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

runCmd("java -jar $gatk -T LeftAlignAndTrimVariants \\
  -R $fr --trimAlleles \\
  --variant $fi -o $fo");
runCmd("rm $fi*.idx $fo*.idx");



__END__
