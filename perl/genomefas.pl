#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  genomefas.pl - process a genome fasta file

=head1 SYNOPSIS
  
  genomefas.pl [-help] [-in input-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input (fasta) file

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

my ($fi) = ('');
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "in|i=s"  => \$fi,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi;

my $dir = dirname($fi);
chdir $dir || die "cannot chdir to $dir\n";

runCmd("seqcheck.pl -i $fi -o 11_genome.fas");

runCmd("seqlen.pl -i 11_genome.fas -o 15.sizes");
runCmd("awk 'BEGIN {FS=\"\\t\"; OFS=\"\\t\"} {print \$1, 0, \$2}' 15.sizes > 15.bed");

runCmd("seqgap.pl -i 11_genome.fa -o 16_gap.bed -m 10");
#runCmd("awk 'BEGIN {FS=\"\\t\"; OFS=\"\\t\"} {if(NR != 1) \
#  {print \$1, \$2-1, \$3}}' 16_gap.tbl > 16_gap.bed");
runCmd("bedToBigBed 16_gap.bed 15.sizes 16_gap.bb");


__END__

