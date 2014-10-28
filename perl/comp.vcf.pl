#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  comp.vcf.pl - call variants in VCF forat

=head1 SYNOPSIS
  
  comp.vcf.pl [-help] [-qry query-genome] [-tgt target-genome]

  Options:
    -h (--help)   brief help message
    -q (--qry)    query genome (def: HM056)
    -t (--tgt)    target genome (def: HM101)

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
use Location;
use List::Util qw/min max sum/;

my ($qry, $tgt) = ('HM056', 'HM101');
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "qry|q=s" => \$qry,
  "tgt|t=s" => \$tgt,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$qry || !$tgt;

my $data = $ENV{'data'};
my $qry_fas = "$data/genome/$qry/11_genome.fas";
my $tgt_fas = "$data/genome/$tgt/11_genome.fas";

my $dir = "$data/misc3/$qry\_$tgt/23_blat/31.9";
chdir $dir or die "cannot chdir to $dir\n";

runCmd("snp2vcf.pl -i snp -o snp.vcf -s $qry");
runCmd("idm2vcf.pl -q $qry -t $tgt");
runCmd("vcf-concat snp.vcf idm.vcf | vcf-sort > vnt.1.vcf");
runCmd("vcf.fix.indel.pl -i vnt.1.vcf -o vnt.vcf");
runCmd("vcf2tbl.pl -i vnt.vcf -o vnt.tbl");
runCmd("rm vnt.*.vcf");


__END__
