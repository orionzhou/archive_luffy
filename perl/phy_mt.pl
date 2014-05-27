#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  phy_mt.pl - extract SNPs and construct phylogeny

=head1 SYNOPSIS
  
  phy_mt.pl [-help] [-opt option]

  Options:
    -h (--help)   brief help message
    -t (--opt)    option

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

my ($opt) = ('ril2');
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "opt|t=s" => \$opt,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$opt;

my $dir = "/home/youngn/zhoup/Data/misc3/hapmap/31_phylogeny/$opt";
-d $dir || make_path($dir);
chdir $dir || die "cannot chdir to $dir\n";
my $fv = "../../30_vnt/acc319.vcf.gz";
-s $fv || die "$fv not there\n";
my $fc = "../../30_vnt/$opt.txt";
-s $fc || die "$fc not there\n";

my $d01 = "01_snp";
my $d11 = "11_snp_sample";
my $d12 = "12_phy";
my $d13 = "13_aln";
my $d21 = "21_phynj";
my $d22 = "22_phyml";
make_path($d01, $d11, $d12, $d13, $d21, $d22);

#run_clustalw_nj($d16, $d21, $chrs); 
#run_phyml($d15, $d22, $chrs);

my $regions = ["chr5"];
for my $reg (@$regions) {
  my $f01 = "$d01/$reg.snp";
  runCmd("bcftools view -U -m2 -M2 -c2 -O u -v snps -s $fc $fv $reg | 
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%SAMPLE=%GT]\\n' - |
    snpfilter.pl -n 8 -m 1 -o $f01");

  my $f11 = "$d11/$reg.snp";
  runCmd("samplelines.pl -i $f01 -o $f11 -n 10000");
  my $f12 = "$d12/$reg.phy";
  runCmd("snp2phy.pl -i $f11 -o $f12");
  
  my $f13 = "$d13/$reg.aln";
#    runCmd("seqret $f12 $f13 -osformat aln");
#    runCmd("clustalw2 -INFILE=$f13 -BOOTSTRAP=1000 -OUTORDER=INPUT -OUTPUTTREE=phylip -BOOTLABELS=node -CLUSTERING=NJ -KIMURA");
#    runCmd("mv $d13/$reg.phb $d21");
  
  my $f22 = "$d22/$reg.nwk";
  runCmd("/home/youngn/zhoup/Source/PhyML-3.1/phyml -i $f12 -d nt");
  runCmd("mv $f12\_phyml_tree.txt $f22");
  runCmd("rm $f12\_phyml*");
}

sub snp_convert {
    my ($dirI, $dirO, $chrs, $format) = @_;
    system("mkdir -p $dirO") unless -d $dirO;
    for my $chr (@$chrs) {
#    next if $chr eq "chr1";
        my $fo = "$dirO/$chr.$format";
        if($format =~ /^(aln)|(nexus)$/) {
            my $fi = "$dirI/$chr.phylip";
            runCmd("seqret $fi $fo -osformat $format");
        } elsif($format =~/^(phylip)|(structure)$/) {
            my $fi = "$dirI/$chr.ssp";
            runCmd("sspConvert -i $fi -o $fo -f $format");
        }
    }
}
sub run_phyml {
  my ($dirI, $dirO, $chrs) = @_;
  system("mkdir -p $dirO") unless -d $dirO;
  for my $chr (@$chrs) {
#    next unless $chr ge "chr7";
  }
}



