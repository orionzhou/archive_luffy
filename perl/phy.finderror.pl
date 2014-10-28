#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  phy.finderror.pl - 

=head1 SYNOPSIS
  
  phy.finderror.pl [-help]

  Options:
    -h (--help)   brief help message

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

my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
) or pod2usage(2);
pod2usage(1) if $help_flag;

my $dir = "$ENV{'misc3'}/phy_finderror";
-d $dir || make_path($dir);
chdir $dir || die "cannot chdir to $dir\n";

my $fv = "$ENV{'misc3'}/hapmap/30_vnt/acc319.vcf.gz";
my $fc = "$ENV{'misc3'}/hapmap/30_vnt/ingroup.txt";
#runCmd("bcftools view -U -m2 -M2 -c2 -O z -v snps -S $fc -o 01.vcf.gz $fv chr5");

#runCmd("bcftools view -O z -o 05.hm018.dn.chr5.vcf.gz 04.hm018.dn.vcf.gz chr5");
#runCmd("bcftools merge 01.vcf.gz 05.hm018.dn.chr5.vcf.gz -o 11.vcf");

#vcf2bed("11.vcf", "11.bed");
my $fg = "$ENV{'misc3'}/HM018_HM101/23_blat/31.9/gax.bed";
#runCmd("coverageBed -a $fg -b 11.bed | cut -f1,3,4 > 15.hm018.bed");
#vcf_fill("11.vcf", "15.hm018.bed", "21.vcf");
#runCmd("bgzip -f 21.vcf");
#runCmd("tabix -p vcf 21.vcf.gz");

sub vcf2bed {
  my ($fi, $fb) = @_;
  open(my $fhi, "<$fi") or die "cannot read $fi\n";
  open(my $fhb, ">$fb") or die "cannot write $fb\n";
  while(<$fhi>) {
    chomp;
    /^\#/ && next;
    my @ps = split "\t";
    @ps >= 10 || die "not 10 fields: $_\n";
    my ($chr, $pos) = @ps[0..1];
    print $fhb join("\t", $chr, $pos-1, $pos)."\n";
  }
  close $fhi;
  close $fhb;
}
sub vcf_fill {
  my ($fi, $fb, $fo) = @_;
  open(my $fhi, "<$fi") or die "cannot read $fi\n";
  open(my $fho, ">$fo") or die "cannot write $fo\n";

  open(my $fhb, "<$fb") or die "cannot read $fb\n";
  while(<$fhi>) {
    chomp;
    if(/^\#/) {
      print $fho "$_\n";
      next;
    }
    my @ps = split "\t";
    @ps >= 10 || die "not 10 fileds: $_\n";
    my ($chr, $pos) = @ps[0,1];

    my $idx = $#ps; 
    my $gt = $ps[$idx];
    $gt =~ /^([\d\.])\/\1/ or die "unknown gt: $gt\n";
    my $line = readline($fhb);
    
    if($1 ne ".") {
      print $fho join("\t", @ps)."\n";
      next;
    }
    chomp $line;
#    print $orgs[$idx]." ".$line."\n";
    my ($chr2, $pos2, $gt2) = split("\t", $line);
    ($chr eq $chr2 && $pos == $pos2) || 
      die "error $idx $chr:$pos $chr2:$pos2\n";
    $ps[$idx] = "0/0" if $gt2;
    print $fho join("\t", @ps)."\n";
  }
  close $fhi;
  close $fho;
  close $fhb;;
}

#runCmd("bcftools view -m2 -M2 -O u -v snps 21.vcf.gz chr5 | \\
#  bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%SAMPLE=%GT]\\n' - | \\
#  snpfilter.pl -n 80 -m 1 -o 51.snp");

runCmd("samplelines.pl -i 51.snp -o 52.snp -n 10000");
runCmd("snp2phy.pl -i 52.snp -o 52.phy -r HM101_ref");
runCmd("phyml -i 52.phy -d nt");
runCmd("mv 52.phy_phyml_tree.txt 52.nwk");
runCmd("rm 52.phy_phyml*");




