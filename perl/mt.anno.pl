#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  mt.anno.pl - annotate an Medicago genome

=head1 SYNOPSIS
  
  mt.anno.pl [-help] [-org organism]

  Options:
    -h (--help)   brief help message
    -g (--org)    genmome ID (organism) to process

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin";
use Common;
use Location;
use Data::Dumper;
use File::Path qw/make_path remove_tree/;
use File::Basename;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;

my ($org, $ncpu) = ('', 16);
my $help_flag;
#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "org|g=s"  => \$org,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$org;

my $dir = "$ENV{'genome'}/$org";
-d $dir || make_path($dir);
chdir $dir || die "cannot chdir $dir\n";

my $fg = "11_genome.fas";
-s $fg || die "$fg is not there\n";

my $f_aug = "augustus/41.gtb";
-s $f_aug or die "$f_aug is not there\n";
runCmd("ln -sf $f_aug 41.gtb");

my $f_nbs = "42.nbs/11.gtb";
-s $f_nbs or die "$f_nbs is not there\n";
runCmd("awk 'BEGIN {FS=\"\\t\"; OFS=\"\\t\"} {if(NR>1) \\
  {\$15=\"mRNA\"; \$16=\"NBS-LRR\"} print}' $f_nbs | \\
  cut -f1-18 > 42.nbs.gtb");

my $f_crp = "$ENV{'misc4'}/spada.crp.$org/61_final.gtb";
-s $f_crp or die "$f_crp is not there\n";
runCmd("awk 'BEGIN {FS=\"\\t\"; OFS=\"\\t\"} {if(NR>1) \\
  {\$15=\"mRNA\"; \$16=\"CRP\"; \$18=\"\"} print}' $f_crp | \\
  cut -f1-18 > 43.crp.gtb");

#runCmd("gtb.merge.pl -a 41.gtb -b 42.nbs.gtb -o 49.1.gtb");
#runCmd("gtb.merge.pl -a 49.1.gtb -b 43.crp.gtb -o 49.gtb");
#runCmd("gtb.dedup.pl -i 49.gtb -o 50.1.dedup.gtb");
#runCmd("gtb.pickalt.pl -i 50.1.dedup.gtb -o 50.2.pickalt.gtb");

my $f_ctm = "\$misc3/$org\_HM101/41_novseq/15.foreign.scf.txt";
if($org eq "HM101") {
  runCmd("gtb.fill.len.pl -i 50.2.pickalt.gtb | gtb.filter.pl -l 30 -o 51.gtb");
} else {
  runCmd("gtb.fill.len.pl -i 50.2.pickalt.gtb | gtb.filter.pl -l 30 -c $f_ctm -o 51.gtb");
}

runCmd("awk 'BEGIN {FS=\"\\t\"; OFS=\"\\t\"} \\
  {if(NR==1 || tolower(\$16) != \"te\") print}' 51.gtb > 55_noTE.gtb");

runCmd("gtb2gff.pl -i 51.gtb -o 51.gff");
runCmd("gtb.idx.pl -i 51.gtb -s 15.sizes");
runCmd("gtb2tbl.pl -i 51.gtb -o 51.tbl");
runCmd("gtb2fas.pl -i 51.gtb -d 11_genome.fas -o 51.fas");


__END__

