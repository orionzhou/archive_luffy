#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  comp.pl - pipeline of pairwise comparison between 2 genomes

=head1 SYNOPSIS
  
  comp.pl [-help] [-qry query-genome] [-tgt target-genome]

  Options:
    -h (--help)   brief help message
    -q (--qry)    query genome
    -t (--tgt)    target genome

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

my ($qry, $tgt) = ('HM056', 'HM101');
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "qry|q=s" => \$qry,
  "tgt|t=s" => \$tgt,
) or pod2usage(2);
pod2usage(1) if $help_flag;

my $data = '/home/youngn/zhoup/Data';
my $qry_fas = "$data/genome/$qry/11_genome.fas";
my $tgt_fas = "$data/genome/$tgt/11_genome.fas";
my $qry_2bit = "$data/db/blat/$qry.2bit";
my $tgt_2bit = "$data/db/blat/$tgt.2bit";
my $qry_size = "$data/genome/$qry/15.sizes";
my $tgt_size = "$data/genome/$tgt/15.sizes";
my $qry_size_bed = "$data/genome/$qry/15.bed";
my $tgt_size_bed = "$data/genome/$tgt/15.bed";
my $qry_gap = "$data/genome/$qry/16_gap.bed";
my $tgt_gap = "$data/genome/$tgt/16_gap.bed";

my $dir = "$data/misc3/$qry\_$tgt";

my $d23 = "$dir/23_blat";
-d $d23 || make_path($d23);
chdir $d23 || die "cannot chdir to $d23\n";

#process_blat1();
#process_blat2();
process_blat3();

sub process_blat1 { # blat 11_genome.fa -> 11.psl
  runCmd("psl2gal.pl -i 11.psl -o 11.gal");
  runCmd("galfix.pl -i 11.gal -o - | \\
    galrmgap.pl -i - -q $qry_gap -t $tgt_gap -o - | \\
    galfillstat.pl -i - -q $qry_fas -t $tgt_fas -o 12.fixed.gal");
  runCmd("gal2psl.pl -i 12.fixed.gal -o 12.fixed.psl");

  runCmd("axtChain -linearGap=medium -psl 12.fixed.psl \\
    $tgt_2bit $qry_2bit 21.chain");
  runCmd("chainPreNet 21.chain $tgt_size $qry_size 23.chain");
  runCmd("chain2gal.pl -i 23.chain -o - | \\
    galfillstat.pl -i - -q $qry_fas -t $tgt_fas -o 23.gal");
  runCmd("gal2gax.pl -i 23.gal -o 23.gax");
  runCmd("gax2bed.pl -i 23.gax -p qry -o - | sortBed -i stdin | \\
    mergeBed -i stdin > 23.bed");
  runCmd("subtractBed -a $qry_size_bed -b $qry_gap | \\
    subtractBed -a stdin -b 23.bed | \\
    awk '(\$3-\$2) >= 50' - > 24.nov.bed");
  runCmd("seqret.pl -d $qry_fas -b 24.nov.bed -o 24.nov.fa");
  runCmd("rm 11.gal 23.chain 23.gax");
}

sub pipe_chain_net {
  my ($pre, $qFas, $tFas, $q2bit, $t2bit, $qSize, $tSize) = @_;
  runCmd("axtChain -linearGap=medium -psl $pre.1.psl \\
    $t2bit $q2bit $pre.2.chain");
  runCmd("chainPreNet $pre.2.chain $tSize $qSize $pre.3.chain");
  runCmd("chain2gal.pl -i $pre.3.chain -o - | \\
    galfillstat.pl -i - -q $qFas -t $tFas -o $pre.3.gal");
  runCmd("chainSwap $pre.3.chain $pre.3.swap.chain");
  gal_expand("$pre.3.gal", $qFas, $tFas, $qSize, $tSize);

  runCmd("chainNet $pre.3.chain $tSize $qSize stdout /dev/null | \\
    netSyntenic stdin $pre.5.net");
  runCmd("netChainSubset $pre.5.net $pre.3.chain stdout | \\
    chainSort stdin $pre.5.chain");
#  runCmd("chainSwap $pre.5.chain stdout | \\
#    chainSort stdin $pre.5.swap.chain");
  runCmd("chain2gal.pl -i $pre.5.chain -o - | \\
    galfillstat.pl -i - -q $qFas -t $tFas -o $pre.5.gal");
  gal_expand("$pre.5.gal", $qFas, $tFas, $qSize, $tSize);

  runCmd("chainNet $pre.5.chain $tSize $qSize /dev/null stdout | \\
    netSyntenic stdin $pre.8.swap.net");
  runCmd("netChainSubset $pre.8.swap.net $pre.3.swap.chain \\
    $pre.8.swap.chain");
  runCmd("chainSwap $pre.8.swap.chain $pre.8.chain");
  runCmd("chain2gal.pl -i $pre.8.chain -o - | \\
    galfillstat.pl -i - -q $qFas -t $tFas -o $pre.8.gal");
  runCmd("galfilter.pl -i $pre.8.gal -m 100 -p 0.6 -o $pre.9.gal");
  gal_expand("$pre.9.gal", $qFas, $tFas, $qSize, $tSize);

  runCmd("rm $pre.3.chain $pre.3.swap.chain $pre.5.net $pre.5.chain \\
    $pre.8.swap.net $pre.8.swap.chain $pre.8.chain");
}
sub process_blat2 { # blat 24_nov.fa -> 24.nov.psl
  runCmd("psl2gal.pl -i 24.nov.psl -o - | \\
    galcoord.pl -i - -p qry -q $qry_size -o - | \\
    galfix.pl -i - -o 25.gal");
  runCmd("gal2psl.pl -i 25.gal -o 25.psl");
  runCmd("cat 12.fixed.psl 25.psl > 31.1.psl");
  runCmd("pslSwap 31.1.psl 41.1.psl");
  pipe_chain_net("31", $qry_fas, $tgt_fas, $qry_2bit, $tgt_2bit, $qry_size, $tgt_size);
  pipe_chain_net("41", $tgt_fas, $qry_fas, $tgt_2bit, $qry_2bit, $tgt_size, $qry_size);
  runCmd("rm 25.gal");
}

sub process_blat3 {
  gal_expand("31.3.gal", $qry_fas, $tgt_fas, $qry_size, $tgt_size);
  gal_expand("31.5.gal", $qry_fas, $tgt_fas, $qry_size, $tgt_size);
  gal_expand("31.9.gal", $qry_fas, $tgt_fas, $qry_size, $tgt_size);
  gal_expand("41.3.gal", $tgt_fas, $qry_fas, $tgt_size, $qry_size);
  gal_expand("41.5.gal", $tgt_fas, $qry_fas, $tgt_size, $qry_size);
  gal_expand("41.9.gal", $tgt_fas, $qry_fas, $tgt_size, $qry_size);
}
sub gal_expand {
  my ($fi, $qFas, $tFas, $qSize, $tSize) = @_;
  my $base = basename($fi, ".gal");
  runCmd("gal2bb.pl -i $fi -s $tSize -o $base.bb");
  runCmd("rm $base.gax*");
  
  -d $base || make_path($base);
  chdir $base || die "cannod chdir to $base";
  
  runCmd("idxgal.pl -i ../$base.gal -o gax");
  runCmd("gal2snp.pl -i ../$base.gal -o snp -q $qFas -t $tFas");
  runCmd("idxsnp.pl -i snp -s $tSize");

#  runCmd("gal2idm.pl -i ../$base.gal -o idm -q $qFas -t $tFas");
#  runCmd("idxidm.pl -i idm -s $tSize");
  chdir "..";
}

__END__

