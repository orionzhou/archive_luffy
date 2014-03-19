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
my $qry_fas = "$data/genome/$qry/11_genome.fa";
my $tgt_fas = "$data/genome/$tgt/11_genome.fa";
my $qry_2bit = "$data/db/blat/$qry.2bit";
my $tgt_2bit = "$data/db/blat/$tgt.2bit";
my $qry_size = "$data/genome/$qry/15.sizes";
my $tgt_size = "$data/genome/$tgt/15.sizes";
my $qry_gap = "$data/genome/$qry/16_gap.tbl";
my $tgt_gap = "$data/genome/$tgt/16_gap.tbl";

my $dir = "$data/misc3/$qry\_$tgt";

my $d23 = "$dir/23_blat";
-d $d23 || make_path($d23);
chdir $d23 || die "cannot chdir to $d23\n";

#pipe_blat();
#pipe_prepare_idx();

sub pipe_blat {
  runCmd("psl2gal.pl -i 11.psl -o 11.gal");
  runCmd("galfix.pl -i 11.gal -o - | \\
    galrmgap.pl -i - -q $qry_gap -t $tgt_gap -o - | \\
    galfillstat.pl -i - -q $qry_fas -t $tgt_fas -o 12.fixed.gal");
  runCmd("gal2psl.pl -i 12.fixed.gal -o 12.fixed.psl");
  runCmd("pslSwap 12.fixed.psl 14.swap.psl");

  runCmd("axtChain -linearGap=medium -psl 12.fixed.psl \\
    $tgt_2bit $qry_2bit 21.chain");
  runCmd("chainPreNet 21.chain $tgt_size $qry_size 23.chain");
  runCmd("chain2gal.pl -i 23.chain -o - | \\
    galfillstat.pl -i - -q $qry_fas -t $tgt_fas -o 23.gal");
  runCmd("chainNet 23.chain $tgt_size $qry_size stdout /dev/null | \\
    netSyntenic stdin 25.net");
  runCmd("netChainSubset 25.net 23.chain 25.chain");
  runCmd("chain2gal.pl -i 25.chain -o - | \\
    galfillstat.pl -i - -q $qry_fas -t $tgt_fas -o 25.gal");
  runCmd("galfilter.pl -i 25.gal -m 100 -p 0.6 -o 26.gal");

  runCmd("axtChain -linearGap=medium -psl 14.swap.psl \\
    $qry_2bit $tgt_2bit 31.chain");
  runCmd("chainPreNet 31.chain $tgt_size $tgt_size 33.chain");
  runCmd("chain2gal.pl -i 33.chain -o - | \\
    galfillstat.pl -i - -q $tgt_fas -t $qry_fas -o 33.gal");
  runCmd("chainNet 33.chain $qry_size $tgt_size stdout /dev/null | \\
    netSyntenic stdin 35.net");
  runCmd("netChainSubset 35.net 33.chain 35.chain");
  runCmd("chain2gal.pl -i 35.chain -o - | \\
    galfillstat.pl -i - -q $tgt_fas -t $qry_fas -o 35.gal");
  runCmd("galfilter.pl -i 35.gal -m 100 -p 0.6 -o 36.gal");
}

sub gal_callvnt {
  my ($fi, $qFas, $tFas, $qSize, $tSize) = @_;
  my $base = basename($fi, ".gal");
  runCmd("galexpand.pl -i $fi -o $base.gax");
  runCmd("sort -k1,1 -k2,2n -k3,3n $base.gax -o $base.gax");
  runCmd("bgzip -c $base.gax > $base.gax.gz");
  runCmd("tabix -s 1 -b 2 -e 3 $base.gax.gz");
  
  runCmd("gal2psl.pl -i $base.gal -o $base.psl");
  runCmd("pslToBed $base.psl $base.bed");
  runCmd("bedSort $base.bed $base.bed");
  runCmd("bedToBigBed -tab $base.bed $tSize $base.bb");
  runCmd("rm $base.bed $base.psl");

  -d "$base.vnt" || make_path("$base.vnt");
  chdir "$base.vnt" || die "cannod chdir to $base.vnt";
  runCmd("gal2snp.pl -i ../$base.gal -o snp -q $qFas -t $tFas");
  runCmd("sort -k1,1 -k2,2n snp -o snp");
  runCmd("bgzip -c snp > snp.gz");
  runCmd("tabix -s 1 -b 2 -e 2 snp.gz");

  runCmd("gal2indel.pl -i ../$base.gal -o idm");
  runCmd("sort -k1,1 -k2,2n -k3,3n idm -o idm");
  runCmd("bgzip -c idm > idm.gz");
  runCmd("tabix -s 1 -b 2 -e 3 idm.gz");
  chdir "..";
}
sub pipe_prepare_idx {
  runCmd("ln -sf 23.gal 51.gal");
  gal_callvnt("$d23/51.gal", $qry_fas, $tgt_fas, $qry_size, $tgt_size);

  runCmd("chain2gal.pl -i 25.chain -o - | \\
    galfillstat.pl -i - -q $qry_fas -t $tgt_fas -o 35.gal");
  runCmd("galfilter.pl -i 25.gal -m 100 -p 0.6 -o 26.gal");
  gal_callvnt("$d23/26.gal", $qry_fas, $tgt_fas, $qry_size, $tgt_size);

  runCmd("galswap.pl -i 51.gal -o 56.gal");
  gal_callvnt("$d23/56.gal", $tgt_fas, $qry_fas, $tgt_size, $qry_size);

  runCmd("chain2gal.pl -i 35.chain -o - | \\
    galfillstat.pl -i - -q $tgt_fas -t $qry_fas -o 35.gal");
  runCmd("galfilter.pl -i 35.gal -m 100 -p 0.6 -o 36.gal");
  gal_callvnt("$d23/36.gal", $tgt_fas, $qry_fas, $tgt_size, $qry_size);
}

__END__

