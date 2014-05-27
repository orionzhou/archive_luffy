#!/usr/bin/perl -w
use strict;
use FindBin;
use lib $FindBin::Bin;
use File::Path qw/make_path remove_tree/;
use Data::Dumper;
use List::Util qw/min max sum/;
use Common;
use Bam;

my $org = "HM340.AP";
my $dir = "/home/youngn/zhoup/Data/genome/$org/augustus";
-d $dir || make_path($dir);
chdir $dir || die "cannot chdir $dir\n";

my $bam_in = "\$misc2/rnaseq/mt/24_genome/$org.bam";
#runCmd("bamtools filter -isPrimaryAlignment 1 -isMapped 1 -isMateMapped 1 -in $bam_in -out 11.f.bam");
#runCmd("samtools view -H 11.f.bam > 12.header.txt");

#runCmd("samtools sort 11.f.bam 14.sf");
runCmd("bam2hints --intronsonly --in=14.sf.bam --out=15.hints.gff");

runCmd("augustus --species=medicago --extrinsicCfgFile=\$AUGUSTUS_CONFIG_PATH/extrinsic/extrinsic.M.RM.E.W.cfg --alternatives-from-evidence=true --hintsfile=15.hints.gff --allow_hinted_splicesites=atac --introns=on --genemodel=complete --strand=both --gff3=on --outfile=21.gff ../11_genome.fas");



