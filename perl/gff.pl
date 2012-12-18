#!/usr/bin/perl
use strict;
use Init;
use Common;
use Data::Dumper; 
use Seq;
use Gff;
use Gtb;
use Bed;
sub seqPreProcess {
    my ($dir, $db) = @_;
    if( $db eq "mt_30" ) {
        my $f01 = "$dir/01_all_bac.fa";
        my $f02 = "$dir/02_chr0-8.fa";
        my $f06 = "$dir/06_chr1-8.tbl";

        my $f09 = "$dir/09_bac_info.txt";
#    writeSeqDesc($f01, $f09);

        my $f03 = "$dir/03_chr0.fa";
        my $f11 = "$dir/11_chr1-8.fa";
#    splitBacFromChr($f02, $f11, $f03);
        my $f21 = "$dir/21_chr1-8.tbl";
#    pgp2Tbl($f06, $f21);
        
        my $f12 = "$dir/12_chr0.fa";
        my $f22 = "$dir/22_chr0.tbl";
#    mergeShortSeqs(-in=>$f03, -seqid=>"chr0", -gap=>50000, -type=>'BAC', -out1=>$f12, -out2=>$f22);
    } elsif($db eq "mt_35") {
        my $f00 = "$dir/00_all_bac.fa";
        my $f01 = "$dir/01_chr1-8.fa";
        my $f02 = "$dir/02_chr0.fa";
        my $f03 = "$dir/03_chrU.fa";
        my $f04 = "$dir/04_chrT.fa";
        my $f06 = "$dir/06_chr1-8.tbl";

        my $f11 = "$dir/11_chr1-8.fa";
#    changeChrId($f01, $f11);
        my $f21 = "$dir/21_chr1-8.tbl";
#    makeAsblTbl($f06, $f21);

        my $f12 = "$dir/12_chr0.fa";
        my $f22 = "$dir/22_chr0.tbl";
#    mergeShortSeqs(-in=>$f02, -seqid=>"chr0", -gap=>50000, -out1=>$f12, -out2=>$f22, -type=>'BAC');
        my $f13 = "$dir/13_chrU.fa";
        my $f23 = "$dir/23_chrU.tbl";
#    mergeShortSeqs(-in=>$f03, -seqid=>"chrU", -gap=>1000, -out1=>$f13, -out2=>$f23, -type=>"contig";
        my $f14 = "$dir/14_chrT.fa";
        my $f24 = "$dir/24_chrT.tbl";
#    mergeShortSeqs(-in=>$f04, -seqid=>"chrT", -gap=>1000, -out1=>$f14, -out2=>$f24, -type=>"tc";
    }
}
sub seqPostProcess {
    my ($dir, $db) = @_;
    my $f11 = "$dir/11_chr1-8.fa";
    my $f12 = "$dir/12_chr0.fa";
    my $f13 = "$dir/13_chrU.fa";
    my $f14 = "$dir/14_chrT.fa";
#  cat $f11-$f14 > 19.fa;
    my $f19 = "$dir/19.fa";
    
    my $f21 = "$dir/21_chr1-8.tbl";
    my $f22 = "$dir/22_chr0.tbl";
    my $f23 = "$dir/23_chrU.tbl";
    my $f24 = "$dir/24_chrT.tbl";
    my $f27 = "$dir/27_cen.tbl";
    
#  cat 2* > 29.tbl;
    my $f29 = "$dir/29.tbl";
    my $f31 = "$dir/31.gff";
    my $f32 = "$dir/32.bed";
    makeAsblGff($f29, $f31, $f32);
}
sub geneProcess {
    my ($dir, $fm, $f_seq) = rearrange(['dir', 'fm', 'f_seq'], @_);
    my $f01 = "$dir/01_gene.gff";
    my $f02 = "$dir/02_te.gff";
    my $f03 = "$dir/03_tRNA.gff";
#  cat $f01 $f02 $f03 > $f11 and remove ##gff headers;
    my $f11 = "$dir/11.gff";
    my $f12 = "$dir/12.gff";
#  format_gff_jcvi($f11, $f12);

#  cat $f12 $f23 $f33 > $f51 and remove ##gff headers;
    my $f51 = "$dir/51.gff";
    my $f52 = "$dir/52.gff";
#  format_gff_loc(-fi=>$f51, -fo=>$f52, -fm=>$fm);

    my $f61 = "$dir/61.gtb";
#  gff2Gtb($f52, $f61);
    my $f62 = "$dir/62_phase_fixed.gtb";
#  gtb_fix_phase($f61, $f62, $f_seq);
    my $f69 = $f62;

    my $f71 = "$dir/71.gff";
#  gtb2Gff($f69, $f71);
    my $f81 = "$dir/81.bed";
#  gtb2Bed($f69, $f81);
    my $f82 = "$dir/82_refined.bed";
#  bed_refine($f81, $f82);
    my $f86 = "$dir/86.tbl";
#  gtb2Tbl($f69, $f86);
    my $f91 = "$dir/91_proteome.fa";
#  gtb2Seq(-in=>$f69, -out=>$f91, -seq=>$f_seq, -opt=>2);
}
sub processCpMt {
    my ($dir) = @_;
    my $f01 = "$dir/01_cp.gb";
    my $f02 = "$dir/02.fa";
    my $f03 = "$dir/03.gff";
#  gb2Gff(-in=>$f01, -out=>$f03, -outseq=>$f02, -type=>'chloroplase_sequence');
    my $f11 = "$dir/11_mt.gb";
    my $f12 = "$dir/12.fa";
    my $f13 = "$dir/13.gff";
#  gb2Gff(-in=>$f11, -out=>$f13, -outseq=>$f12, -type=>'mitochondrial_sequence');
    my $f41 = "$dir/41.fa";
    my $f51 = "$dir/51.gff";
}
sub merge_gene_genome {
    my ($f_seq, $f_gff, $fo) = @_;
    my $firstline = `head -n 1 $f_seq`;
    if($firstline eq "##FASTA\n") {
        system("cat $f_gff $f_seq > $fo");
    } else {
        system("echo '##FASTA' | cat $f_gff - $f_seq > $fo");
    }
    runCmd('sed -i \'s/\\t\w*_gene\\t/\\tgene\\t/g\' '.$fo);
}

my $db = "mt_35";
my $dir = "$DIR_Genome/$db";
my $d00 = "$dir/00_seq";
#seqPreProcess($d00, $db);
#seqPostProcess($d00, $db);
my $f_mapping = "$d00/29.tbl";
my $f_seq_genome = "$d00/19.fa";
my $f_gff_genome = "$d00/31.gff";
my $d10 = "$dir/10_model_Mt3.5v5";
geneProcess(-dir=>$d10, -fm=>$f_mapping, -f_seq=>$f_seq_genome);
my $f_gff_gene = "$d10/71.gff";

my $d80 = "$dir/80_cp_mt";
#processCpMt(-dir=>$d80);
my $f_seq_cpmt = "$d80/41.fa";
my $f_gff_cpmt = "$d80/51.gff";

#my $ld = Localdb->new(-db=>$db);
#$ld->loadGff(-files=>[$f_gff_gene, $f_seq_cpmt, $f_gff_cpmt, $f_seq_genome, $f_gff_genome], -empty=>1);

my $f41 = "$dir/41_genome.fa";
#cat $f_seq_genome $f_seq_cpmt > $f41

my $f43 = "$dir/43_genome_snpEff.gff";
#merge_gene_genome($f41, $f_gff_gene, $f43);


