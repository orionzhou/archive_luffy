#!/usr/bin/perl
use strict;
use lib ($ENV{"SCRIPT_HOME_PERL"});
use Data::Dumper; 
use InitPath;
use Common;
use Seq;
use Gff;
use Gtb;
use Bed;

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

my $org = "Mtruncatula_4.0";
my $dir = "$DIR_genome/$org";
my $f01 = "$dir/01_refseq.fa";
print "seqconcat.pl\n";
print "seqlen.pl\n";
my $f61 = "$dir/61_global_loc.gff";
print "gff_convert_loc.pl\n";
my $f62 = "$dir/62_global_loc.gtb";
#gff2Gtb($f61, $f62);
my $f65 = "$dir/65_phase_fixed.gtb";
#gtb_fix_phase($f62, $f65, $f01);
my $f69 = $f65;

my $f71 = "$dir/71.gff";
#gtb2Gff($f69, $f71);
my $f81 = "$dir/81.bed";
#gtb2Bed($f69, $f81);
my $f82 = "$dir/82_refined.bed";
#bed_refine($f81, $f82);
my $f86 = "$dir/86.tbl";
#gtb2Tbl($f69, $f86);
my $f91 = "$dir/91_proteome.fa";
#gtb2Seq(-in=>$f69, -out=>$f91, -seq=>$f01, -opt=>2);


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


