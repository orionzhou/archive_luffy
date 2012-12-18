#!/usr/bin/perl -w
use strict; use Path::Class; use Data::Dumper;
use Init; use Medicago; 

my $genome = "mt_35";
my $f_refseq = file($DIR_Genome, $genome, "41_genome.fa");
my $dir = dir($DIR_Repo, $genome, "30_vnt");
my $d01 = dir($dir, "01_bamlist");
my $d11 = dir($dir, "11_bcf");
my $d13 = dir($dir, "13_vcf");

my $opt = "acc26";
my $f01 = file($d01, "$opt.txt");
writeBamList($opt, $f01);

sub writeBamList {
    my ($opt, $fo) = @_;
    my $dir_bam = dir($DIR_Repo, $genome, "11_pipe_bwa/13_dup_removed");
    my $ids = get_acc_ids($opt);
    my $fh = new IO::File $fo, "w";
    for my $id (@$ids) {
        print $fh "$dir_bam/$id.bam\n";
    }
    close $fh;
}


