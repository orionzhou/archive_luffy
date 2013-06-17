#!/usr/bin/perl
use strict;
use FindBin;
use lib "$FindBin::Bin";
use InitPath;
use Common;
use Data::Dumper;
use Hmm;

my $ps = {
    "01_crp_hmmsearchP" => {opt=>'p', min_e=>1,  min_len=>10},
    "02_crp_hmmsearchX" => {opt=>'x', min_e=>10, min_len=>10},
};
my $tag = "02_crp_hmmsearchX";
my $p = $ps->{$tag};
my $dir = "$DIR_Misc2/hmmsearch/$tag";
system("mkdir -p $dir") unless -d $dir;

my $f_qry = "$dir/00_seq.fa"; 
my $f_hmm = "$dir/00.hmm";
my $f_ref = "$DIR_Genome/mt_35/41_genome.fa";
my $f_gtb = "$DIR_Genome/mt_35/10_model_Mt3.5v5/62_frame_fixed.gtb";

my $f01 = "$dir/01_raw.txt";
#hmmRun($f_hmm, $f_qry, $f01);
my $f07 = "$dir/07.bes";
hmmPostProcess($f01, $f_qry, $f_ref, $f_gtb, $p);


