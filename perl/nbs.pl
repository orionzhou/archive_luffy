#!/usr/bin/perl
use strict;
use Init;
use Common;
use Gtb;
use Align;

my $f_ref = "$DIR_Genome/mt_35/41_genome.fa";
my $f_gtb = "$DIR_Genome/mt_35/10_model_Mt3.5v5/62_phase_fixed.gtb";

my $dir = "$DIR_Misc2/nbs/mt_35";
my $f00 = "$dir/00_id_gene.tbl";
my $f01 = "$dir/01.gtb";
#extract_gtb_by_ids($f00, $f_gtb, $f01);
my $f06 = "$dir/06.gff";
#gtb2Gff($f01, $f06);
my $f07 = "$dir/07.bed";
#gtb2Bed($f01, $f07);

my $f11 = "$dir/11_seq.fa";
#gtb2Seq(-in=>$f01, -out=>$f11, -seq=>$f_ref, -opt=>2);

my $f16 = "$dir/16.aln";
#run_clustalo(-in=>$f11, -out=>$f16);

my $f17 = "$dir/17.phb";
#convertTree(-in=>$f021, -out=>$f022, -informat=>'newick', -outformat=>'newick');

=cut
#run_hmmscan(-in=>$f03, -out=>$f04);
my $f05 = file($dir, "05_domain.txt");
#hmmParse(-in=>$f04, -out=>$f05);
#hmmscanTiling(-in=>$f04, -e=>1e-5, -inc=>"LRR", -out=>$f05);

my $f06 = dir($dir, "06_meme");
my $param_meme = {flag=>["protein"], mod=>"anr", nmotifs=>3, maxsize=>1000000};
#run_meme(-in=>$f03, -out=>$f06, -param=>$param_meme);

my $f11 = file($dir, "../../11_mt1_anno.txt");
my $f12  = file($dir, '12_anno.tlf');
#annForNBS(-in=>$f11, -out=>$f12, -db=>'mt_nbs', -ids=>\@ids, -opt=>2);

my $f081  = file($dir, '08_anno.txt');
my $f082  = file($dir, '08_anno.tlf');
#annProbe(-out=>$f081, -ids=>$param->{ids}, -opt=>1);
#annProbe(-out=>$f082, -ids=>$param->{ids}, -opt=>2);
my $f091  = file($dir, '09_anno.txt');
my $f09   = file($dir, '09_anno.tlf');

#annHclust(-out=>$f09, -in1=>$f081, -in2=>$f091);
my $f101 = file($dir, '10_exp.nwk');
=cut

sub extract_gtb_by_ids {
    my ($f00, $f_gtb, $f01) = @_;
    my $t0 = readTable(-in=>$f00, -header=>1);
    my $h = readGtb(-in=>$f_gtb, -opt=>2);

    open(FH, ">$f01");
    print FH join("\t", qw/id parent chr beg end strand locE locI locC loc5 loc3 phase source conf cat1 cat2 cat3 note/)."\n";
    for my $i (0..$t0->nofRow-1) {
        my ($id_g, $cat) = $t0->row($i);
        my $id = "$id_g.1";
        if(exists $h->{$id}) {
            my @ps = @{$h->{$id}};
            $ps[16] = $cat;
            print FH join("\t", @ps)."\n";
        } else {
            print "$id_g not found in [$f_gtb]\n";
        }
    }
    close FH;
}



