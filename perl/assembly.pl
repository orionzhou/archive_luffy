#!/usr/bin/perl -w
use strict;
use FindBin;
use lib "$FindBin::Bin";

use Common;
use Seq;
use Align;
use Blast;
use Gal;
use Mapping;
use Data::Dumper;
use File::Path qw/make_path remove_tree/;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;

my $data = "/home/youngn/zhoup/Data";

my $acc = "hm056";
#my $acc = "hm340";
my $dir = "$data/misc3/$acc";
my $f_seq = "$dir/01_assembly.fa";
my $f_len = "$dir/11_seqlen.tbl";
my $f_gap = "$dir/12_gaploc.tbl";
#print "seqlen.pl -out $f_len $f_seq\n";
#print "seqgap.pl -out $f_gap -min 100 $f_seq\n";

my $d21 = "$dir/21_blastn";
make_path($d21) unless -d $d21;
my $f21_02 = "$d21/02_raw.tbl";
my $f21_05 = "$d21/05_tiled.tbl";
#runCmd("blastn -db \$data/db/blast/mt4.0 -outfmt 6 -query $f_seq -out $f21_02\n", 1);
#runCmd("blastTiling -i $f21_02 -o $f21_05");

my $d23 = "$dir/23_blat";
-d $d23 || make_path($d23);
my $d23_01 = "$d23/01_seq";
-d $d23_01 || make_path($d23_01);
print "ln -sf ../../01_assembly.fa part.fa\n";
print "pyfasta split -n 160 part.fa\n";
my $d23_03 = "$d23/03_raw";
-d $d23_03 || make_path($d23_03);
my $f_ref = "$data/genome/Mtruncatula_4.0/11_genome.fa";
my $f_ref_2bit = "$data/db/blat/Mtruncatula_4.0.2bit";
#run comp.hm056

print "pslSort dirs 04.psl tmp 03_raw\n";
print "axtChain -linearGap=loose -psl 04.psl $data/db/blat/Mtruncatula_4.0.2bit $data/db/blat/HM056.2bit 53.chain\n";
print "chainMergeSort 53.chain > 55_mergesort.chain\n";
print "chainPreNet 55_mergesort.chain ../hm101.sizes ../$acc.sizes 57_pre.chain\n";
chain2Gal("$d23/57_pre.chain", "$d23/61.gal");
gal_filter("$d23/61.gal", "$d23/62_chain.gal", "$d23/63_block.gal", 300, 0.05, 0.05);


#print "pslReps 04.psl 05.psl 05.psr\n";
my $f23_05 = "$d23/05.psl";
my $f23_11 = "$d23/11.gal";
#psl2Gal($f23_05, $f23_11);
my $f23_12 = "$d23/12_tiled.gal";
#gal_tiling($f23_11, $f23_12);
my $f23_13 = "$d23/13_checked.gal";
#gal_rm_inv($f23_12, $f23_13);
my $f23_15 = "$d23/15.gal";
#gal_breakdown($f23_13, $f23_15, 100000);
my $f23_17 = "$d23/17_chain.gal";
my $f23_18 = "$d23/18_block.gal";
#gal_filter($f23_15, $f23_17, $f23_18, 300, 0.05, 0.05);
my $f23_19 = "$d23/19_block.gal";
#gal_validate($f23_18, $f23_19, $f_seq, $f_ref);
my $f23_21 = "$d23/21_indel.gal";
#gal_indel($f23_18, $f23_21);
my $f23_23a = "$d23/23_qnov.tbl";
my $f23_23b = "$d23/23_hnov.tbl";
#get_new_loc($f23_17, $f_len, $f23_23a, $f23_23b);
sub get_new_loc {
    my ($fi, $fl, $fo1, $fo2) = @_;
    my $t = readTable(-in=>$fi, -header=>1);
    my $tl = readTable(-in=>$fl, -header=>1);
    my $hl = { map {$tl->elm($_, "id") => $tl->elm($_, "length")} (0..$tl->nofRow-1) };
   
    open(FHO, ">$fo1") or die "cannot write $fo1\n";
    print FHO join("\t", qw/id beg end length/)."\n";
    my $hq;
    my $hh;
    for my $i (0..$t->nofRow-1) {
        my ($idc, $qId, $qBeg, $qEnd, $qSrd, $qLen,
            $hId, $hBeg, $hEnd, $hSrd, $hLen) = $t->row($i);
        $hq->{$qId} ||= [];
        push @{$hq->{$qId}}, [$qBeg, $qEnd];
        $hh->{$hId} ||= [];
        push @{$hh->{$hId}}, [$hBeg, $hEnd];
    }
    for my $id (keys(%$hl)) {
        my $locA = [[1, $hl->{$id}]];
        my $loc = $locA;
        if(exists $hq->{$id}) {
            ($loc) = posDiff($locA, posMerge($hq->{$id}));
        }
        for (@$loc) {
            my ($beg, $end) = @$_;
            print FHO join("\t", $id, $beg, $end, $end-$beg+1)."\n";
        }
    }
    close FHO; 
}

# run R script assembly.R
my $f21_21 = "$d21/21_scaffold_status.tbl";
my $f21_22 = "$d21/22_unmapped.tbl";

my $f21_31 = "$d21/31_unmapped.fa";
#print "seqextract.pl -out $f21_31 -id $f21_22 $f_seq\n";
my $f21_32 = "$d21/32_raw.tbl";
my $f21_35 = "$d21/35_tiled.tbl";
#print "blastn -db /project/db/blast/current/nt -outfmt 6 -query $f21_31 -out $f21_32\n";
#print "blastTiling -i $f21_32 -o $f21_35\n";
my $f21_36 = "$d21/36_annotated.tbl";
#annotate_nr_hits($f21_35, $f21_36);

# run R script assembly.R
my $f21_41 = "$d21/41_scaffold_status.tbl";
my $f21_42 = "$d21/42_unknown.tbl";

my $f21_51 = "$d21/51_unknown.fa";
#print "seqextract.pl -out $f21_51 -id $f21_42 $f_seq\n";

sub annotate_nr_hits {
    my ($fi, $fo) = @_;
    my $t = readTable(-in=>$fi, -header=>1);
    my @gis;
    for my $i (0..$t->nofRow-1) {
        my $id = $t->elm($i, "hId");
        if($id =~ /gi\|(\d+)\|/) {
            push @gis, $1;
        } else {
            die "unknown hId: $id\n";
        }
    }
    my ($cats, $species) = get_gi_taxonomy(@gis);
    $t->addCol($cats, "category");
    $t->addCol($species, "species");
    open(FH, ">$fo") or die "cannot write to $fo\n";
    print FH $t->tsv(1);
    close FH;
}


__END__

