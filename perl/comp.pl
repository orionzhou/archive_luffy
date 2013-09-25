#!/usr/bin/perl -w
use strict;
use FindBin;
use lib "$FindBin::Bin";

use Common;
use Align;
use Blast;
use Gal;
use Mapping;
use Data::Dumper;
use File::Path qw/make_path remove_tree/;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;

my $data = "/home/youngn/zhoup/Data";

my $org_q = "HM056";
#my $org_q = "HM340";
my $org_t = "HM101";

my $f_seq_q = "$data/genome/$org_q/11_genome.fa";
my $f_2bit_q = "$data/db/blat/$org_q.2bit";
my $f_len_q = "$data/genome/$org_q/15_seqlen.tbl";
my $f_gap_q = "$data/genome/$org_q/16_gaploc.tbl";

my $f_seq_t = "$data/genome/$org_t/11_genome.fa";
my $f_2bit_t = "$data/db/blat/$org_t.2bit";
my $f_len_t = "$data/genome/$org_t/15_seqlen.tbl";
my $f_gap_t = "$data/genome/$org_t/16_gaploc.tbl";

my $dir = "$data/misc3/$org_q\_$org_t";

my $d21 = "$dir/21_blastn";
-d $d21 || make_path($d21);

my $d23 = "$dir/23_blat";
-d $d23 || make_path($d23);
#run comp.hm056

#chain2Gal("$d23/17.chain", "$d23/17.gal");
#chain2Gal("$d23/27.chain", "$d23/27.gal");
#gal_filter("$d23/61.gal", "$d23/62_chain.gal", "$d23/63_block.gal", 300, 0.05, 0.05);
#net2Gal("$d23/61.net", "$d23/71.gal");

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

__END__

