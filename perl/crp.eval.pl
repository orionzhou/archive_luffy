#!/usr/bin/perl
use strict;
use lib $ENV{"SCRIPT_HOME_PERL"};
use Data::Dumper;
use List::Util qw/min max sum/;
use POSIX qw/ceil floor/;
use Time::HiRes qw/gettimeofday tv_interval/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;

use InitPath;
use Common;
use CompareModel;

my $org = "Athaliana";
$org = "Osativa";
my $dir = "$DIR_Misc4/spada.crp.$org";
my $f_gtb_gs = "$dir/41_perf_eval/01_model.gtb";

my $f_gtb_cur = "$dir/01_preprocessing/61_gene.gtb";
my $f_gtb_una = "$DIR_Misc2/crp.ssp/$org/08.gtb";
#eval_gtb($f_gtb_cur, $f_gtb_gs);
eval_gtb($f_gtb_una, $f_gtb_gs);

sub model_eval {
    my ($f_eval, $f_ref_gtb, $fo) = @_;
    my $tv = readTable(-in=>$f_eval, -header=>1);
    my $tr = readTable(-in=>$f_ref_gtb, -header=>1);

    my $hg;
    for my $i (0..$tr->nofRow-1) {
        my ($id, $locS) = map {$tr->elm($i, $_)} qw/id locC/;
        my $loc = locStr2Ary($locS);
        my $len = locAryLen($loc);
        my $n_cds = @$loc;
        $hg->{$id} = [0, $len, $n_cds];
    }

    open(FH, ">$fo") || die "cannot open $fo for writing\n";
    print FH join("\t", qw/id tag gene lenTP lenFP lenFN exonTP exonFP exonFN/)."\n";
    for my $i (0..$tv->nofRow-1) {
        my ($id, $gene, $tag, $lenTP, $lenFP, $lenFN, $exonTP, $exonFP, $exonFN) = $tv->row($i);
        $tag = 5 if $tag == 7 || $tag == 8 || ($tag==2 && $lenFP+$lenFN>30);

        print FH join("\t", $id, $tag, $gene, $lenTP, $lenFP, $lenFN, $exonTP, $exonFP, $exonFN)."\n";
        $hg->{$gene}->[0] ++ if $gene ne "";
    }
    for my $gene (keys(%$hg)) {
        my ($cnt, $len, $n_cds) = @{$hg->{$gene}};
        if($cnt == 0) {
            print FH join("\t", '', 10, $gene, 0, 0, $len, 0, 0, $n_cds)."\n";
        } elsif($cnt > 1) {
            print "  $gene hit $cnt times\n";
        }
    }
    close FH;
}
sub get_sn_sp {
    my ($fi) = @_;
    my $t = readTable(-in=>$fi, -header=>1);
    my $sn_nt = sum($t->col("lenTP")) / ( sum($t->col("lenTP")) + sum($t->col("lenFN")) );
    my $sp_nt = sum($t->col("lenTP")) / ( sum($t->col("lenTP")) + sum($t->col("lenFP")) );
    my $sn_ex = sum($t->col("exonTP")) / ( sum($t->col("exonTP")) + sum($t->col("exonFN")) );
    my $sp_ex = sum($t->col("exonTP")) / ( sum($t->col("exonTP")) + sum($t->col("exonFP")) );
    return ($sn_nt, $sp_nt, $sn_ex, $sp_ex);
}
sub eval_gtb {
    my ($f_gtb_qry, $f_gtb_ref) = @_;
    my $f_eval = "eval.tbl";
    compare_models($f_gtb_qry, $f_gtb_gs, $f_eval);
    my $f_stat = "stat.tbl";
    model_eval($f_eval, $f_gtb_gs, $f_stat);
    my ($sn_nt, $sp_nt, $sn_ex, $sp_ex) = get_sn_sp($f_stat);
    print join(" ", $sn_nt, $sp_nt, $sn_ex, $sp_ex)."\n";
    runCmd("rm -f $f_eval $f_stat");
}



