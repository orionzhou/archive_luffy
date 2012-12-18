#!/usr/bin/perl -w
use strict;
use Init;
use Common;
use Data::Dumper; 

my $dir = "$DIR_Misc2/crp";


sub cmp_hits {
    my ($f1, $f2, $fo) = @_;
    my $t1 = readTable(-in=>$f1, -header=>1);
    my $t2 = readTable(-in=>$f2, -header=>1);
    my $h2;
    for my $j (0..$t2->nofRow-1) {
        my ($chr, $beg, $end, $strand) = map {$t2->elm($j, $_)} qw/chr start end strand/;
        $h2->{$chr}->{$strand} ||= [];
        push @{$h2->{$chr}->{$strand}}, [$beg, $end];
    }

    open(FH, ">$fo");
    print FH join("\t", "id1", "tag", "id2")."\n";
    for my $i (0..$t1->nofRow-1) {
        my ($chr, $beg, $end, $strand) = map {$t1->elm($j, $_)} qw/chr start end strand/;
        my $locQ = [[$beg, $end]];
        my $locT = $h2->{$chr}->{$strand};
    }
    close FH;
}
sub crp_hit_pick_longest_cds {
#my $fi = "$dir/04_hits_picked/20_final.tbl";
#my $fo = "$dir/03_hits/05_picked_201112.tbl";
#tmp_fun($fi, $fo);
    my ($fi, $fo) = @_;
    my $ti = readTable(-in=>$fi, -header=>1);
    
    open(FH, ">$fo");
    print FH join("\t", qw/id1 id2 family chr beg end strand location e source/)."\n";
    for my $i (0..$ti->nofRow-1) {
        my ($id1, $id2, $fam, $chr, $locStr, $source, $e) = $ti->row($i);
        my $locObj = locStr2Obj($locStr);
        my $strand = $locObj->strand;
        my $locAry_o = locObj2Ary($locObj);
        my ($beg, $end) = get_longest_cds($locAry_o, $strand);
        my $locStr_n = "$beg..$end";
        $locStr_n = "complement($locStr_n)" if $strand == -1;
        print FH join("\t", $id1, $id2, $fam, $chr, $beg, $end, $strand, $locStr_n, $e, $source)."\n";
    }
    close FH;
}




