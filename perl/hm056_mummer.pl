#!/usr/bin/perl -w
use strict;
use lib ($ENV{"SCRIPT_HOME_PERL"});
use Common;
use Seq;
use Align;
use Blast;
use Data::Dumper;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;

my $dir = "/home/youngn/zhoup/Data/misc3/hm056/mummer";
my $f_seq = "$dir/../hm056.fa";
my $fl = "$dir/../hm056_seqlen.tbl";
#$src/MUMmer3.23/nucmer -p chr1-4 ../mt4_chr1-4.fa ../HM056.ALLPATHSLG.Jan2013.fasta
#$src/MUMmer3.23/nucmer -p chr5-8 ../mt4_chr5-8.fa ../HM056.ALLPATHSLG.Jan2013.fasta
#$src/MUMmer3.23/show-coords chr1-4.delta > chr1-4.coords
#$src/MUMmer3.23/show-coords chr5-8.delta > chr5-8.coords
#cat chr1-4.coords chr5-8.coords > 01_mummer.coords
my $f01 = "$dir/01_mummer.coords";
my $f02 = "$dir/02_coords.tbl";
#mummer_coords2tbl($f01, $f02);
my $f05 = "$dir/05_tiled.tbl";
#mummer_tiling($f02, $f05);
my $f06 = "$dir/06.tbl";
#find_block($f05, $f06);
my $f08 = "$dir/08_long.tbl";
my $f09 = "$dir/09_wide.tbl";
#filter_block($f06, $fl, $f08, $f09);

my $f11 = "$dir/11_status.tbl";
my $f21 = "$dir/21_unmapped.fa";
#extract_unmapped($f11, $f_seq, $f21);
my $d23 = "$dir/23_nr";
blast_nr_batch($f21, $d23);

sub mummer_coords2tbl {
    my ($fi, $fo) = @_;
    open(FHI, "<$fi") or die "cannot open $fi for reading\n";
    open(FHO, ">$fo") or die "cannot open $fo for writing\n";
    print FHO join("\t", qw/qId qBeg qEnd strand hId hBeg hEnd qLen hLen pct_idty/)."\n";
    my $i = 1;
    while(<FHI>) {
        chomp;
        my $line = $_;
        my @ps = split /\|/, $line;
        next if @ps != 5 || $ps[0] =~ /\[S1\]/;
        my ($hBeg, $hEnd) = split(" ", $ps[0]);
        my ($qBeg, $qEnd) = split(" ", $ps[1]);
        my $strd = $qBeg > $qEnd ? "-" : "+";
        ($qBeg, $qEnd) = ($qEnd, $qBeg) if $strd eq "-";
        my ($hLen, $qLen) = split(" ", $ps[2]);
        my ($hId, $qId) = split(" ", $ps[4]);
        my $pct = $ps[3];
        $pct =~ s/(^ +)|( +$)//g;
        die "$line\n" if $hEnd-$hBeg+1 != $hLen || $qEnd-$qBeg+1 != $qLen;
        print FHO join("\t", $qId, $qBeg, $qEnd, $strd, $hId, $hBeg, $hEnd, $qLen, $hLen, $pct)."\n";
    }
    close FHI;
    close FHO;
}
sub mummer_tiling {
    my ($fi, $fo) = @_;
    my $t = readTable(-in=>$fi, -header=>1);
    $t->sort("qId", 1, 0, "qBeg", 0, 0);
    my $ref = group($t->colRef("qId"));

    open(FHO, ">$fo") or die "cannot write to $fo\n";
    print FHO join("\t", $t->header)."\n";
    for my $qId (sort keys(%$ref)) {
        my ($idx, $cnt) = @{$ref->{$qId}};
        my @locs = map {[$t->elm($_, "qBeg"), $t->elm($_, "qEnd")]} ($idx..$idx+$cnt-1);
        my @stats = map {$t->elm($_, "qLen")} ($idx..$idx+$cnt-1);
        my $refs = tiling(\@locs, \@stats, 2);
        for (@$refs) {
            my ($beg, $end, $idxR) = @$_;
            next if ($end - $beg + 1) < 1000;
            my $idxA = $idxR + $idx;
            my ($qId2, $qBeg, $qEnd, $strand, $hId, $hBeg, $hEnd, $qLen, $hLen, $pct) 
                = $t->row($idxA);
            my $begH = sprintf "%d", $hBeg + ($beg-$qBeg) * ($hEnd-$hBeg)/($qEnd-$qBeg);
            my $endH = sprintf "%d", $hEnd - ($qEnd-$end) * ($hEnd-$hBeg)/($qEnd-$qBeg);
            $qLen = $end - $beg + 1;
            $hLen = $endH - $begH + 1;
            print FHO join("\t", $qId, $beg, $end, $strand, $hId, $begH, $endH, $qLen, $hLen, $pct)."\n";
        }
    }
    close FHO;
}
sub assign_block {
    my ($t) = @_;

    my $clu = 0;
    my ($qBegP, $qEndP, $strdP, $hIdP, $hBegP, $hEndP);
    for my $i (0..$t->nofRow-1) {
        next if $t->elm($i, "clu") > 0;
        $t->setElm($i, "clu", ++$clu);
        ($qBegP, $qEndP, $strdP, $hIdP, $hBegP, $hEndP) = 
            map {$t->elm($i, $_)} qw/qBeg qEnd strand hId hBeg hEnd/;
        for my $j ($i+1..$t->nofRow-1) {
            my ($qBeg, $qEnd, $strd, $hId, $hBeg, $hEnd) = 
                map {$t->elm($j, $_)} qw/qBeg qEnd strand hId hBeg hEnd/;
            if( abs($qBeg-$qEndP)<10000 && $hId eq $hIdP && $strd eq $strdP &&
                ( ($strd eq "+" && abs($hBeg-$hEndP)<1000000) ||
                  ($strd eq "-" && abs($hEnd-$hBegP)<1000000) ) ) {
                $t->setElm($j, "clu", $clu);
                ($qBegP, $qEndP, $hIdP, $hBegP, $hEndP) = ($qBeg, $qEnd, $hId, $hBeg, $hEnd);
            }
        }
    }
    return $t;
}
sub find_block {
    my ($fi, $fo) = @_;
    my $ti = readTable(-in=>$fi, -header=>1);
    $ti->sort("qId", 1, 0, "qBeg", 0, 0);
    $ti->addCol([(0)x$ti->nofRow], "clu");

    open(FHO, ">$fo") or die "cannot open $fo for writing\n";
    print FHO join("\t", $ti->header)."\n";
    my $ref = group($ti->colRef("qId"));
    for my $qId (sort(keys(%$ref))) {
        my ($idx, $cnt) = @{$ref->{$qId}};
        my $t = $ti->subTable([$idx..$idx+$cnt-1], [$ti->header]);
        $t = assign_block($t);
        print FHO $t->tsv(0);
    }
    close FHO;
}
sub filter_block {
    my ($fi, $fl, $fo1, $fo2) = @_;
    
    my $tl = readTable(-in=>$fl, -header=>1);
    my $hl = { map {$tl->elm($_, "id") => $tl->elm($_, "length")} (0..$tl->nofRow-1) };
    
    my $t = readTable(-in=>$fi, -header=>1);
    $t->sort("qId", 1, 0, "clu", 0, 0, "qBeg", 0, 0);
    
    open(FHO1, ">$fo1") or die "cannot write to $fo1\n";
    print FHO1 join("\t", $t->header)."\n";
    open(FHO2, ">$fo2") or die "cannot write to $fo2\n";
    print FHO2 join("\t", qw/qId qBeg qEnd strand hId hBeg hEnd qLen hLen clu lenS cov/)."\n";
    my $ref1 = group($t->colRef("qId"));
    for my $qId (sort keys %$ref1) {
        my ($idx1, $cnt1) = @{$ref1->{$qId}};
        my $t1 = $t->subTable([$idx1..$idx1+$cnt1-1], [$t->header]);
        my $lenS = $hl->{$qId};

        my $ref2 = group($t1->colRef("clu"));
        for my $clu (sort {$a<=>$b} keys %$ref2) {
            my ($idx2, $cnt2) = @{$ref2->{$clu}};
            my $t2 = $t1->subTable([$idx2..$idx2+$cnt2-1], [$t1->header]);
            my $qBeg = min($t2->col("qBeg"));
            my $qEnd = max($t2->col("qEnd"));
            my $hBeg = min($t2->col("hBeg"));
            my $hEnd = max($t2->col("hEnd"));
            my $qLoc = [map {[$t2->elm($_, "qBeg"), $t2->elm($_, "qEnd")]} (0..$t2->nofRow-1)];
            my $hLoc = [map {[$t2->elm($_, "hBeg"), $t2->elm($_, "hEnd")]} (0..$t2->nofRow-1)];
            my $qLen = locAryLen(posMerge($qLoc));
            my $hLen = locAryLen(posMerge($hLoc));
            my ($strand, $hId) = map {$t2->elm(0, $_)} qw/strand hId/;
            my $cov = sprintf "%.03f", $qLen / $lenS;
            next if $cov < 0.2;
            print FHO1 $t2->tsv(0);
            print FHO2 join("\t", $qId, $qBeg, $qEnd, $strand, $hId, $hBeg, $hEnd,
                $qLen, $hLen, $clu, $lenS, $cov)."\n";
        }
    }
}
sub extract_unmapped {
    my ($fi, $f_seq, $fo) = @_;
    my $t = readTable(-in=>$fi, -header=>1);
    $t = $t->match_pattern("\$_->[2] < 2");
    my %h = map {$_=>1} $t->col("id");
   
    my $cnt = 0;
    my $seqHI = Bio::SeqIO->new(-file=>"<$f_seq", -format=>'fasta');
    my $seqHO = Bio::SeqIO->new(-file=>">$fo", -format=>'fasta');
    while(my $seqO = $seqHI->next_seq()) {
        my ($id) = ($seqO->id);
        if(exists $h{$id}) {
            $seqHO->write_seq($seqO);
            $cnt ++;
        }
    }
    $seqHI->close();
    $seqHO->close();
    printf "  %4d sequences extracted\n", $cnt;
}

__END__

