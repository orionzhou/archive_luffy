#!/usr/bin/perl -w
use strict;
use lib ($ENV{"SCRIPT_HOME_PERL"});
use Bio::SeqIO;
use Common;
use Seq;
use Align;
use Data::Dumper;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;

my $dir = "/home/youngn/zhoup/Data/misc3/hm056/blastn";
my $f_seq = "$dir/../hm056.fa";
my $fl = "$dir/../hm056_seqlen.tbl";
my $f01 = "$dir/01_blast.tbl";
#blastn -db $data/db/blast/mt4.0 -outfmt 6 -query ../hm056.fa -out 01_blast.tbl
my $d02 = "$dir/02_raw";
#split_blast_output($f01, $d02);
my $f05 = "$dir/05_tiled.tbl";
#blast_tiling($fl, $d02, $f05);
my $f06 = "$dir/06.tbl";
#find_block($f05, $f06);
my $f08 = "$dir/08_long.tbl";
my $f09 = "$dir/09_wide.tbl";
#filter_block($f06, $fl, $f08, $f09);
sub write_blast {
    my ($id, $rows, $fo) = @_;
    open(my $fho, ">$fo") or die "cannot open $fo for writing\n";
    print $fho join("\t", qw/qId qBeg qEnd strand hId hBeg hEnd pct e score/)."\n";
    print $fho join("\n", map {join("\t", @$_)} @$rows)."\n";
    close $fho;
}
sub split_blast_output {
    my ($fi, $dirO) = @_;
    
    open(my $fhi, "<$fi") or die "cannot open $fi for reading\n";
    my ($idP, $lineP) = ("", "");
    my @rows;
    while(my $line = <$fhi>) {
        chop($line);
        my ($qId, $hId, $pct, $alnLen, $mm, $gap, $qBeg, $qEnd, $hBeg, $hEnd, $e, $score) = split("\t", $line);
        my $strand = $hBeg > $hEnd ? "-" : "+";
        ($hBeg, $hEnd) = ($hEnd, $hBeg) if $strand eq "-";
        my $row = [$qId, $qBeg, $qEnd, $strand, $hId, $hBeg, $hEnd, $pct, $e, $score];
        if($qId eq $idP) {
            push @rows, $row;
        } else {
            write_blast($idP, \@rows, "$dirO/$idP.tbl") if $idP;
            @rows = ($row);
            $idP = $qId;
        }
    }
    write_blast($idP, \@rows, "$dirO/$idP.tbl") if @rows > 0;
}
sub blast_tiling {
    my ($fl, $dirI, $fo) = @_;
    open(my $fho, ">$fo") or die "cannot open $fo for writing\n";
    print $fho join("\t", qw/qId qBeg qEnd strand hId hBeg hEnd qLen hLen pct e score/)."\n";
   
    my $tl = readTable(-in=>$fl, -header=>1);
    my @ids = $tl->col("id");
    my ($cntB, $cntG) = (0, 0);
    for my $id (sort @ids) {
        my $fi = "$dirI/$id.tbl";
        if( ! -s $fi ) {
            $cntB ++;
            next;
        }
       
        my $ti = readTable(-in=>$fi, -header=>1);
        my @locs = map { [$ti->elm($_, "qBeg"), $ti->elm($_, "qEnd")] } (0..$ti->nofRow-1);
        my @es = $ti->col("e");
        my $refs = tiling(\@locs, \@es, 1);
        for (@$refs) {
            my ($beg, $end, $idx) = @$_;
            next if ($end - $beg + 1) < 100;
            my ($qId, $qBeg, $qEnd, $strd, $hId, $hBeg, $hEnd, $pct, $e, $score) = $ti->row($idx);
            my $begH = sprintf "%d", $hBeg + ($beg-$qBeg) * ($hEnd-$hBeg)/($qEnd-$qBeg);
            my $endH = sprintf "%d", $hEnd - ($qEnd-$end) * ($hEnd-$hBeg)/($qEnd-$qBeg);
            my $qLen = $end - $beg + 1;
            my $hLen = $endH - $begH + 1;
            print $fho join("\t", $id, $beg, $end, $strd, $hId, $begH, $endH, $qLen, $hLen, $pct, $e, $score)."\n";
        }
        printf " %4d: %10d\n", (++$cntG)+$cntB, $ti->nofRow;
    }
    close $fho;
    printf "\n %4d / %4d tiled\n", $cntG, $cntG+$cntB;
}
