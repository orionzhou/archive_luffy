#!/usr/bin/perl -w
use strict;
use lib ($ENV{"SCRIPT_HOME_PERL"});
use Common;
use Seq;
use Align;
use Blast;
use GenomeCompare;
use Eutils;
use Data::Dumper;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;

=begin
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
=end
=cut

my $dir = "/home/youngn/zhoup/Data/misc3/hm056/blastn";
my $f_seq = "$dir/../hm056.fa";
my $fl = "$dir/../hm056_seqlen.tbl";
my $f01 = "$dir/01_blast.tbl";
#blastn -db $data/db/blast/mt4.0 -outfmt 6 -query ../hm056.fa -out 01_blast.tbl
my $d02 = "$dir/02_raw";
my $f05 = "$dir/05_tiled.tbl";
#split_blast_output($f01, $d02);
#blast_tiling($fl, $d02, $f05);

my $f06 = "$dir/06.tbl";
#find_block($f05, $f06);
my $f08 = "$dir/08_long.tbl";
my $f09 = "$dir/09_wide.tbl";
#filter_block($f06, $fl, $f08, $f09);

my $f11 = "$dir/11_status.tbl";
my $f21 = "$dir/21_unmapped.fa";
#extract_unmapped($f11, $f_seq, $f21);
my $f23 = "$dir/23_nr.tbl";
my $d25 = "$dir/25_nr_raw";
my $f26 = "$dir/26_nr_tiled.tbl";
#split_blast_output($f23, $d25);
#blast_tiling($fl, $d25, $f26);
my $f31 = "$dir/31_block.tbl";
#find_block($f26, $f31);
my $f33 = "$dir/33_long.tbl";
my $f34 = "$dir/34_wide.tbl";
#filter_block($f31, $fl, $f33, $f34);
my $f36 = "$dir/36_desc.tbl";
#annotate_nr_hits($f34, $f36);

sub annotate_nr_hits {
    my ($fi, $fo) = @_;
    my $t = readTable(-in=>$fi, -header=>1);
    my @ids;
    for my $i (0..$t->nofRow-1) {
        my $id = $t->elm($i, "hId");
        if($id =~ /gi\|(\d+)\|/) {
            push @ids, $1;
        } else {
            die "unknown hId: $id\n";
        }
    }
    my @orgs = get_gi_note(@ids);
    $t->addCol(\@orgs, "desc");
    my @flags = map {$_ =~ /(medicago)|(truncatula)|(soybean)|(glycine)|(lotus)/ ? 1 : 0} @orgs;
    $t->addCol(\@flags, "org");
    open(FH, ">$fo") or die "cannot write to $fo\n";
    print FH $t->tsv(1);
    close FH;
}

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

__END__

