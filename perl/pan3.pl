#!/usr/bin/perl -w
use strict;
use FindBin;
use lib "$FindBin::Bin";
use Data::Dumper;
use File::Path qw/make_path remove_tree/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
use Common;
use Seq;

my $data = "/home/youngn/zhoup/Data";
my $dir = "$data/misc3/pan3";
my $org = "HM340";

make_chr("$dir/$org.tbl", "$data/genome/$org/11_genome.fa", $org, "$dir/$org");
make_chr("$dir/$org\_unc.tbl", "$data/genome/$org/11_genome.fa", "$org\_unc", "$dir/$org\_unc");

sub make_chr {
    my ($fi, $fs, $tag, $fop) = @_;
    my $foa = "$fop.agp";
    
    open(FHO, ">$foa") or die "cannot write to $foa\n";
    print FHO join("\t", qw/id beg end srd sid sbeg send len/)."\n";
    my $t = readTable(-in=>$fi, -header=>1);
    my $seq = "";
    my $l = 0;
    for my $i (0..$t->nofRow-1) {
        my ($id, $beg, $end, $len) = $t->row($i);
        my $seq1 = seqRet([[$beg, $end]], $id, "+", $fs);
        if($i > 0) {
            $seq .= "N" x 50;
            $l += 50;
        }
        print FHO join("\t", $tag, $l+1, $l+$len, "+", $id, $beg, $end, $len)."\n";
        $seq .= $seq1;
        $l += $len;
    }
    close FHO;

    my $seqH = Bio::SeqIO->new(-file=>">$fop.fa", -format=>"fasta");
    $seqH->write_seq(Bio::Seq->new(-id=>$tag, -seq=>$seq));
    $seqH->close();
}



