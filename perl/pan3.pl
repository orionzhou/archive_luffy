#!/usr/bin/perl -w
use strict;
use FindBin;
use lib "$FindBin::Bin";

use Common;
use Seq;
use Eutils;
use Data::Dumper;
use File::Path qw/make_path remove_tree/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;

my $data = "/home/youngn/zhoup/Data";
my $dir = "$data/misc3/pan3";
my $org = "HM340";
#make_chr("$dir/$org/01.tbl", "$data/genome/$org/11_genome.fa", $org, "$dir/$org");

# blastn -db /project/db/blast/current/nt -outfmt '6 qseqid qstart qend sseqid sstart send length nident mismatch gaps evalue bitscore qseq sseq' -num_threads 4 -query 01.fa -out 11_blast.tbl

annotate_blast_nr("$dir/$org/14.gal", "$dir/$org/15.gal", "$dir/$org/16_cat.tbl");

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
sub annotate_blast_nr {
    my ($fi, $fo1, $fo2) = @_;

    my $t = readTable(-in=>$fi, -header=>1);
    my @gis;
    for my $i (0..$t->nofRow-1) {
        my $id = $t->elm($i, "tId");
        if($id =~ /gi\|(\d+)\|/) {
            push @gis, $1;
            $t->setElm($i, "tId", $1);
        } else {
            die "unknown tId: $id\n";
        }
    }
    my $h1 = gi2Taxid(@gis);
    my $h2 = annotate_taxid(values(%$h1));
    
    open(FH1, ">$fo1") or die "cannot write to $fo1\n";
    print FH1 $t->tsv(1);
    close FH1;
    
    open(FH2, ">$fo2") or die "cannot write to $fo2\n";
    print FH2 join("\t", qw/id superkingdom kingdom family species/)."\n";
    for my $id (sort(keys(%$h1))) { 
        my $taxid = $h1->{$id};
        my @cats = @{$h2->{$taxid}};
        print FH2 join("\t", $id, @cats)."\n";
    }
    close FH2;
}


