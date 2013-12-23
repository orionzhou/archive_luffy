#!/usr/bin/perl -w
use strict;
use FindBin;
use lib "$FindBin::Bin";
use Bio::DB::Fasta;
use Time::HiRes qw/gettimeofday tv_interval/;

my $fi = "/home/youngn/zhoup/Data/misc3/pan3probe/02_wins.tbl";
my $fs = "/home/youngn/zhoup/Data/genome/pan3/11_genome.fa";

my $t0 = [gettimeofday];

open (my $fhi, $fi) || die "Can't open file $fi: $!\n";
my $db = Bio::DB::Fasta->new($fs);
my $cnt = 0;
while(<$fhi>) {
    chomp;
    my @ps = split "\t";
    my ($seqid, $beg, $end) = @ps;
    my $seq = $db->seq($seqid, $beg, $end);
    $cnt ++;
    printf "%d: %.01f min\n", $cnt, tv_interval($t0, [gettimeofday]) / 60 if $cnt % 100000 == 0;
}
close $fhi;

