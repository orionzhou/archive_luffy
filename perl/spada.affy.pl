#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use FindBin;
use lib $FindBin::Bin;
use File::Path qw/make_path remove_tree/;
use Data::Dumper;
use List::Util qw/min max sum/;

my $org = "Athaliana";
#$org = "Mtruncatula_3.5";
my $dir = "/home/youngn/zhoup/Data/misc2/spada.affy/$org";

my $f01 = "$dir/01_probe.tbl";
my $f02 = "$dir/02_probe.fa";
#get_probe_seq("$dir/00.tbl", $f01, $f02);

sub get_probe_seq {
    my ($fi, $fo1, $fo2) = @_;
    my ($fhi, $fho);
    open($fhi, "<$fi") or die "cannot read $fi\n";
    open($fho, ">$fo1") or die "cannot write $fo1\n";
    my $seqHO = Bio::SeqIO->new(-file=>">$fo2", -format=>'fasta');
    
    print $fho join("\t", qw/id set seq/)."\n";
    my $h;
    while(<$fhi>) {
        chomp;
        my ($set, $num1, $num2, $num3, $seq) = split "\t";
        $h->{$set} ||= 0;
        $h->{$set} ++;
        my $id = sprintf "%s.%02d", $set, $h->{$set};
        $seqHO->write_seq( Bio::Seq->new(-id=>$id, -seq=>$seq) );
        print $fho join("\t", $id, $set, $seq)."\n";
    }
    $seqHO->close();
    close $fho;
}


