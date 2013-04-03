#!/usr/bin/perl -w
use strict;

my ($fi, $fo) = @ARGV;
open(FHI, "<$fi") or die "cannot open $fi for reading\n";
open(FHO, ">$fo") or die "cannot open $fo for writing\n";
my $line = <FHI>;
print FHO join("\t", qw/probe chrom loc/)."\n";
while(<FHI>) {
    chomp;
    my ($xb, $yb, $x, $y, $id2, $id, $type, $org, $sym, $loc) = split("\t", $_);
    $loc =~ /[\d\.]+\:([\w\d]+)\:(\d+)\-(\d+)\:/;
    my ($chr, $beg, $end) = ($1, $2, $3);
    my $pos = int( ($beg+$end)/2 );
    print FHO join("\t", $id, $chr, $pos)."\n";
}
close FHI;
close FHO;

