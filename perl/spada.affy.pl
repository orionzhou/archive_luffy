#!/usr/bin/perl -w
use strict;
use FindBin;
use lib $FindBin::Bin;
use File::Path qw/make_path remove_tree/;
use InitPath;
use Mapping;
use Common;
use Data::Dumper;

my $pre = "Mtruncatula_3.5";
my $dir = "$DIR_misc2/spada.affy/$pre";
my $f01 = "$dir/01_probeset.tbl";
my $f05 = "$dir/05_seq.fa";
#get_probeseq($f01, $f05);

sub get_probeseq {
    my ($fi, $fo) = @_;
    my $t = readTable(-in=>$fi, -header=>1);
    my $h;
    open(FHO, ">$fo") or die "cannot open $fo for writing\n";
    for my $i (0..$t->nofRow-1) {
        my ($ids, $seq) = $t->row($i);
        $h->{$ids} ||= 0;
        $h->{$ids} ++;
        my $id = sprintf "%s.%02d", $ids, $h->{$ids};
        print FHO join("\n", ">$id", $seq)."\n";
    }
    close FHO;
}

