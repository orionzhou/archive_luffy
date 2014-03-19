#!/usr/bin/perl -w
use strict;
use FindBin;
use lib $FindBin::Bin;
use File::Path qw/make_path remove_tree/;
use Common;
use Data::Dumper;
use List::Util qw/min max sum/;

my $dir = "/home/youngn/zhoup/Data/misc2/rnaseq/mt";
my $f01 = "$dir/01.tbl";
my $d11 = "$dir/11_fastq";
my $d21 = "$dir/21_tophat";


run_tophat($f01, $d11, $d21);

sub run_tophat {
  my ($f01, $d11, $d21) = @_;
  my $t = readTable(-in=>$f01, -header=>1);
  -d $d21 || make_path($d21);
  for my $i (0..$t->lastRow) {
#    next if $i >= 27;
    my ($sam, $id, $label, $note1, $note2) = $t->row($i);
    my $fa = "$d11/$id\_$label\_R1_001.fastq";
    my $fb = "$d11/$id\_$label\_R2_001.fastq";
    -s $fa || die "$fa is not there\n";
    -s $fb || die "$fb is not there\n";
    print "mapping $id\_$label to $sam\n";
    runCmd("tophat2 --num-threads 16 --mate-inner-dist 0 \\
      --min-intron-length 45 --max-intron-length 5000 \\
      --min-segment-intron 45 --max-segment-intron 5000 \\
      -o $d21/$id \$data/db/bowtie2/$sam $fa $fb", 1);
  }
}


