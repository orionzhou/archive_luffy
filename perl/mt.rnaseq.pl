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
my $d12 = "$dir/12_fastq_gz";
my $d21 = "$dir/21_tophat";

#gzip_fastq($f01, $d11, $d12);
run_tophat($f01, $d12, $d21);

sub gzip_fastq {
  my ($f01, $d11, $d12) = @_;
  -d $d12 || make_path($d12);
  my $t = readTable(-in=>$f01, -header=>1);
  for my $i (0..$t->lastRow) {
#    next if $i < 13;
    my ($sam, $id, $label, $note1, $note2) = $t->row($i);
    my $fa = "$d11/$id\_$label\_R1_001.fastq";
    my $fb = "$d11/$id\_$label\_R2_001.fastq";
    -s $fa || die "$fa is not there\n";
    -s $fb || die "$fb is not there\n";
    print "compressing $id\_$label\n";
    runCmd("gzip -c $fa > $d12/$id.1.fq.gz", 1); 
    runCmd("gzip -c $fb > $d12/$id.2.fq.gz", 1); 
  }
}
sub run_tophat {
  my ($f01, $d12, $d21) = @_;
  my $t = readTable(-in=>$f01, -header=>1);
  -d $d21 || make_path($d21);
  for my $i (0..$t->lastRow) {
#    next if $i >= 27;
    my ($sam, $id, $label, $note1, $note2) = $t->row($i);
    my $fa = "$d12/$id.1.fq.gz";
    my $fb = "$d12/$id.2.fq.gz";
    -s $fa || die "$fa is not there\n";
    -s $fb || die "$fb is not there\n";
    print "mapping $id\_$label to $sam\n";
    runCmd("tophat2 --num-threads 16 --mate-inner-dist 0 \\
      --min-intron-length 45 --max-intron-length 5000 \\
      --min-segment-intron 45 --max-segment-intron 5000 \\
      -o $d21/$id \$data/db/bowtie2/$sam $fa $fb", 1);
  }
}


