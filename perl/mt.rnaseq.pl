#!/usr/bin/perl -w
use strict;
use FindBin;
use lib $FindBin::Bin;
use File::Path qw/make_path remove_tree/;
use Data::Dumper;
use List::Util qw/min max sum/;
use Common;
use Bam;

my $dir = "/home/youngn/zhoup/Data/misc2/rnaseq/mt";
-d $dir || make_path($dir);
chdir $dir || die "cannot chdir $dir\n";

#gzip_fastq("01.tbl", "11_fastq", "12_fastq_gz");
#run_tophat("01.tbl", "12_fastq_gz", "21_tophat");
#merge_bam_treatment("01.tbl", "21_tophat", "23_treatment");
#merge_bam_sample("01.tbl", "21_tophat", "24_sample");

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
  my ($fi, $ds, $do) = @_;
  my $t = readTable(-in=>$fi, -header=>1);
  -d $do || make_path($do);
  for my $i (0..$t->lastRow) {
#    next if $i >= 27;
    my ($sam, $id, $label, $note1, $note2) = $t->row($i);
    if(check_bam("$do/$id/unmapped.bam")) {
      printf "%s passed\n", $id;
      next;
    }
    printf "%s: working on %s\n", $id, $sam;
    my $fa = "$ds/$id.1.fq.gz";
    my $fb = "$ds/$id.2.fq.gz";
    -s $fa || die "$fa is not there\n";
    -s $fb || die "$fb is not there\n";
    runCmd("tophat2 --num-threads 16 --mate-inner-dist 0 \\
      --min-intron-length 45 --max-intron-length 5000 \\
      --min-segment-intron 45 --max-segment-intron 5000 \\
      -o $do/$id \$data/db/bowtie2/$sam $fa $fb", 1);
  }
}
sub merge_bam_sample {
  my ($fi, $di, $do) = @_;
  my $t = readTable(-in=>$fi, -header=>1);
  -d $do || make_path($do);

  my $h;
  for my $i (0..$t->lastRow) {
    my ($sam, $id, $label, $note1, $note2) = $t->row($i);
    $h->{$sam} ||= [];
    push @{$h->{$sam}}, $id;
  }

  for my $sm (keys(%$h)) {
    my @ids = @{$h->{$sm}};
    for my $id (@ids) {
      check_bam("$di/$id/accepted_hits.bam") || 
        die "no $di/$id/accepted_hits.bam\n";
    }
    my $str_rg = "\@RG\\tID:$sm\\tSM:$sm\\tLB:$sm\\tPL:ILLUMINA\\tPU:lane";
    my $str_in = join(" ", map {"$di/$_/accepted_hits.bam"} @ids);
    runCmd("samtools cat $str_in -o $do/$sm.bam");
  }
}
  
    

