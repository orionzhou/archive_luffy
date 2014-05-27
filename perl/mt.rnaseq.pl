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
#run_tophat("21.tbl", "12_fastq_gz", "22_tophat");
#merge_bam_genome("21.tbl", "22_tophat", "24_genome");

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
#    next if $i+1 < 37;
    my ($sam, $id, $genome) = $t->row($i);
    if(check_bam("$do/$id\_$genome/unmapped.bam")) {
      printf "%s|%s passed\n", $id, $genome;
      next;
    }
    printf "%s: working on %s\n", $id, $genome;
    my $fa = "$ds/$id.1.fq.gz";
    my $fb = "$ds/$id.2.fq.gz";
    -s $fa || die "$fa is not there\n";
    -s $fb || die "$fb is not there\n";
    runCmd("tophat2 --num-threads 16 --mate-inner-dist 0 \\
      --min-intron-length 45 --max-intron-length 5000 \\
      --min-segment-intron 45 --max-segment-intron 5000 \\
      -o $do/$id\_$genome \$data/db/bowtie2/$genome $fa $fb", 1);
  }
}
sub merge_bam_genome {
  my ($fi, $di, $do) = @_;
  my $t = readTable(-in=>$fi, -header=>1);
  -d $do || make_path($do);

  my $h;
  for my $i (0..$t->lastRow) {
    my ($sam, $id, $genome) = $t->row($i);
    $h->{$genome} ||= [];
    push @{$h->{$genome}}, $id;
  }

  for my $gm (keys(%$h)) {
    my @ids = @{$h->{$gm}};
    for my $id (@ids) {
      check_bam("$di/$id\_$gm/accepted_hits.bam") || 
        die "no $di/$id\_$gm/accepted_hits.bam\n";
    }
    my $str_rg = "\@RG\\tID:$gm\\tSM:$gm\\tLB:$gm\\tPL:ILLUMINA\\tPU:lane";
    my $str_in = join(" ", map {"$di/$_\_$gm/accepted_hits.bam"} @ids);
    runCmd("samtools cat $str_in -o $do/$gm.bam");
  }
}
  
    

