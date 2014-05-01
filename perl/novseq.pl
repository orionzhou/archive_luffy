#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  novseq.pl - identify and characterize novel sequences

=head1 SYNOPSIS
  
  novseq.pl [-help] [-qry query-genome] [-tgt target-genome]

  Options:
    -h (--help)   brief help message
    -q (--qry)    query genome
    -t (--tgt)    target genome

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin";
use Common;
use Data::Dumper;
use File::Path qw/make_path remove_tree/;
use File::Basename;
use List::Util qw/min max sum/;

my ($qry, $tgt) = ('HM056', 'HM101');
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "qry|q=s" => \$qry,
  "tgt|t=s" => \$tgt,
) or pod2usage(2);
pod2usage(1) if $help_flag;

my $data = '/home/youngn/zhoup/Data';
my $qry_fas = "$data/genome/$qry/11_genome.fa";
my $tgt_fas = "$data/genome/$tgt/11_genome.fa";
my $qry_2bit = "$data/db/blat/$qry.2bit";
my $tgt_2bit = "$data/db/blat/$tgt.2bit";
my $qry_size = "$data/genome/$qry/15.sizes";
my $tgt_size = "$data/genome/$tgt/15.sizes";
my $qry_size_bed = "$data/genome/$qry/15.bed";
my $tgt_size_bed = "$data/genome/$tgt/15.bed";
my $qry_gap = "$data/genome/$qry/16_gap.bed";
my $tgt_gap = "$data/genome/$tgt/16_gap.bed";

my $dir = "$data/misc3/$qry\_$tgt/41_novseq";
-d $dir || make_path($dir);
chdir $dir || die "cannot chdir to $dir\n";

pipe_novseq1();

sub pipe_novseq1 {
  runCmd("gal2gax.pl -i ../23_blat/41.3.gal -o 41.3.gax");
  runCmd("gax2bed.pl -i 41.3.gax -p tgt -o - | sortBed -i stdin | \\
    mergeBed -i stdin > 41.3.bed");
  runCmd("subtractBed -a $qry_size_bed -b $qry_gap | \\
    subtractBed -a stdin -b 41.3.bed | bedfilter.pl -l 50 -o 01.bed");
  runCmd("rm 41.3.gax 41.3.bed");
  runCmd("seqret.pl -d $qry_fas -b 01.bed -o 01.fas");

  runCmd("dustmasker -in 01.fas -outfmt interval -out 03.dust.txt");
  dust2bed("03.dust.txt", "03.dust.bed");  
  runCmd("bedfilter.pl -i 03.dust.bed -l 10 -o 04.dust.filter.bed");
  runCmd("subtractBed -a 01.bed -b 04.dust.filter.bed | \\
    bedfilter.pl -l 50 -o 11.bed");
  runCmd("seqret.pl -d $qry_fas -b 11.bed -o 11.fas");
  runCmd("usearch.pl -i 11.fas -o 21");
}

sub dust2bed {
  my ($fi, $fo) = @_;
  open(my $fhi, "<$fi") or die "cannot read $fi\n";
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  my ($id, $beg, $end);
  while(<$fhi>) {
    chomp;
    if(/^\>(\S+)$/) {
      ($id, $beg, $end) = split("-", $1);
    } else {
      my ($rb, $re) = split "-";
      print $fho join("\t", $id, $beg + $rb, $beg + $re + 1)."\n";
    }
  }
  close $fhi;
  close $fho;
}
sub cdhit2bed {
  my ($fi, $fo) = @_;
  open(my $fhi, "<$fi") or die "cannot read $fi\n";
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  my ($h, $hl);
  while(<$fhi>) {
    chomp;
    if(/^\>(\S+)$/) {
      my ($id, $beg, $end) = split("-", $1);
      print $fho join("\t", $id, $beg - 1, $end)."\n";
    }
  }
  close $fhi;
  close $fho;
}




