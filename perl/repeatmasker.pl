#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  repeatmasker.pl - run RepeatMasker and parse output

=head1 SYNOPSIS
  
  repeatmasker.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     RepeatMasker output file (*.out)
    -o (--out)    output prefix
    -s (--size)   chromosome size file

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

my ($fi, $fo, $fs) = ('') x 2;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "in|i=s"  => \$fi,
  "out|o=s" => \$fo,
  "size|s=s" => \$fs,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo || !$fs;

parse_rm($fi, "$fo.tbl");
runCmd("awk 'BEGIN{OFS=\"\\t\"} {print \$1, \$2-1, \$3, \$9 \" | \" (\$5)}' $fo.tbl > $fo.raw.bed");
runCmd("sortBed -i $fo.raw.bed > $fo.bed");
runCmd("bedToBigBed -tab -type=bed4 $fo.bed $fs $fo.bb");
runCmd("rm $fo.raw.bed");

sub parse_rm {
  my ($fi, $fo) = @_;
  open(my $fhi, "<$fi") or die "cannot read $fi\n";
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  while(<$fhi>) {
    chomp;
    my @ps = split " ";
    next if @ps == 0 || $ps[0] !~ /^[0-9e\.]+/i;
    my ($score, $div, $del, $ins, $qid, $qbeg, $qend, $qleft,
      $srd, $tid, $tfam, $tbeg, $tend, $tleft, $id) = @ps;
    if($srd eq "C") {
      $srd = "-";
      ($tbeg, $tleft) = ($tleft, $tbeg);
    }
    print $fho join("\t", $qid, $qbeg, $qend, $id, 
      $tid, $tbeg, $tend, $srd, $tfam, $score, $div, $del, $ins)."\n";
  }
  close $fhi;
  close $fho;
}


exit 0;
