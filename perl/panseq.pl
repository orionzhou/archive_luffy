#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  panseq.pl - construct pan-genome fasta

=head1 SYNOPSIS
  
  panseq.pl [-help]

  Options:
    -h (--help)   brief help message
    -s (--stat)   print statistics

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin";
use Bio::SeqIO;
use Bio::AlignIO;
use Data::Dumper;
use Common;
use File::Path qw/make_path remove_tree/;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;

my $help_flag;
my $stat_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "stat|s"  => \$stat_flag,
) or pod2usage(2);
pod2usage(1) if $help_flag;

my $dir = "$ENV{'misc3'}/panseq";
chdir $dir || die "cannot chdir to $dir\n";

my @orgs = qw/HM004 HM010 HM018 HM022 HM034 HM050/;# HM056 HM058 HM060 
#  HM095 HM125 HM129 HM185 HM324 HM340/;

if($stat_flag) {
  exit;
}

for my $org (@orgs) {
  runCmd("ln -sf $ENV{'misc3'}/$org\_HM101/41_novseq/21.fas $org.fas");
}
my $fstr = join(" ", map {"$_.fas"} @orgs);
runCmd("mugsy -p out --directory $dir $fstr");

maf2tbl('out.maf', 'out.tbl');
sub maf2tbl {
  my ($fi, $fo) = @_;
  my $ai = Bio::AlignIO->new(-file => $fi, -format => 'maf');
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho join("\t", qw/cid org chr beg end srd/)."\n";
  my $i = 1;
  while(my $aln = $ai->next_aln()) {
    my $n = $aln->num_sequences;
    for my $seq ($aln->each_seq) {  
      my ($id, $beg, $end, $srd) = ($seq->id, $seq->start, $seq->end,
        $seq->strand);
      my ($org, $nid) = split(/\./, $id);
      $srd = $srd == -1 ? "-" : "+";
      print $fho join("\t", $i, $org, $nid, $beg, $end, $srd)."\n";
    }
    $i ++;
  }
  close $fho;
}

