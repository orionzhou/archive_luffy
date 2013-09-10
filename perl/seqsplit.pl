#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  seqsplit.pl - split a fasta file into N files of equal sizes

=head1 SYNOPSIS
  
  seqlen.pl [-help] [-n number-of-output-files] [-in input-file] [-out output-file-prefix]

  Options:
      -help   brief help message
      -in     input file
      -out    output file prefix
      -n      number of output files [default=1]

=head1 DESCRIPTION

  This program split a multi-fasta file into N output files (each of equal size)

=cut
  
#### END of POD documentation.
#-----------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::Seq;
use Bio::SeqIO;

my ($fi, $fo) = ('') x 2;
my $n = 1;
my $help_flag;

#----------------------------------- MAIN -----------------------------------#
GetOptions(
    "help|h"  => \$help_flag,
    "in|i=s"  => \$fi,
    "out|o=s" => \$fo,
    "n=i" => \$n,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

open(my $fhi, "<$fi") || die "Can't open file $fi for reading";
my $seqHOs = [];
for my $i (0..$n-1) {
    my $fon = sprintf "%s\_%03d.fa", $fo, $i;
    my $seqHO = Bio::SeqIO->new(-file=>">$fon", -format=>'fasta');
    push @$seqHOs, $seqHO;
}

my $seqHI = Bio::SeqIO->new(-fh=>$fhi, -format=>'fasta');
my $cnt = 0;
while(my $seqO = $seqHI->next_seq()) {
    my $i = $cnt % $n;
    $seqHOs->[$i]->write_seq($seqO);
    $cnt ++;
}
$seqHI->close();

for my $i (0..$n-1) {
    $seqHOs->[$i]->close();
}

exit 0;
