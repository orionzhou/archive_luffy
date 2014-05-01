#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  seqgap.pl - report gap ('N's) locations in an input fasta file

=head1 SYNOPSIS
  
  seqgap.pl [-help] [-min minimum-len] [-in input-file] [-out output-file] 

  Options:
    -h (--help)   brief help message
    -i (--in)     input (Fasta) file
    -o (--out)    output (BED) file
    -m (--min)    minimum length of gaps to report (default: 10)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::Seq;
use Bio::SeqIO;

my ($fi, $fo) = ('') x 2;
my $fhi;
my $fho;
my $len_min = 1;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "in|i=s"  => \$fi,
  "out|o=s" => \$fo,
  "min|m=i" => \$len_min,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

if ($fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "Can't open file $fi: $!\n";
}

if ($fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

my $seqHI = Bio::SeqIO->new(-fh=>$fhi, -format=>'fasta');
open($fho, ">$fo") or die "cannot open $fo for writing\n";
#print $fho join("\t", qw/id beg end len/)."\n";

while(my $seqO = $seqHI->next_seq()) {
  my ($id, $seq, $seqlen) = ($seqO->id, $seqO->seq, $seqO->length);
  while($seq =~ /N+/ig) {
    my ($beg, $end) = ($-[0]+1, $+[0]);
    my $len = $end - $beg + 1;
    next if $len < $len_min;
    print $fho join("\t", $id, $beg - 1, $end)."\n";
  }
}
$seqHI->close();
close $fho;

exit 0;
