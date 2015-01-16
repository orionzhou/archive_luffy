#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  seqlen.pl - report sequence lengths in an input sequence file 

=head1 SYNOPSIS
  
  seqlen.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input (fasta) file (default: stdin)
    -o (--out)    output (tabular) file (default: stdout)

=head1 DESCRIPTION

  This program reports lengths of each sequence record in the input file

=cut
  
#### END of POD documentation.
#------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::Seq;
use Bio::SeqIO;

my ($fi, $fo) = ('') x 2;
my $help_flag;

#------------------------------ MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "in|i=s"  => \$fi,
  "out|o=s" => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;

my ($fhi, $fho);
if ($fi eq "" || $fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, "<$fi") || die "cannot read $fi\n";
}

if ($fo eq "" || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

my $seqHI = Bio::SeqIO->new(-fh => $fhi, -format => 'fasta');
while(my $seqO = $seqHI->next_seq()) {
  my ($id, $len) = ($seqO->id, $seqO->length);
  print $fho join("\t", $id, $len)."\n";
}
$seqHI->close();
close $fho;

exit 0;
