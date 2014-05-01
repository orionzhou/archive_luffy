#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  equalindel.pl

=head1 SYNOPSIS
  
  equalindel.pl [-help] [-in input-file] [-out output-file]

  Options:
      -help   brief help message
      -in     input file
      -out    output file

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::Seq;
use Bio::SeqIO;

my ($fi, $fo) = ('') x 2;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "in|i=s"  => \$fi,
  "out|o=s" => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

my ($fhi, $fho);
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
my $seqHO = Bio::SeqIO->new(-fh=>$fho, -format=>'fasta');
while(my $seqO = $seqHI->next_seq()) {
  my ($id, $len) = ($seqO->id, $seqO->length);
  my $newid = sprintf("x:x:%s:%d:%d:x", $id, 1, $len);
  $seqO->id($newid);
  $seqHO->write_seq($seqO);
}
$seqHI->close();
$seqHO->close();

exit 0;
