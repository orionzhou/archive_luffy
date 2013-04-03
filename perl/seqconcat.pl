#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  seqconcat.pl - concatenate all sequences from input file to one sequence

=head1 SYNOPSIS
  
  seqconcat.pl [-help] [-out output-file] [-id seq-id] [-gap length of gap] <input-file>

  Options:
      -help   brief help message
      -id     output sequene ID (default: "chrU")
      -gap    gap length (default: 1000)
      -out    output file, instead of stdout

=head1 DESCRIPTION

  This program concatenates all sequences from an input file to one sequence

=head1 OPTIONS

=over 6
  
=item B<-help>
  
  Print a usage summary.

=item B<input-file>
  
  Needs to be fasta format

=item B<output-file>

  To write to stdout, the user could either specify 'stdout' or simply leave this
  augument empty.

=back
  
=head1 BUGS
  
=head1 REFERENCES
  
=head1 VERSION
  
  0.1
  
=cut
  
#### END of POD documentation.
#-----------------------------------------------------------------------------



use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::Seq;
use Bio::SeqIO;

my $fi = '';
my $fo = '';
my $gapLen = 1000;
my $seqid = "chrU";
my $fhi;
my $fho;
my $help_flag;

#----------------------------------- MAIN -----------------------------------#
GetOptions(
    "help|h"  => \$help_flag,
    "in|i=s"  => \$fi,
    "out|o=s" => \$fo,
    "gap|g=i" => \$gapLen,
    "id|d=s"  => \$seqid,
) or pod2usage(2);
pod2usage(1) if $help_flag;

($fi)= @ARGV;
if(!$fi) {
    pod2usage(2);
} elsif ($fi eq '-' || $fi eq "stdin") {
    $fhi = \*STDIN;
} else {
    open ($fhi, "<$fi") || die "Can't open file $fi for reading: $!\n";
}

if(!$fo || $fo eq "stdout") {
    $fho = \*STDOUT;
} else {
    open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

my $seqHI = Bio::SeqIO->new(-fh=>$fhi, -format=>'fasta');
my ($endPrev, $cnt) = (0, 0);
my $seqStr = "";
while(my $seqO = $seqHI->next_seq()) {
    my ($seqLen, $seq) = ($seqO->length, $seqO->seq);
    if($endPrev > 0) {
        $seqStr .= 'N' x $gapLen;
        $endPrev += $gapLen;
    }
    $seqStr .= $seq;
    $endPrev += $seqLen;
}
$seqHI->close();

my $seqHO = Bio::SeqIO->new(-fh=>$fho, -format=>'fasta');
$seqHO->write_seq(Bio::Seq->new(-id=>$seqid, -seq=>$seqStr));
$seqHO->close();

exit 0;


