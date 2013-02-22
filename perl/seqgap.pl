#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  seqgap.pl - report gap ('N's) locations in an input fasta file

=head1 SYNOPSIS
  
  seqgap.pl [-help] [-min minimum-length] [-out output-file] <input-file>

  Options:
      -help   brief help message
      -out    output file, instead of stdout
      -min    minimum length of gaps to report (default: 1)

=head1 DESCRIPTION

  This program report locations of all gaps ('N's) from the input fasta file

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
my $fhi;
my $fho;
my $len_min = 1;
my $help_flag;

#----------------------------------- MAIN -----------------------------------#
GetOptions(
    "help|h"  => \$help_flag,
    "out|o=s" => \$fo,
    "min|i=i" => \$len_min,
) or pod2usage(2);
pod2usage(1) if $help_flag;

my @ids;
($fi)= @ARGV;
if(!$fi) {
    pod2usage(2);
} elsif ($fi eq '-' || $fi eq "stdin") {
    $fhi = \*STDIN;
} else {
    open ($fhi, $fi) || die "Can't open file $fi: $!\n";
}

if(!$fo || $fo eq "stdout") {
    $fho = \*STDOUT;
} else {
    open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

my $seqHI = Bio::SeqIO->new(-fh=>$fhi, -format=>'fasta');
open($fho, ">$fo") or die "cannot open $fo for writing\n";
print $fho join("\t", qw/id beg end len/)."\n";

while(my $seqO = $seqHI->next_seq()) {
    my ($id, $seq, $seqlen) = ($seqO->id, $seqO->seq, $seqO->length);
    while($seq =~ /N+/ig) {
        my ($beg, $end) = ($-[0]+1, $+[0]);
        my $len = $end - $beg + 1;
        next if $len < $len_min;
        print $fho join("\t", $id, $beg, $end, $len)."\n";
    }
}
$seqHI->close();
close $fho;

exit 0;
