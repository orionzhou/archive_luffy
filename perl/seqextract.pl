#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  seqextract.pl - extract one or more sequence records from an input fasta file

=head1 SYNOPSIS
  
  seqextract.pl [-help] [-out output-file] <input-file> <IDs>

  Options:
      -help   brief help message
      -out    output file, instead of stdout

=head1 DESCRIPTION

  This program extracts sequences from the input fasta file

=head1 OPTIONS

=over 6
  
=item B<-help>
  
  Print a usage summary.

=item B<input-file>
  
  Needs to be fasta format

=item B<IDs>
  
  Fasta IDs

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
my $help_flag;

#----------------------------------- MAIN -----------------------------------#
GetOptions(
    "help|h"  => \$help_flag,
    "out|o=s" => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;

my @ids;
($fi, @ids)= @ARGV;
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
my $seqHO = Bio::SeqIO->new(-fh=>$fho, -format=>'fasta');
my $h = { map {$_=>1} @ids };
while(my $seqO = $seqHI->next_seq()) {
    my ($id) = ($seqO->id);
    $seqHO->write_seq($seqO) if exists $h->{$id};
}
$seqHI->close();
$seqHO->close();

exit 0;



