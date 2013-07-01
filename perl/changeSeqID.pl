#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  changeSeqID.pl - replace IDs in a fasta file

=head1 SYNOPSIS
  
  seqlen.pl [-help] [-opt convert-option] [-in input-file] [-out output-file]

  Options:
      -help   brief help message
      -in     input file
      -out    output file
      -opt    convert option (default=1)

=head1 DESCRIPTION

  This program convert sequence identifiers in a Fasta file to standard IDs

=cut
  
#### END of POD documentation.
#-----------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::SeqIO;

my ($fi, $fo) = ('') x 2;
my $opt = 1;
my $help_flag;

#----------------------------------- MAIN -----------------------------------#
GetOptions(
    "help|h"  => \$help_flag,
    "in|i=s"  => \$fi,
    "out|o=s" => \$fo,
    "opt|p=i" => \$opt,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

open(my $fhi, "<$fi") || die "Can't open file $fi for reading";
open(my $fho, ">$fo") || die "Can't open file $fo for writing";

my $seqHI = Bio::SeqIO->new(-fh=>$fhi, -format=>'fasta');
my $seqHO = Bio::SeqIO->new(-fh=>$fho, -format=>'fasta');
while(my $seqO = $seqHI->next_seq()) {
    my $id = $seqO->id;
    if($opt == 1) {
        $id = convert_id($id);
    } elsif($opt == 7) { # Zmays
        $id = convert_id_zmays($id);
    }
    $seqHO->write_seq(Bio::Seq->new(-id=>$id, -seq=>$seqO->seq()));
}
$seqHI->close();
$seqHO->close();

sub convert_id {
    my ($id) = @_;
    if($id =~ /chr(\w+)/i) {
        $id = "chr$1";
    } elsif($id =~ /\|?([A-Z]{2}\d{6}\.[0-9DF]{1,2})\|?/) {
        $id = $1;
    } elsif($id =~ /([A-Z]{2}\d{6})([^\.]|$)/) {
        $id = $1;
    }
    return $id;
}
sub convert_id_zmays {
    my ($id) = @_;
    if($id =~ /chromosome (\w+)$/) {
        $id = $1;
        if($id =~ /^\d+$/) {
            $id = "chr$id";
        } elsif($id eq "UNKNOWN") {
            $id = "chrU";
        } elsif($id eq "mitochondrion") {
            $id = "Mt";
        } elsif($id eq "chloroplast") {
            $id = "Pt"
        } else {
            die "unknonw chr ID: $id\n";
        }
    } else {
        die "unknonw chr ID: $id\n";
    }
    return $id;
}
exit 0;

