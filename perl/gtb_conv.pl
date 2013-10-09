#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gtb_conv.pl - convert a Gtb file using specified format 

=head1 SYNOPSIS
  
  gtb_conv.pl [-help] [-in input-file] [-fmt output-format] [-out output-file]

  Options:
      -help   brief help message
      -in     input file
      -out    output file
      -fmt    output format (default: gff)
      -ref    reference sequence file (optional)
      -size   chromosome size file (optional)

=head1 DESCRIPTION

  This program converts an input Gtb file to an output file of specified format

=cut
  
#### END of POD documentation.
#-----------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;
use Gtb;

my ($fi, $fo, $fr, $fs) = ('') x 4;
my $fmt = "gff";
my $help_flag;

#----------------------------------- MAIN -----------------------------------#
GetOptions(
    "help|h"   => \$help_flag,
    "in|i=s"   => \$fi,
    "out|o=s"  => \$fo,
    "fmt|f=s"  => \$fmt,
    "ref|r=s"  => \$fr,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

if($fmt =~ /^gff3?$/i) {
    gtb2Gff($fi, $fo);
} elsif($fmt =~ /^bed$/i) {
    gtb2Bed($fi, $fo);
} elsif($fmt =~ /^bigbed$/i) {
    -s $fs || die "$fs is not there\n";
    gtb2BigBed($fi, $fo, $fs);
} elsif($fmt =~ /^tbl$/i) {
    gtb2Tbl($fi, $fo);
} elsif($fmt =~ /^seq1$/i) {
    gtb2Seq(-in=>$fi, -out=>$fo, -seq=>$fr, -opt=>1);
} elsif($fmt =~ /^seq2$/i) {
    gtb2Seq(-in=>$fi, -out=>$fo, -seq=>$fr, -opt=>2);
} else {
    die "unsupported format: $fmt\n";
}


__END__
