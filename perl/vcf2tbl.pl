#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  vcf2tbl.pl - convert VCF file to TBL file

=head1 SYNOPSIS
  
  vcf2tbl.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input (VCF) file
    -o (--out)    output (TBL) file

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;

my ($fi, $fo) = ('') x 2;
my ($fhi, $fho);
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

if ($fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, "<$fi") || die "cannot read $fi\n";
}

if ($fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

while( <$fhi> ) {
  chomp;
  next if /(^\#)|(^\s*$)/s;
  my ($chr, $pos, $id, $ref, $alt, $qual, $fil, $info, $fmt, @sams) = 
    split "\t";
  my @alts = split(",", $alt);
  @alts == 1 || die "$chr:$pos $ref-$alt >1 alts\n";
  print $fho join("\t", $chr, $pos, $ref, $alt, $qual)."\n";
}
close $fhi;
close $fho;


__END__
