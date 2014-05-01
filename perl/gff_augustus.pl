#!/usr/bin/perl -w
use strict;
use FindBin;
use lib $FindBin::Bin;
use Getopt::Long;
use Pod::Usage;
use Common;

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
if ($fi eq "" || $fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, "<$fi") || die "Can't open file $fi for reading: $!\n";
}

if ($fo eq "" || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

print $fho "##gff-version 3\n";
while(<$fhi>) {
  chomp;
  next if !$_ || /^\#/;
  my @ps = split "\t";
  next if $ps[2] !~ /^(gene)|(transcript)|(CDS)$/;
  $ps[2] = "mRNA" if $ps[2] eq "transcript";
  print $fho join("\t", @ps)."\n";
}
close $fhi;
close $fho;


