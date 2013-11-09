#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  winslide.pl - create sliding windows for given intervals

=head1 SYNOPSIS
  
  winslide.pl [-help] [-in input-file] [-step win-step] [-size win-size] [-out output-file]

  Options:
      -help   brief help message
      -in     input file
      -out    output file
      -step   sliding window step (default: 5)
      -size   sliding window size (default: 60)

=cut
  
#### END of POD documentation.
#-----------------------------------------------------------------------------

use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;

my ($fi, $fo) = ('') x 2;
my ($size, $step) = (60, 5);
my $help_flag;

#----------------------------------- MAIN -----------------------------------#
GetOptions(
    "help|h"  => \$help_flag,
    "in|i=s"  => \$fi,
    "out|o=s" => \$fo,
    "step|t=i" => \$step,
    "size|z=i" => \$size,
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

print $fho join("\t", qw/chr beg end/)."\n";
while(<$fhi>) {
    chomp;
    next if /(^#)|(^id\s)|(^chr\s)/;
    my @ps = split "\t";
    my ($chr, $beg, $end);
    next unless @ps >= 2;
    if(@ps >= 3) {
        ($chr, $beg, $end) = @ps;
    } else {
        ($chr, $beg, $end) = ($ps[0], 1, $ps[1]);
    }
    my $wins = sliding_windows($beg, $end, $step, $size);
    for (@$wins) {
        print $fho join("\t", $chr, $_->[0], $_->[1])."\n";
    }
}
close $fhi;
close $fho;

exit 0;
