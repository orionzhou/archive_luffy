#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  seqgc.pl - calculate GC content in a sliding window approach for a fasta file

=head1 SYNOPSIS
  
  seqlen.pl [-help] [-in input-tbl] [-seq seq-fasta] [-out output-tbl]

  Options:
      -help   brief help message
      -in     input file
      -out    output file
      -seq    sequence file (fasta)

=cut
  
#### END of POD documentation.
#-----------------------------------------------------------------------------

use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Seq;

my ($fi, $fo, $fs) = ('') x 3;
my $help_flag;

#----------------------------------- MAIN -----------------------------------#
GetOptions(
    "help|h"  => \$help_flag,
    "in|i=s"  => \$fi,
    "out|o=s" => \$fo,
    "seq|s=s" => \$fs,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo || !$fs;

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

print $fho join("\t", qw/chr beg end gc/)."\n";
while(<$fhi>) {
    chomp;
    next if /(^#)|(^id\s)|(^chr\s)/;
    my @ps = split "\t";
    next unless @ps >= 3;
    my ($chr, $beg, $end) = @ps;
    my $seq = seqRet([[$beg, $end]], $chr, "+", $fs);
    print $fho join("\t", $chr, $beg, $end, calc_gc($seq))."\n";
}
close $fhi;
close $fho;

sub calc_gc {
    my ($str) = @_;
    my ($cntGC, $cntN) = (0, 0);
    while($str =~ /([GCN])/ig) {
        if($1 eq "N") {
            $cntN ++;
        } else {
            $cntGC ++;
        }
    }
    my $len = length($str) - $cntN;
    return 0 if $len == 0;
    return sprintf "%.03f", $cntGC / $len;
}



exit 0;
