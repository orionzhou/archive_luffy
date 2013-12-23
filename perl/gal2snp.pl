#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gal2snp.pl - Call SNPs from a Gal file

=head1 SYNOPSIS
  
  gal2snp.pl [-help] [-in input-file] [-qry qry-fasta] [-tgt tgt-fasta] [-out output-file]

  Options:
      -help   brief help message
      -in     input file
      -out    output file
      -qry    query-seq file 
      -tgt    target-seq file

=cut
  
#### END of POD documentation.
#-----------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Location;
use Gal;

my ($fi, $fo) = ('') x 2;
my ($fq, $ft) = ('') x 2; 
my ($fhi, $fho);
my $help_flag;

#----------------------------------- MAIN -----------------------------------#
GetOptions(
    "help|h"   => \$help_flag,
    "in|i=s"   => \$fi,
    "out|o=s"  => \$fo,
    "qry|q=s"  => \$fq,
    "tgt|t=s"  => \$ft,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;
pod2usage(2) if !$fq || !$ft;

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

while( <$fhi> ) {
    chomp;
    next if /(^id)|(^\#)|(^\s*$)/;
    my $ps = [ split "\t" ];
    next unless @$ps == 19;
    my ($id, $qId, $qBeg, $qEnd, $qSrd, $qSize, $tId, $tBeg, $tEnd, $tSrd, $tSize,
        $match, $misMatch, $baseN, $ident, $e, $score, $qLocS, $tLocS) = @$ps;
    
    my ($qPoss, $tPoss, $qNts, $tNts) = gal_get_mismatch($ps, $fq, $ft);
    for my $i (0..@$qPoss-1) {
        my ($qPos, $tPos, $qNt, $tNt) = ($qPoss->[$i], $tPoss->[$i], $qNts->[$i], $tNts->[$i]);
        print $fho join("\t", $qId, $qPos, $tId, $tPos, $qNt, $tNt)."\n";
    }
}
close $fhi;
close $fho;


__END__
