#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  galCoord.pl - transform the coordinates of a Gal file

=head1 SYNOPSIS
  
  galCoord.pl [-help] [-in input-file] [-opt option] [-out output-file]

  Options:
      -help   brief help message
      -in     input file
      -out    output file
      -opt    option (coordq / coordt)

=head1 DESCRIPTION

  This program transforms the coordinates in an GAL file

=cut
  
#### END of POD documentation.
#-----------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;
use Seq;

my ($fi, $fo, $opt) = ('') x 3;
my ($fhi, $fho);
my $help_flag;

#----------------------------------- MAIN -----------------------------------#
GetOptions(
    "help|h"   => \$help_flag,
    "in|i=s"   => \$fi,
    "out|o=s"  => \$fo,
    "opt|p=s"  => \$opt,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo || !$opt;

if ($fi eq "stdin") {
    $fhi = \*STDIN;
} else {
    open ($fhi, $fi) || die "Can't open file $fi: $!\n";
}

if ($fo eq "stdout") {
    $fho = \*STDOUT;
} else {
    open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

print $fho join("\t", qw/id qId qBeg qEnd qSrd qLen tId tBeg tEnd tSrd tLen
    match misMatch baseN ident e score/)."\n";

while(<$fhi>) {
    chomp;
    next if /(^id)|(^\#)|(^\s*$)/;
    my ($id, $qId, $qBeg, $qEnd, $qSrd, $qLen, $tId, $tBeg, $tEnd, $tSrd, $tLen,
        $match, $misMatch, $baseN, $ident, $e, $score) = split "\t";
    if($opt =~ /^coordq$/i) {
        if($qId =~ /^(\w+)\-([0-9e\+]+)\-([0-9e\+]+)$/) {
            my ($qi, $beg, $end) = ($1, $2, $3);
            $qId = $qi;
            $qBeg = $beg + $qBeg - 1;
            $qEnd = $beg + $qEnd - 1;
        }
    } elsif($opt =~ /^coordt$/i) {
        if($tId =~ /^(\w+)\-([0-9e\+]+)\-([0-9e\+]+)$/) {
            my ($ti, $beg, $end) = ($1, $2, $3);
            $tId = $ti;
            $tBeg = $beg + $tBeg - 1;
            $tEnd = $beg + $tEnd - 1;
        }
    } else {
        die "unknown opt: $opt\n";
    }
    print $fho join("\t", $id, $qId, $qBeg, $qEnd, $qSrd, $qLen,
        $tId, $tBeg, $tEnd, $tSrd, $tLen,
        $match, $misMatch, $baseN, $ident, $e, $score)."\n";
}
close $fhi;
close $fho;


__END__
