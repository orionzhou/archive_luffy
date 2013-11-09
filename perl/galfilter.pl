#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  galfix.pl - filter records in a GAL file

=head1 SYNOPSIS
  
  galfix.pl [-help] [-in input-file] [-out output-file]

  Options:
      -help   brief help message
      -in     input file
      -out    output file
      -len    minimum aln length (default: 1) 

=cut
  
#### END of POD documentation.
#-----------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Location;

my ($fi, $fo) = ('') x 2;
my $len = 1;
my ($fhi, $fho);
my $help_flag;

#----------------------------------- MAIN -----------------------------------#
GetOptions(
    "help|h"   => \$help_flag,
    "in|i=s"   => \$fi,
    "out|o=s"  => \$fo,
    "len|l=i"  => \$len,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

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

print $fho join("\t", qw/id qId qBeg qEnd qSrd qSize tId tBeg tEnd tSrd tSize
    match misMatch baseN ident e score qLoc tLoc/)."\n";

my $cnt = 0;
while( <$fhi> ) {
    chomp;
    next if /(^id)|(^\#)|(^\s*$)/;
    my $ps = [ split "\t" ];
    next unless @$ps == 19;
    my ($id, $qId, $qBeg, $qEnd, $qSrd, $qSize, $tId, $tBeg, $tEnd, $tSrd, $tSize,
        $match, $misMatch, $baseN, $ident, $e, $score, $qLocS, $tLocS) = @$ps;
    my $qLoc = locStr2Ary($qLocS);
    my $qLen = locAryLen($qLoc);
    next if $qLen < $len;

    print $fho join("\t", @$ps)."\n";
    $cnt ++;
}
print STDERR "$cnt rows passed filter\n";
close $fhi;
close $fho;


__END__
