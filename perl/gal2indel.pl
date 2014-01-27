#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gal2indel.pl - Call InDels from a Gal file

=head1 SYNOPSIS
  
  gal2indel.pl [-help] [-in input-Gal] [-out output-file]

  Options:
      -h (--help)   brief help message
      -i (--in)     input Gal
      -o (--out)    output indel

=cut
  
#### END of POD documentation.
#-----------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Location;
use Common;
use Gal;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;

my ($fi, $fo) = ('') x 2;
my ($fn) = (''); 
my ($fhi, $fho);
my $help_flag;

#----------------------------------- MAIN -----------------------------------#
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
    open ($fhi, $fi) || die "Can't open file $fi: $!\n";
}

if ($fo eq "stdout" || $fo eq "-") {
    $fho = \*STDOUT;
} else {
    open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

print $fho join("\t", qw/qid qbeg qend qins tid tbeg tend tins id/)."\n";
while( <$fhi> ) {
    chomp;
    next if /(^id)|(^\#)|(^\s*$)/;
    my $ps = [ split "\t" ];
    next unless @$ps == 19;
    my ($id, $qId, $qBeg, $qEnd, $qSrd, $qSize, $tId, $tBeg, $tEnd, $tSrd, $tSize,
        $match, $misMatch, $baseN, $ident, $e, $score, $qLocS, $tLocS) = @$ps;

    my ($rqLoc, $rtLoc) = (locStr2Ary($qLocS), locStr2Ary($tLocS));
    $rqLoc = [ sort {$a->[0] <=> $b->[0]} @$rqLoc ];
    $rtLoc = [ sort {$a->[0] <=> $b->[0]} @$rtLoc ];
    @$rqLoc == @$rtLoc || die "unequal pieces\n";
    my $nBlock = @$rqLoc;
    my @lens = map {$_->[1] - $_->[0] + 1} @$rqLoc;
   
    $nBlock > 1 || next;
    for my $i (1..$nBlock-1) {  
        my ($rqb, $rqe) = ($rqLoc->[$i-1]->[1], $rqLoc->[$i]->[0]);
        my ($rtb, $rte) = ($rtLoc->[$i-1]->[1], $rtLoc->[$i]->[0]);
        my $qins = $rqe - $rqb - 1;
        my $tins = $rte - $rtb - 1;
        my ($qb, $qe) = $qSrd eq "-" ? ($qEnd - $rqe + 1, $qEnd - $rqb + 1) : 
            ($qBeg + $rqb - 1, $qBeg + $rqe - 1);
        my ($tb, $te) = $tSrd eq "-" ? ($tEnd - $rte + 1, $tEnd - $rtb + 1) : 
            ($tBeg + $rtb - 1, $tBeg + $rte - 1);
        print $fho join("\t", $qId, $qb, $qe, $qins, $tId, $tb, $te, $tins, $id)."\n";
    }
}
close $fhi;
close $fho;


__END__
