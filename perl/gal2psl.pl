#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gal2psl.pl - convert GAL file to PSL format

=head1 SYNOPSIS
  
  gal2psl.pl [-help] [-in input-file] [-out output-file]

  Options:
      -help   brief help message
      -in     input file
      -out    output file

=cut
  
#### END of POD documentation.
#-----------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;
use Location;
use Gal;
use List::Util qw/min max sum/;

my ($fi, $fo) = ('', '');
my ($fhi, $fho);
my $help_flag;

#----------------------------------- MAIN -----------------------------------#
GetOptions(
    "help|h"   => \$help_flag,
    "in|i=s"   => \$fi,
    "out|o=s"  => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;

if ($fi eq '' || $fi eq "stdin" || $fi eq "-") {
    $fhi = \*STDIN;
} else {
    open ($fhi, $fi) || die "Can't open file $fi: $!\n";
}

if ($fo eq '' || $fo eq "stdout" || $fo eq "-") {
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
    my $srd = ($qSrd eq $tSrd) ? "+" : "-";

    my ($qLoc, $tLoc) = (locStr2Ary($qLocS), locStr2Ary($tLocS));
    @$qLoc == @$tLoc || die "unequal pieces\n";
    my $nBlock = @$qLoc;
    
    my (@blockSizes, @qBegs, @tBegs);
    my (@qIns, @tIns);
    my ($rqe_p, $rte_p);
    for my $i (0..$nBlock-1) {
        my ($rqb, $rqe) = @{$qLoc->[$i]};
        my ($rtb, $rte) = @{$tLoc->[$i]};
        my ($len, $len2) = ($rqe-$rqb+1, $rte-$rtb+1);
        die "block size unequal: $qId-$tId $rqb-$rqe : $rtb-$rte\n" if $len != $len2;
        my $tb = $tBeg + $rtb - 1;
        my $qb = $srd eq "-" ? $qSize-$qEnd+1 + $rqb-1 : $qBeg + $rqb - 1;
        
        push @blockSizes, $len;
        push @tBegs, $tb-1;
        push @qBegs, $qb-1;
        if($i > 0) {
            my $tIns = $rtb - $rte_p - 1;
            my $qIns = $rqb - $rqe_p - 1;
            push @tIns, $tIns if $tIns > 0;
            push @qIns, $qIns if $qIns > 0;
        }
        ($rqe_p, $rte_p) = ($rqe, $rte);
    }
    my $repMatch = 0;
    my ($qNumIns, $tNumIns) = (scalar(@qIns), scalar(@tIns));
    my ($qBaseIns, $tBaseIns) = (0, 0);
    $qBaseIns = sum(@qIns) if $qNumIns > 0;
    $tBaseIns = sum(@tIns) if $tNumIns > 0;
    my $blockSizes = join(",", @blockSizes).",";
    my $qBegs = join(",", @qBegs).",";
    my $tBegs = join(",", @tBegs).",";
    print $fho join("\t", $match, $misMatch, $repMatch, $baseN, 
        $qNumIns, $qBaseIns, $tNumIns, $tBaseIns, $srd, 
        $qId, $qSize, $qBeg-1, $qEnd, $tId, $tSize, $tBeg-1, $tEnd, 
            $nBlock, $blockSizes, $qBegs, $tBegs)."\n";
}
close $fhi;
close $fho;

__END__
