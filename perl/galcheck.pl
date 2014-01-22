#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  galcheck.pl - Check a GAL file

=head1 SYNOPSIS
  
  galcheck.pl [-help] [-in input-file]

  Options:
      -help   brief help message
      -in     input file

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
my ($fhi, $fho);
my $help_flag;

#----------------------------------- MAIN -----------------------------------#
GetOptions(
    "help|h"   => \$help_flag,
    "in|i=s"   => \$fi,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi;

if ($fi eq "stdin" || $fi eq "-") {
    $fhi = \*STDIN;
} else {
    open ($fhi, $fi) || die "Can't open file $fi: $!\n";
}

while( <$fhi> ) {
    chomp;
    next if /(^id)|(^\#)|(^\s*$)/;
    my $ps = [ split "\t" ];
    next unless @$ps == 19;
    
    my ($id, $qId, $qBeg, $qEnd, $qSrd, $qSize, $tId, $tBeg, $tEnd, $tSrd, $tSize,
        $match, $misMatch, $baseN, $ident, $e, $score, $qLocS, $tLocS) = @$ps;
    my ($qLoc, $tLoc) = (locStr2Ary($qLocS), locStr2Ary($tLocS));
    my $srd = ($qSrd eq $tSrd) ? "+" : "-";

    my ($rqb, $rqe) = ($qLoc->[0]->[0], $qLoc->[-1]->[1]);
    my ($rtb, $rte) = ($tLoc->[0]->[0], $tLoc->[-1]->[1]);
    print "qLen err: $id $qId($qBeg-$qEnd): $rqb-$rqe\n" unless $qEnd-$qBeg==$rqe-$rqb;
    print "tLen err: $id $tId($tBeg-$tEnd): $rtb-$rte\n" unless $tEnd-$tBeg==$rte-$rtb;

    my @rqloc = $qLoc->[0];
    my @rtloc = $tLoc->[0];
    for my $i (0..@$qLoc-1) {
        my ($rqb, $rqe) = @{$qLoc->[$i]};
        my ($rtb, $rte) = @{$tLoc->[$i]};
        my ($rql, $rtl) = ([$rqb, $rqe], [$rtb, $rte]);
        print "blk err: $id $qId($rqb-$rqe) <> $tId($rtb-$rte)\n" unless $rqe-$rqb==$rte-$rtb;
        
        next if $i == 0;
        my ($prqb, $prqe) = @{$rqloc[-1]};
        my ($prtb, $prte) = @{$rtloc[-1]};
        if($rqb <= $prqe || $rtb <= $prte) {
            my $plen = $prqe - $prqb + 1;
            my $len = $rqe - $rqb + 1;
            if($plen < $len) {
                $rqloc[-1] = $rql;
                $rtloc[-1] = $rtl;
            }
            print "$id: $qId\[$prqb-$prqe, $rqb-$rqe] $tId\[$prtb-$prte, $rtb-$rte]\n", 
        } else {
            push @rqloc, $rql;
            push @rtloc, $rtl;
        }
    }
}
close $fhi;


__END__
