#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  galtiling.pl - coordinate tiling of a Gal file

=head1 SYNOPSIS
  
  galtiling.pl [-help] [-in input-file] [-out output-file]

  Options:
      -help   brief help message
      -in     input file - needs to be sorted by qId
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

my ($fi, $fo) = ('') x 2;
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

print $fho join("\t", @HEAD_GAL)."\n";

my @rows;
my $tag = '';
while(<$fhi>) {
    chomp;
    next if /(^id)|(^\#)|(^\s*$)/;
    my @ps = split "\t";
    $ps[-2] = locStr2Ary($ps[-2]);
    $ps[-1] = locStr2Ary($ps[-1]);
    if($ps[1] ne $tag && $tag ne "") {
        gal_tiling(\@rows, $fho);
        @rows = ();
    }
    push @rows, \@ps;
    $tag = $ps[1];
}
gal_tiling(\@rows, $fho) if @rows > 0;
close $fhi;
close $fho;

sub gal_tiling {
    my ($rows, $fho) = @_;
    my $loc1 = [ map {[$_->[2], $_->[3]]} @$rows ];
    my $scores = [ map {$_->[-3]} @$rows ];
    my $ref = tiling($loc1, $scores, 2);
    my $qId = $rows->[0]->[1];
    
    for (@$ref) {
        my ($qBeg, $qEnd, $i) = @$_;
        my $row = $rows->[$i];
        my ($id, $qSize, $tId, $tSize, $score) = @$row[0,5,6,10,16];
        my ($qb, $qe, $qSrd) = @$row[2..4];
        my ($tb, $te, $tSrd) = @$row[7..9];
        my $srd = $qSrd eq $tSrd ? "+" : "-";
       
        my ($rqloc, $rtloc) = @$row[17..18];
        my ($rqb, $rqe) = ($qBeg-$qb+1, $qEnd-$qb+1);
        my $nrqloc = trimLoc($rqloc, $rqb, $rqe);
        ($rqb, $rqe) = ($nrqloc->[0]->[0], $nrqloc->[-1]->[1]);
        ($qBeg, $qEnd) = ($qb+$rqb-1, $qb+$rqe-1);
    
        my ($rtb, $rte) = map {coordTransform($_, $rqloc, "+", $rtloc, "+")} ($rqb, $rqe);
        my $nrtloc = trimLoc($rtloc, $rtb, $rte);
        my $tBeg = $srd eq "-" ? $te-$rte+1 : $tb+$rtb-1;
        my $tEnd = $srd eq "-" ? $te-$rtb+1 : $tb+$rte-1;
        
        my $nrqlocS = locAry2Str( [map {[$_->[0]-$rqb+1, $_->[1]-$rqb+1]} @$nrqloc] );
        my $nrtlocS = locAry2Str( [map {[$_->[0]-$rtb+1, $_->[1]-$rtb+1]} @$nrtloc] );
        print $fho join("\t", $id, $qId, $qBeg, $qEnd, "+", $qSize,
            $tId, $tBeg, $tEnd, $srd, $tSize, 
            ('') x 5, $score, $nrqlocS, $nrtlocS)."\n";
    }
}


__END__
