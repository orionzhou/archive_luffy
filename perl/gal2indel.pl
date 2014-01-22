#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gal2indel.pl - Call InDels from a Gal file and a Neti file

=head1 SYNOPSIS
  
  gal2indel.pl [-help] [-in input-file] [-neti neti-file] [-out output-file]

  Options:
      -help   brief help message
      -in     input-Gal
      -out    output-Tbl
      -neti   Neti-file

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
    "neti|q=s"  => \$fn,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;
pod2usage(2) if !$fn;

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

my $hf;
my $hr;
open (my $fhn, "<$fn") || die "Can't read $fn: $!\n";
while(<$fhn>) {
    chomp;
    next if /(^id)|(^\#)|(^\s*$)/;
    my $ps = [ split "\t" ];
    my ($id, $par, $lev, $type) = @$ps;
    if($lev == 1) {
        $hf->{$id} = [];
    } elsif(exists($hf->{$par})) {
        $hr->{$id} = $par;
    }
}

my @rows;
my $hi;
while( <$fhi> ) {
    chomp;
    next if /(^id)|(^\#)|(^\s*$)/;
    my $ps = [ split "\t" ];
    next unless @$ps == 19;
    my ($id, $qId, $qBeg, $qEnd, $qSrd, $qSize, $tId, $tBeg, $tEnd, $tSrd, $tSize,
        $match, $misMatch, $baseN, $ident, $e, $score, $qLocS, $tLocS) = @$ps;

    $hi->{$id} ||= 0;
    $hi->{$id} ++;
    $id .= ".".$hi->{$id} if $hi->{$id} > 1;

    if(exists($hf->{$id})) {
        $ps->[0] = $id;
        push @rows, \@$ps;
    } elsif(exists($hr->{$id})) {
        my $par = $hr->{$id};
        push @{$hf->{$par}}, [$tBeg, $tEnd, $tSrd] if exists $hf->{$par};
    }
}
close $fhi;

for (@rows) {
    my ($id, $qId, $qBeg, $qEnd, $qSrd, $qSize, $tId, $tBeg, $tEnd, $tSrd, $tSize,
        $match, $misMatch, $baseN, $ident, $e, $score, $qLocS, $tLocS) = @$_;
    
    my ($rqLoc, $rtLoc) = (locStr2Ary($qLocS), locStr2Ary($tLocS));
    @$rqLoc == @$rtLoc || die "unequal pieces\n";
    my $nBlock = @$rqLoc;
    my @lens = map {$_->[1] - $_->[0] + 1} @$rqLoc;

    my $tLoc = $tSrd eq "-" ? [ map {[$tEnd-$_->[1]+1, $tEnd-$_->[0]+1]} @$rtLoc ]
        : [ map {[$tBeg+$_->[0]-1, $tBeg+$_->[1]-1]} @$rtLoc ]; 
    my $qLoc = $qSrd eq "-" ? [ map {[$qEnd-$_->[1]+1, $qEnd-$_->[0]+1]} @$rqLoc ]
        : [ map {[$qBeg+$_->[0]-1, $qBeg+$_->[1]-1]} @$rqLoc ];
    
    my ($tgLoc) = posDiff([[$tBeg, $tEnd]], $tLoc);
    
    defined $hf->{$id} || die "$id\n";
    my @cLoc = @{$hf->{$id}};
    if(@cLoc == 0) {
        for (@$tgLoc) {
            my ($tb, $te) = @$_;
            print $fho join("\t", $tId, $tb, $te, $te-$tb+1)."\n";
        }
    } else {
        my ($dLoc) = posDiff($tgLoc, \@cLoc);
        for (@$dLoc) {
            my ($tb, $te) = @$_;
            print $fho join("\t", $tId, $tb, $te, $te-$tb+1)."\n";
        }
    }
}
close $fho;


__END__
