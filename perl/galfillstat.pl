#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  galfillstat.pl - Fills the stat fields (match, misMatch, baseN) in a GAL file

=head1 SYNOPSIS
  
  galfillstat.pl [-help] [-in input-file] [-qry qry-fasta] [-tgt tgt-fasta] [-out output-file]

  Options:
      -h (--help)   brief help message
      -i (--in)     input Gal
      -o (--out)    output Gal
      -q (--qry)    query fasta
      -t (--tgt)    target fasta

=cut
  
#### END of POD documentation.
#-----------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Time::HiRes qw/gettimeofday tv_interval/;
use Location;
use Seq;
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
pod2usage(2) if !$fq || !$ft;

if ($fi eq "" || $fi eq "stdin" || $fi eq "-") {
    $fhi = \*STDIN;
} else {
    open ($fhi, $fi) || die "Can't open file $fi: $!\n";
}

if ($fo eq "" || $fo eq "stdout" || $fo eq "-") {
    $fho = \*STDOUT;
} else {
    open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

my $t0 = [gettimeofday];
print $fho join("\t", @HEAD_GAL)."\n";

my $cnt = 1;
while( <$fhi> ) {
    chomp;
    next if /(^id)|(^\#)|(^\s*$)/;
    my $ps = [ split "\t" ];
    next unless @$ps == 19;
    my ($id, $qId, $qBeg, $qEnd, $qSrd, $qSize, $tId, $tBeg, $tEnd, $tSrd, $tSize,
        $match, $misMatch, $baseN, $ident, $e, $score, $qLocS, $tLocS) = @$ps;
    my ($rqLoc, $rtLoc) = (locStr2Ary($qLocS), locStr2Ary($tLocS));
    @$rqLoc == @$rtLoc || die "unequal pieces\n";
    my $nBlock = @$rqLoc;

#    my $seqT = seqRet([[$tBeg, $tEnd]], $tId, $tSrd, $ft);
#    my $seqQ = seqRet([[$qBeg, $qEnd]], $qId, $qSrd, $fq);
#    my $tSeq = getSubSeq($seqT, $rtLoc);
#    my $qSeq = getSubSeq($seqQ, $rqLoc);
    my $tLoc = $tSrd eq "-" ? [ map {[$tEnd-$_->[1]+1, $tEnd-$_->[0]+1]} @$rtLoc ]
        : [ map {[$tBeg+$_->[0]-1, $tBeg+$_->[1]-1]} @$rtLoc ]; 
    my $qLoc = $qSrd eq "-" ? [ map {[$qEnd-$_->[1]+1, $qEnd-$_->[0]+1]} @$rqLoc ]
        : [ map {[$qBeg+$_->[0]-1, $qBeg+$_->[1]-1]} @$rqLoc ]; 
    my $tSeq = seqRet($tLoc, $tId, $tSrd, $ft);
    my $qSeq = seqRet($qLoc, $qId, $qSrd, $fq);
    ($match, $misMatch, $baseN) = seqCompare($tSeq, $qSeq);
    $ident = ($match+$misMatch==0) ? 0 : sprintf "%.03f", $match/($match+$misMatch);
    @$ps[11..14] = ($match, $misMatch, $baseN, $ident);
    print $fho join("\t", @$ps)."\n";

    printf "%5d: %.01f min\n", $cnt++, tv_interval($t0, [gettimeofday]) / 60 if $cnt % 1000 == 0;
}
close $fhi;
close $fho;


__END__
