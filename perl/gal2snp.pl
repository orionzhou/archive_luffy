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
use Common;
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

print $fho join("\t", qw/qid qpos tid tpos id qbase tbase/)."\n";
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
    my @lens = map {$_->[1] - $_->[0] + 1} @$rqLoc;

    my $tLoc = $tSrd eq "-" ? [ map {[$tEnd-$_->[1]+1, $tEnd-$_->[0]+1]} @$rtLoc ]
        : [ map {[$tBeg+$_->[0]-1, $tBeg+$_->[1]-1]} @$rtLoc ]; 
    my $qLoc = $qSrd eq "-" ? [ map {[$qEnd-$_->[1]+1, $qEnd-$_->[0]+1]} @$rqLoc ]
        : [ map {[$qBeg+$_->[0]-1, $qBeg+$_->[1]-1]} @$rqLoc ]; 
    my $tSeq = seqRet($tLoc, $tId, $tSrd, $ft);
    my $qSeq = seqRet($qLoc, $qId, $qSrd, $fq);

    my $len = 0;
    my (@tNts, @qNts, @tPoss, @qPoss);
    for my $i (0..$nBlock-1) {
        my ($rtb, $rte) = @{$rtLoc->[$i]};
        my ($rqb, $rqe) = @{$rqLoc->[$i]};
        for my $j (0..$lens[$i]-1) {
            my $tNt = uc(substr($tSeq, $len+$j, 1));
            my $qNt = uc(substr($qSeq, $len+$j, 1));
            my $tPos = $tSrd eq "-" ? $tEnd-($rtb+$j)+1 : $tBeg+($rtb+$j)-1;
            my $qPos = $qSrd eq "-" ? $qEnd-($rqb+$j)+1 : $qBeg+($rtb+$j)-1;
            if($tNt ne "N" && $qNt ne "N" && $tNt ne $qNt) {
                push @tNts, $tNt;
                push @qNts, $qNt;
                push @tPoss, $tPos;
                push @qPoss, $qPos;
            }
        }
        $len += $lens[$i];
    }
    @tPoss == $misMatch || die "$id not $misMatch SNPs: ".scalar(@tPoss)."\n";
    
    for my $i (0..@qPoss-1) {
        my ($qPos, $tPos, $qNt, $tNt) = ($qPoss[$i], $tPoss[$i], $qNts[$i], $tNts[$i]);
        print $fho join("\t", $qId, $qPos, $tId, $tPos, $id, $qNt, $tNt)."\n";
    }
}
close $fhi;
close $fho;


__END__
