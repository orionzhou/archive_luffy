package Psl;
use strict;
use Data::Dumper;
use Common;
use Location;
use Seq;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter/;
@EXPORT = qw/psl2Gal chain2Gal net2Gal/;
@EXPORT_OK = qw//;

sub psl2Gal {
    my ($fhi, $fho, $fq, $ft) = @_;
    print $fho join("\t", qw/id qId qBeg qEnd qSrd qLen 
        tId tBeg tEnd tSrd tLen
        match misMatch baseN ident e score qLoc tLoc/)."\n";
    my $id = 0;
    while(<$fhi>) {
        chomp;
        my @ps = split " ";
        next unless @ps == 21;
        next if /^(psLayout)|(match)/;
        $id ++;
        my ($match, $misMatch, $repMatch, $baseN, 
            $qNumIns, $qBaseIns, $tNumIns, $tBaseIns, $qSrd, 
            $qId, $qSize, $qBeg, $qEnd, $tId, $tSize, $tBeg, $tEnd, 
            $blockNum, $blockSizes, $qBegs, $tBegs) = @ps;
        my $tSrd = "+";
        my @qBegs = split(",", $qBegs);
        my @tBegs = split(",", $tBegs);
        my @blockSizes = split(",", $blockSizes);
        
        $blockNum == @qBegs || die "unequal pieces\n";
        $blockNum == @tBegs || die "unequal pieces\n";
        $blockNum == @blockSizes || die "unequal pieces\n";
        my $alnLen = $match + $misMatch + $repMatch + $baseN;
        $alnLen == sum(@blockSizes) || die "block size error:$alnLen/".sum(@blockSizes)."\n";
        $alnLen + $qBaseIns == $qEnd-$qBeg || die "qLen error\n";
        $alnLen + $tBaseIns == $tEnd-$tBeg || die "hLen error\n";
        
#        my $seqT = seqRet([[$tBeg+1, $tEnd]], $tId, $tSrd, $ft);
#        my $seqQ = seqRet([[$qBeg+1, $qEnd]], $qId, $qSrd, $fq);
        my $score = $match;
        my (@qLoc, @tLoc);

        for my $i (0..$blockNum-1) {
            my $len = $blockSizes[$i];
            my $tb = $tBegs[$i] + 1;
            my $te = $tb + $len - 1;
            
            my ($qb, $qe);
            if($qSrd eq "+") {
                $qb = $qBegs[$i] + 1;
                $qe = $qb + $len - 1;
            } else {
                $qSrd eq "-" || die "unknown strand $qSrd\n";
                $qe = $qSize - $qBegs[$i];
                $qb = $qe - $len + 1;
            }
            
            my $rtb = $tBegs[$i] - $tBeg;
            my $rqb = $qSrd eq "-" ? $qBegs[$i]-($qSize-$qEnd) : $qBegs[$i]-$qBeg;
            push @tLoc, [$rtb+1, $rtb+$len];
            push @qLoc, [$rqb+1, $rqb+$len];
        }

#        my $tSeq = getSubSeq($seqT, \@tLoc);
#        my $qSeq = getSubSeq($seqQ, \@qLoc);
#        ($match, $misMatch, $baseN) = seqCompare($qSeq, $tSeq);
        my ($tLocS, $qLocS) = (locAry2Str(\@tLoc), locAry2Str(\@qLoc));
        print $fho join("\t", $id, $qId, $qBeg+1, $qEnd, $qSrd, $qEnd-$qBeg,
            $tId, $tBeg+1, $tEnd, $tSrd, $tEnd-$tBeg,
            ('')x5, $score, $qLocS, $tLocS)."\n";
    }
    close $fhi;
    close $fho;
}

1;
__END__
