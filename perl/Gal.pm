package Gal;
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
@EXPORT = qw/gal2psl
    gal_tiling 
    gal_check gal_fix
    gal_break gal_filter
    gal_complete gal_indel/;
@EXPORT_OK = qw//;

sub gal2psl {
    my ($ps) = @_;
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
    return [ $match, $misMatch, $repMatch, $baseN, 
        $qNumIns, $qBaseIns, $tNumIns, $tBaseIns, $srd, 
        $qId, $qSize, $qBeg-1, $qEnd, $tId, $tSize, $tBeg-1, $tEnd, 
            $nBlock, $blockSizes, $qBegs, $tBegs ];
}

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
sub gal_tiling_o {
    my ($fi, $fo) = @_;
    my $t = readTable(-in=>$fi, -header=>1);
    
    open(FHO, ">$fo") or die "cannot write $fo\n";
    print FHO join("\t", $t->header)."\n";

    my $ref = group( $t->colRef("idc") );
    for my $idc (sort {$a<=>$b} (keys(%$ref))) {
        my ($idx, $cnt) = @{$ref->{$idc}};
        my ($qId, $qSrd, $hId, $hSrd) = map {$t->elm($idx, $_)} qw/qId qSrd hId hSrd/;
        my $qLoc1 = [ map { [$t->elm($_, "qBeg"), $t->elm($_, "qEnd")] } ($idx..$idx+$cnt-1) ];
        my $hLoc1 = [ map { [$t->elm($_, "hBeg"), $t->elm($_, "hEnd")] } ($idx..$idx+$cnt-1) ];
        my $qLen1 = [ map {$t->elm($_, "qLen")} ($idx..$idx+$cnt-1) ];
        my $stat = tiling($qLoc1, $qLen1, 2);
        $stat = [ reverse @$stat ] if $qSrd eq "-";
        my $idb = 0;
        for (@$stat) {
            my ($qBeg, $qEnd, $i) = @$_;
            my ($qb, $qe) = map {$t->elm($idx+$i, $_)} qw/qBeg qEnd/;
            my ($hb, $he) = map {$t->elm($idx+$i, $_)} qw/hBeg hEnd/;
            my ($hBeg, $hEnd) = $qSrd eq "+" ? ($hb+($qBeg-$qb), $hb+($qEnd-$qb)) : 
                ($he-($qEnd-$qb), $he-($qBeg-$qb));
            my $qLen = $qEnd - $qBeg + 1;
            my $hLen = $hEnd - $hBeg + 1;
            $qLen == $hLen || die "len err: $idc-$i\n";
            print FHO join("\t", $idc, ++$idb, 'aln', $qId, $qBeg, $qEnd, $qSrd, $qLen,
                $hId, $hBeg, $hEnd, $hSrd, $hLen, ('') x 8)."\n";
        }
    }
    close FHO;
}

sub gal_check {
    my ($ps) = @_;
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
sub gal_fix {
    my ($ps) = @_;
    my ($id, $qId, $qBeg, $qEnd, $qSrd, $qSize, $tId, $tBeg, $tEnd, $tSrd, $tSize,
        $match, $misMatch, $baseN, $ident, $e, $score, $qLocS, $tLocS) = @$ps;
    my ($qLoc, $tLoc) = (locStr2Ary($qLocS), locStr2Ary($tLocS));
    my $srd = ($qSrd eq $tSrd) ? "+" : "-";
    
    my @qlens = map {$_->[1]-$_->[0]+1} @$qLoc;
    my $ref = tiling($qLoc, \@qlens, 2);
    my (@rqloc, @rtloc);
    for (@$ref) {
        my ($rqb, $rqe, $idx) = @$_;

        my ($qb, $qe) = @{$qLoc->[$idx]};
        my ($tb, $te) = @{$tLoc->[$idx]};
        my $rtb = $rqb - $qb + $tb;
        my $rte = $rqe - $qb + $tb;

        if(@rqloc == 0 || $rtb > $rtloc[-1]->[1]) { 
            push @rqloc, [$rqb, $rqe];
            push @rtloc, [$rtb, $rte];
        }
    }
    $ps->[17] = locAry2Str(\@rqloc);
    $ps->[18] = locAry2Str(\@rtloc);
    return $ps;
}

sub gal_complete {
    my ($ps, $fq, $ft) = @_;
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
    ($match, $misMatch, $baseN) = seqCompare($qSeq, $tSeq);
    @$ps[11..13] = ($match, $misMatch, $baseN);
    return $ps;
}

sub gal_break {
    my ($ps, $gap, $fho) = @_;
    my ($id, $qId, $qBeg, $qEnd, $qSrd, $qSize, $tId, $tBeg, $tEnd, $tSrd, $tSize,
        $match, $misMatch, $baseN, $ident, $e, $score, $qLocS, $tLocS) = @$ps;
    my ($qLoc, $tLoc) = (locStr2Ary($qLocS), locStr2Ary($tLocS));
    @$qLoc == @$tLoc || die "unequal pieces\n";
    my $nBlock = @$qLoc;
    
    my @idxs = [0, @$qLoc-1];
    for my $i (0..@$qLoc-1) {
        next if $i == 0;
        my ($rqb, $rqe) = @{$qLoc->[$i]};
        my ($rtb, $rte) = @{$tLoc->[$i]};
        my ($prqb, $prqe) = @{$qLoc->[$i-1]};
        my ($prtb, $prte) = @{$tLoc->[$i-1]};
        ($prqe < $rqb && $prte < $rtb) ||
            die "error: $id $qId\[$prqb-$prqe, $rqb-$rqe] $tId\[$prtb-$prte, $rtb-$rte]\n";
        if($rqb - $prqe - 1 >= $gap) {
            $idxs[-1]->[1] = $i-1;
            push @idxs, [$i, @$qLoc-1];
        }
    }
    my $cnt = 1;
    
    for (@idxs) {
        my ($idxb, $idxe) = @$_;
        my @rql = @$qLoc[$idxb..$idxe];
        my @rtl = @$tLoc[$idxb..$idxe];
        my ($rqb, $rqe) = ($rql[0]->[0], $rql[-1]->[1]);
        my ($rtb, $rte) = ($rtl[0]->[0], $rtl[-1]->[1]);
        my $qb = $qSrd eq "-" ? $qEnd-$rqe+1 : $qBeg+$rqb-1;
        my $qe = $qSrd eq "-" ? $qEnd-$rqb+1 : $qBeg+$rqe-1;
        my $tb = $tSrd eq "-" ? $tEnd-$rte+1 : $tBeg+$rtb-1;
        my $te = $tSrd eq "-" ? $tEnd-$rtb+1 : $tBeg+$rte-1;
        my @nrql = map {[$_->[0]-$rqb+1, $_->[1]-$rqb+1]} @rql;
        my @nrtl = map {[$_->[0]-$rtb+1, $_->[1]-$rtb+1]} @rtl;
        my $qlen = $rqe - $rqb + 1;
        my $tlen = $rte - $rtb + 1;
        my ($qls, $tls) = (locAry2Str(\@nrql), locAry2Str(\@nrtl));
        print $fho join("\t", $id.".".($cnt++), $qId, $qb, $qe, $qSrd, $qSize, 
            $tId, $tb, $te, $tSrd, $tSize, ('')x6, $qls, $tls);
    }
}


sub gal_filter {
    my ($fi, $foc, $fob, $co_aln, $co_qcov, $co_hcov) = @_;
    $co_aln ||= 100;
    $co_qcov ||= 0.01;
    $co_hcov ||= 0.01;
    my $t = readTable(-in=>$fi, -header=>1);
   
    open(FHOB, ">$fob") or die "cannot write $fob\n";
    print FHOB join("\t", $t->header)."\n";
    
    open(FHOC, ">$foc") or die "cannot write $foc\n";
    print FHOC join("\t", qw/idc qId qBeg qEnd qSrd qLen hId hBeg hEnd hSrd hLen aLen qCov hCov/)."\n";

    my $i = 0;
    my $ref = group( $t->colRef("idc") );
    for my $idc (sort(keys(%$ref))) {
        my ($idx, $cnt) = @{$ref->{$idc}};
        my $ts = $t->subTable([$idx..$idx+$cnt-1], [$t->header]);
        my ($qId, $qSrd, $hId, $hSrd) = map {$ts->elm(0, $_)} qw/qId qSrd hId hSrd/;
        my ($qBeg, $qEnd) = (min($ts->col("qBeg")), max($ts->col("qEnd")));
        my ($hBeg, $hEnd) = (min($ts->col("hBeg")), max($ts->col("hEnd")));
        my $qLen = $qEnd - $qBeg + 1;
        my $hLen = $hEnd - $hBeg + 1;
        my $aLen = sum($ts->col("qLen"));
        my $qCov = $aLen / $qLen;
        my $hCov = $aLen / $hLen;
        next if $aLen < $co_aln || $qCov < $co_qcov || $hCov < $co_hcov;
        $i ++;
        print FHOB $ts->tsv(0);
        print FHOC join("\t", $idc, $qId, $qBeg, $qEnd, $qSrd, $qLen, $hId, $hBeg, $hEnd, 
            $hSrd, $hLen, $aLen, $qCov, $hCov)."\n";
    }
    close FHOB;
    close FHOC;
    printf "%s out of %s aln chains passed\n", $i, scalar(keys(%$ref));
}


sub gal_indel {
    my ($fi, $fo) = @_;
    my $t = readTable(-in=>$fi, -header=>1);
   
    open(FHO, ">$fo") or die "cannot write $fo\n";
    print FHO join("\t", $t->header)."\n";

    my $ref = group( $t->colRef("idc") );
    for my $idc (sort(keys(%$ref))) {
        my ($idx, $cnt) = @{$ref->{$idc}};
        next if $cnt < 2;
        for my $i (2..$cnt) {
            my ($qb1, $qe1, $hb1, $he1) = map {$t->elm($idx+$i-2, $_)} qw/qBeg qEnd hBeg hEnd/;
            my ($qb2, $qe2, $hb2, $he2) = map {$t->elm($idx+$i-1, $_)} qw/qBeg qEnd hBeg hEnd/;
            my ($idc2, $idb, $type, $qId, $qBeg, $qEnd, $qSrd, $qLen,
                $hId, $hBeg, $hEnd, $hSrd, $hLen) = $t->row($idx+$i-2);
            my ($qb, $qe, $hb, $he) = ($qe1, $qb2, $he1, $hb2);
            ($qb, $qe) = ($qe2, $qb1) if $qSrd eq "-";
            $qb < $qe || die "qry error: $idc-$idb\n";
            $hb < $he || die "hit error: $idc-$idb\n";
            $qLen = $qe - $qb - 1;
            $hLen = $he - $hb - 1;
            $type = $qLen == 0 ? "del" : $hLen == 0 ? "ins" : "rpl";
            $idb = ($idb + $t->elm($i-1, "idb")) / 2;
            print FHO join("\t", $idc, $idb, $type, $qId, $qb, $qe, $qSrd, $qLen,
                $hId, $hb, $he, $hSrd, $hLen, ('') x 8)."\n";
        }
    }
    close FHO;
}



1;
__END__
