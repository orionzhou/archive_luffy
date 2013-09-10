package Gal;
use strict;
use Data::Dumper;
use Common;
use Seq;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter/;
@EXPORT = qw/psl2Gal chain2Gal
    gal_tiling gal_rm_inv gal_breakdown gal_filter
    gal_validate gal_indel/;
@EXPORT_OK = qw//;

sub psl2Gal {
    my ($fi, $fo) = @_;
#  my @cols = qw/matches misMatches repMatches nCount qNumInsert qBaseInsert tNumInsert tBaseInsert strand qName qSize qStart qEnd tName tSize tStart tEnd blockCount blockSizes qStarts tStarts/;
    open(FHI, "<$fi") or die "cannot read $fi\n";
    open(FHO, ">$fo") or die "cannot write $fo\n";
    print FHO join("\t", qw/idc idb type qId qBeg qEnd qSrd qLen 
        hId hBeg hEnd hSrd hLen
        match misMatch baseN ident e score qSeq hSeq/)."\n";
    my $idc = 0;
    while(<FHI>) {
        chomp;
        my @ps = split " ";
        next unless @ps == 21;
        next if /^(psLayout)|(match)/;
        $idc ++;
        my ($match, $misMatch, $repMatch, $baseN, 
            $qNumIns, $qBaseIns, $hNumIns, $hBaseIns, $qSrd, 
            $qId, $qLen, $qBeg, $qEnd, $hId, $hLen, $hBeg, $hEnd, 
            $blockNum, $blockSizes, $qBegs, $hBegs) = @ps;
        $qBeg += 1;
        $hBeg += 1;
        my $hSrd = "+";
        my @qBegs = split(",", $qBegs);
        my @hBegs = split(",", $hBegs);
        my @blockSizes = split(",", $blockSizes);
        
        $blockNum == @qBegs || die "unequal pieces\n";
        $blockNum == @hBegs || die "unequal pieces\n";
        $blockNum == @blockSizes || die "unequal pieces\n";
        my $alnLen = $match + $misMatch + $repMatch + $baseN;

        my $idb = 0;
        for my $i (0..$blockNum-1) {
            $idb ++;
            my ($hb, $blockSize) = ($hBegs[$i]+1, $blockSizes[$i]);
            my $he = $hb + $blockSize - 1;
            
            my ($qb, $qe);
            if($qSrd eq "+") {
                $qb = $qBegs[$i] + 1;
                $qe = $qb + $blockSize - 1;
            } else {
                $qSrd eq "-" || die "unknown strand $qSrd\n";
                $qe = $qLen - $qBegs[$i];
                $qb = $qe - $blockSize + 1;
            }
            
            print FHO join("\t", $idc, $idb, "aln", $qId, $qb, $qe, $qSrd, $blockSize,
                $hId, $hb, $he, $hSrd, $blockSize,
                ('') x 8)."\n";
        }
        $alnLen == sum(@blockSizes) || die "block size error:$alnLen/".sum(@blockSizes)."\n";
        $alnLen + $qBaseIns == $qEnd-$qBeg+1 || die "qLen error\n";
        $alnLen + $hBaseIns == $hEnd-$hBeg+1 || die "hLen error\n";
    }
    close FHI;
    close FHO;
}

sub chain2Gal {
    my ($fi, $fob) = @_;
    open(FHI, "<$fi") or die "cannot read $fi\n";
    
    open(FHOB, ">$fob") or die "cannot write $fob\n";
    print FHOB join("\t", qw/idc idb level qId qBeg qEnd qSrd qLen 
        hId hBeg hEnd hSrd hLen
        match misMatch baseN dq dh score qSeq hSeq/)."\n";
    while( <FHI> ) {
        chomp;
        my @ps = split /\s/;
        if($ps[0] eq "chain") {
            my ($score, $hId, $hSize, $hSrd, $hBeg, $hEnd, 
                $qId, $qSize, $qSrd, $qBeg, $qEnd, $idc) = @ps[1..$#ps];
            $hSrd eq "+" || die "$idc: hSrd -\n";

            my ($qd, $hd) = (0, 0); 
            my $idb = 1;
            while( <FHI> ) {
                last if /^\s*\n$/;
                my @pps = split /\s/;
                my $len = $pps[0];
                my ($dh, $dq) = (0, 0);
                ($dh, $dq) = @pps[1..2] if @pps >= 3;

                my ($hb, $he) = ($hBeg+$hd+1, $hBeg+$hd+$len);
                my ($qb, $qe) = ($qBeg+$qd+1, $qBeg+$qd+$len);
                ($qb, $qe) = ($qSize-$qe+1, $qSize-$qb+1) if $qSrd eq "-";
                my $sco = $idb == 1 ? $score : "";
                print FHOB join("\t", $idc, $idb++, '', $qId, $qb, $qe, $qSrd, $len,
                    $hId, $hb, $he, $hSrd, $len, ('') x 3, $dq, $dh, $sco, '', '')."\n";
                $qd += $len + $dq;
                $hd += $len + $dh;
            }
        }
    }
    close FHI;
    close FHOB;
}

sub gal_tiling {
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
sub gal_rm_inv {
    my ($fi, $fo) = @_;
    my $t = readTable(-in=>$fi, -header=>1);

    my @idxs;
    my ($idc_p, $qId_p, $qBeg_p, $qEnd_p, $qSrd_p, $qLen_p, $hId_p, $hBeg_p, $hEnd_p, $hSrd_p)
        = map {$t->elm(0, $_)} qw/idc qId qBeg qEnd qSrd qLen hId hBeg hEnd hSrd/;
    for my $i (1..$t->nofRow-1) {
        my ($idc, $idb, $type, $qId, $qBeg, $qEnd, $qSrd, $qLen,
            $hId, $hBeg, $hEnd, $hSrd, $hLen,
            $match, $misMatch, $baseN, $ident, $e, $score, $qSeq, $hSeq) = $t->row($i);
        if($idc == $idc_p) {
            ($qId eq $qId_p && $qSrd eq $qSrd_p) || die "qId error: $idc-$idb\n";
            ($hId eq $hId_p && $hSrd eq $hSrd_p) || die "hId error: $idc-$idb\n";
            my ($flag, $idx) = (0, $i);
            if($qSrd eq "+") {
                if($qEnd_p >= $qBeg || $hEnd_p >= $hBeg) {
                    $flag = 1;
                }
            } else {
                $qSrd eq "-" || die "unknown strand $qSrd\n";
                if($qEnd >= $qBeg_p || $hEnd_p >= $hBeg) {
                    $flag = 1;
                }
            }
            if($flag == 1) {
                $idx = $i-1 if ($qLen_p < $qLen);
                push @idxs, $idx;
            }
            if($flag == 0 || ($flag && $idx == $i-1)) {
                ($idc_p, $qId_p, $qBeg_p, $qEnd_p, $qSrd_p, $qLen_p, $hId_p, $hBeg_p, $hEnd_p, $hSrd_p) = 
                ($idc, $qId, $qBeg, $qEnd, $qSrd, $qLen, $hId, $hBeg, $hEnd, $hSrd); 
            }
        } else {
            ($idc_p, $qId_p, $qBeg_p, $qEnd_p, $qSrd_p, $qLen_p, $hId_p, $hBeg_p, $hEnd_p, $hSrd_p) = 
            ($idc, $qId, $qBeg, $qEnd, $qSrd, $qLen, $hId, $hBeg, $hEnd, $hSrd); 
        }
    }
    $t->delRows(\@idxs);
    printf "removed %d lines\n", $#idxs+1;
    open(FHO, ">$fo") or die "cannot write to $fo\n";
    print FHO $t->tsv(1);
    close FHO;
}
sub gal_breakdown {
    my ($fi, $fo, $co_gap) = @_;
    $co_gap ||= 1_000_000;
    
    my $t = readTable(-in=>$fi, -header=>1);
    open(FHO, ">$fo") or die "cannot write to $fo\n";
    print FHO join("\t", $t->header)."\n";

    my $id = 1;
    my ($idc_p, $qId_p, $qBeg_p, $qEnd_p, $qSrd_p, $qLen_p, $hId_p, $hBeg_p, $hEnd_p, $hSrd_p)
        = map {$t->elm(0, $_)} qw/idc qId qBeg qEnd qSrd qLen hId hBeg hEnd hSrd/;
    for my $i (1..$t->nofRow-1) {
        my ($idc, $idb, $type, $qId, $qBeg, $qEnd, $qSrd, $qLen,
            $hId, $hBeg, $hEnd, $hSrd, $hLen,
            $match, $misMatch, $baseN, $ident, $e, $score, $qSeq, $hSeq) = $t->row($i);
        if($idc == $idc_p) {
            ($qId eq $qId_p && $qSrd eq $qSrd_p) || die "qId error: $idc-$idb\n";
            ($hId eq $hId_p && $hSrd eq $hSrd_p) || die "hId error: $idc-$idb\n";
            my $hGap = $hBeg - $hEnd_p - 1;
            my $qGap = $qBeg - $qEnd_p - 1;
            $qGap = $qBeg_p - $qEnd - 1 if $qSrd eq "-";
            $id ++ if $qGap > $co_gap || $hGap > $co_gap;
        } else {
            $id ++;
        }
        print FHO join("\t", $id, $idb, $type, $qId, $qBeg, $qEnd, $qSrd, $qLen,
            $hId, $hBeg, $hEnd, $hSrd, $hLen,
            ("") x 8 )."\n";
        ($idc_p, $qId_p, $qBeg_p, $qEnd_p, $qSrd_p, $qLen_p, $hId_p, $hBeg_p, $hEnd_p, $hSrd_p) = 
        ($idc, $qId, $qBeg, $qEnd, $qSrd, $qLen, $hId, $hBeg, $hEnd, $hSrd); 
    }
    close FHO;
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

sub gal_validate {
    my ($fi, $fo, $fq, $ft) = @_;
    my $t = readTable(-in=>$fi, -header=>1);
   
    open(FHO, ">$fo") or die "cannot write $fo\n";
    print FHO join("\t", $t->header)."\n";
    for my $i (0..$t->nofRow-1) {
        my ($idc, $idb, $type, $qId, $qBeg, $qEnd, $qSrd, $qLen,
            $hId, $hBeg, $hEnd, $hSrd, $hLen,
            $match, $misMatch, $baseN, $ident, $e, $score, $qSeq, $hSeq) = $t->row($i);
        ($match, $misMatch, $baseN) = (0, 0, 0);
        my $seqQ = seqRet([[$qBeg, $qEnd]], $qId, $qSrd, $fq);
        length($seqQ) == $qLen || die "$idc: $qLen\n";
        my $seqT = seqRet([[$hBeg, $hEnd]], $hId, $hSrd, $ft);
        length($seqT) == $hLen || die "$idc: $hLen\n";
        for my $j (0..$qLen-1) {
            my $chQ = uc(substr($seqQ, $j, 1));
            my $chT = uc(substr($seqT, $j, 1));
            if($chQ eq "N" || $chT eq "N") {
                $baseN ++;
            } elsif($chQ eq $chT) {
                $match ++;
            } else {
                $misMatch ++;
            }
        }
        print FHO join("\t", $idc, $idb, $type, $qId, $qBeg, $qEnd, $qSrd, $qLen,
            $hId, $hBeg, $hEnd, $hSrd, $hLen,
            $match, $misMatch, $baseN, $ident, $e, $score, $qSeq, $hSeq)."\n";
    }
    close FHO;
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

sub mtb_expand {
    my ($fi, $fo, $fq, $ft) = @_;
    my $t = readTable(-in=>$fi, -header=>1);
    open(FHO, ">$fo") or die "cannot write to $fo\n";
    print FHO join("\t", qw/qId qBeg qEnd qSrd/)."\n";
    for my $i (0..$t->nofRow-1) {
        my ($id, $qId, $qBeg, $qEnd, $qSrd, $qLen, $qNumIns, $qBaseIns,
            $hId, $hBeg, $hEnd, $hSrd, $hLen, $hNumIns, $hBaseIns, 
            $match, $misMatch, $repMatch, $baseN, 
            $qLocS, $hLocS, $qIndelS, $hIndelS) = $t->row($i);
        my ($qLoc, $hLoc) = map {locStr2Ary($_)} ($qLocS, $hLocS);
        my ($match2, $misMatch2, $baseN2) = (0, 0, 0);
        my $qSeq = seqRet($qLoc, $qId, $qSrd, $fq);
        my $hSeq = seqRet($hLoc, $hId, $hSrd, $ft);
        for my $j (0..length($qSeq)-1) {
            my $qCh = uc(substr($qSeq, $j, 1));
            my $hCh = uc(substr($hSeq, $j, 1));
            if($qCh !~ /[ATCG]/ || $hCh !~ /[ATCG]/) {
                $baseN2 ++;
            } elsif($qCh eq $hCh) {
                $match2 ++;
            } else {
                $misMatch2 ++;
            }
        }
        if( $match != $match2 || $misMatch != $misMatch2 || $baseN != $baseN2 ) {
            print join("\t", $id, $qId, $hId, $match, $match2, $repMatch, $misMatch, $misMatch2, $baseN, $baseN2)."\n";
        }
    } 
}

1;
__END__
