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
@EXPORT = qw/psl2Gal chain2Gal net2Gal
    gal_tiling gal_rm_inv gal_breakdown gal_filter
    gal_validate gal_indel/;
@EXPORT_OK = qw//;

sub psl2Gal {
    my ($fhi, $fho, $fq, $ft) = @_;
    print $fho join("\t", qw/id qId qBeg qEnd qSrd qLen 
        tId tBeg tEnd tSrd tLen
        match misMatch baseN ident e score/)."\n";
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
        
        my $seqT = seqRet([[$tBeg+1, $tEnd]], $tId, $tSrd, $ft);
        my $seqQ = seqRet([[$qBeg+1, $qEnd]], $qId, $qSrd, $fq);
        my $score = $match;

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
            my $tSeq = substr($seqT, $rtb, $len);
            my $rqb = $qSrd eq "-" ? $qBegs[$i]-($qSize-$qEnd) : $qBegs[$i]-$qBeg;
            my $qSeq = substr($seqQ, $rqb, $len);
            my ($match, $misMatch, $baseN) = seqCompare($qSeq, $tSeq);
            
            print $fho join("\t", $id, $qId, $qb, $qe, $qSrd, $len,
                $tId, $tb, $te, $tSrd, $len,
                $match, $misMatch, $baseN, '', '', $score)."\n";
        }
        $alnLen == sum(@blockSizes) || die "block size error:$alnLen/".sum(@blockSizes)."\n";
        $alnLen + $qBaseIns == $qEnd-$qBeg || die "qLen error: $alnLen + $qBaseIns <> $qEnd-$qBeg\n";
        $alnLen + $tBaseIns == $tEnd-$tBeg || die "hLen error\n";
    }
    close $fhi;
    close $fho;
}
sub chain2Gal {
    my ($fhi, $fho, $fq, $ft) = @_;
    print $fho join("\t", qw/id qId qBeg qEnd qSrd qLen 
        tId tBeg tEnd tSrd tLen
        match misMatch baseN ident e score/)."\n";
    while( <$fhi> ) {
        chomp;
        my @ps = split /\s/;
        if($ps[0] eq "chain") {
            my ($score, $tId, $tSize, $tSrd, $tBeg, $tEnd, 
                $qId, $qSize, $qSrd, $qBeg, $qEnd, $id) = @ps[1..$#ps];
            $tSrd eq "+" || die "$id: tSrd -\n";
            
            my $seqT = seqRet([[$tBeg+1, $tEnd]], $tId, $tSrd, $ft);
            my $seqQ = $qSrd eq "-" ? 
                  seqRet([[$qSize-$qEnd+1, $qSize-$qBeg]], $qId, $qSrd, $fq)
                : seqRet([[$qBeg+1, $qEnd]], $qId, $qSrd, $fq);

            my ($td, $qd) = (0, 0); 
            while( <$fhi> ) {
                last if /^\s*\n$/;
                my @pps = split /\s/;
                my $len = $pps[0];
                my ($dt, $dq) = (0, 0);
                ($dt, $dq) = @pps[1..2] if @pps >= 3;

                my ($tb, $te) = ($tBeg+$td+1, $tBeg+$td+$len);
                my ($qb, $qe) = ($qBeg+$qd+1, $qBeg+$qd+$len);
                ($qb, $qe) = ($qSize-$qe+1, $qSize-$qb+1) if $qSrd eq "-";

                my $tSeq = substr($seqT, $td, $len);
                my $qSeq = substr($seqQ, $qd, $len);
                my ($match, $misMatch, $baseN) = seqCompare($qSeq, $tSeq);
                
                my $sco = "";
                print $fho join("\t", $id, $qId, $qb, $qe, $qSrd, $len,
                    $tId, $tb, $te, $tSrd, $len, $match, $misMatch, $baseN, ('')x3)."\n";
                $td += $len + $dt;
                $qd += $len + $dq;
            }
        }
    }
    close $fhi;
    close $fho;
}
sub net2Gal {
    my ($fhi, $fho) = @_;
    print $fho join("\t", qw/level id type qId qBeg qEnd qSrd qLen 
        hId hBeg hEnd hSrd hLen
        match misMatch baseN dq dh score qSeq hSeq/)."\n";
    
    my ($hId, $hSize, @lines);
    while(<$fhi>) {
        chomp;
        if(/^\#/ || /^\s*$/) {
            next;
        } elsif(/^net (\w+) (\d+)/) {
            @lines == 0 || write_net( parse_net(\@lines, $hId, $hSize), $fho );
            ($hId, $hSize) = ($1, $2);
            @lines = ();
        } else {
            push @lines, $_;
        }
    }
    @lines == 0 || write_net( parse_net(\@lines, $hId, $hSize), $fho );
    close $fhi;
    close $fho;
}
sub write_one_net {
    my ($af, $idx, $fho) = @_;
    my ($lev, $pa, $stat, $statE, $gaps) = @{$af->[$idx]};
    my ($hId, $hb, $hl, $qId, $qSrd, $qb, $ql) = @$stat;
    my @keys = grep {exists($statE->{$_})} qw/id score ali qOver qFar qDup type/;
    my @strs = map {$_." ".$statE->{$_}} @keys;
    my ($id, $score, $type) = map {$statE->{$_}} qw/id score type/;
#    print $fho join('', (' ')x(2*$lev-1)).join(" ", 'fill', $hb, $hl, $qId, $qSrd, $qb, $ql, @strs)."\n";
    for (@$gaps) {
        my ($ghb, $ghl, $gqb, $gql, $chi) = @$_;
        my ($hBeg, $hEnd, $hLen) = ($hb+1, $ghb, $ghb-$hb);
        my ($qBeg, $qEnd, $qLen) = $qSrd eq "-" ?
            ($gqb+$gql+1, $qb+$ql, $qb+$ql-$gqb-$gql) : ($qb+1, $gqb, $gqb-$qb);
        $qLen >= 0 || print "$qb:$ql $gqb:$gql [$qLen] $qSrd $hb:$hl $ghb:$ghl [$hLen]\n"; 
        $hLen >= 0 || print "$qb:$ql $gqb:$gql [$qLen] $qSrd $hb:$hl $ghb:$ghl [$hLen]\n"; 
#        print $fho join('', (' ')x($lev*2)).join(" ", 'gap', $ghb, $ghl, $qId, $qSrd, $gqb, $gql)."\n";
        print $fho join("\t", $lev, $id, $type, $qId, $qBeg, $qEnd, $qSrd, $qLen,
            $hId, $hBeg, $hEnd, "+", $hLen, ('')x3, $gql, $ghl, $score, '', '')."\n";
        for my $idC (@$chi) { 
            write_one_net($af, $idC, $fho);
        }
        ($hb, $hl) = ($ghb+$ghl, $hl-($ghb+$ghl-$hb));
        ($qb, $ql) = $qSrd eq "-" ? 
            ($qb, $gqb-$qb) : ($gqb+$gql, $ql-($gqb+$gql-$qb));
    }
    my ($hBeg, $hEnd, $hLen) = ($hb+1, $hb+$hl, $hl);
    my ($qBeg, $qEnd, $qLen) = ($qb+1, $qb+$ql, $ql); 
    $qLen >= 0 || print "$qb:$ql [$qLen] $qSrd $hb:$hl [$hLen]\n"; 
    $hLen >= 0 || print "$qb:$ql [$qLen] $qSrd $hb:$hl [$hLen]\n"; 
    print $fho join("\t", $lev, $id, $type, $qId, $qBeg, $qEnd, $qSrd, $qLen,
        $hId, $hBeg, $hEnd, "+", $hLen, ('')x3, 0, 0, $score, '', '')."\n";
}
sub write_net {
    my ($af, $fho) = @_;
    my @ids = grep {$af->[$_]->[0] == 1} (0..@$af-1);
    for my $id (@ids) {
        write_one_net($af, $id, $fho);
    }
}
sub parse_net {
    my ($lines, $hId, $hSize) = @_;
    my (@af, @al);
    my $idx = 0;
    for (@$lines) {
        if(/^(\s+)fill (.+)$/) {
            my $level = (length($1)+1) / 2;
            my @ps = split(" ", $2);
            my ($hb, $hl, $qId, $qSrd, $qb, $ql) = @ps[0..5];
            my $stat = [$hId, $hb, $hl, $qId, $qSrd, $qb, $ql];
            my @pps = @ps[6..$#ps];
            my $statE = {map {$pps[$_*2] => $pps[$_*2+1]} (0..@pps/2-1)};
            
            if($level > @al) {
                push @al, $idx;
            } elsif($level == @al) {
                $al[-1] = $idx;
            } else {
                for ($level+1..@al) { pop @al; };
                $al[-1] = $idx;
            }
            my $pa = @al>=2 ? $al[-2] : "-1";
            push @af, [$level, $pa, $stat, $statE, []];
            if($pa >= 0) {
                push @{$af[$pa]->[4]->[-1]->[4]}, $idx;
            }
            $idx ++;
        } elsif(/^(\s+)gap (.+)$/) {
            my $level = length($1) / 2;
            my @ps = split(" ", $2);
            $level <= @al || die "gap level[$level] wrong: $al[-1]\n";
            if($level < @al) {
                for ($level+1..@al) { pop @al; };
            }
            my $id = $al[-1];
            my ($hb, $hl, $qId, $qSrd, $qb, $ql) = @ps;
            my ($qId0, $qSrd0) = @{$af[$id]->[2]}[3..4];
            $qId eq $qId0 || die "fill[$id] $qId0 <> gap[$hb-$hl] $qId\n";
            $qSrd eq $qSrd0 || die "fill[$id] $qSrd0 <> gap[$hb-$hl] $qSrd\n";
            push @{$af[$id]->[4]}, [$hb, $hl, $qb, $ql, []]; 
        } else {
            die "unknonw line: $_\n";
        }
    }
    return \@af;
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
