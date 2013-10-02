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
@EXPORT = qw/gal2Psl chain2Gal net2Gal
    gal_tiling gal_check_fix
    gal_breakdown gal_filter
    gal_complete gal_indel/;
@EXPORT_OK = qw//;

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

sub gal2Psl {
    my ($fhi, $fho, $hq, $ht) = @_;
    while(<$fhi>) {
        chomp;
        next if /(^id)|(^\#)|(^\s*$)/;
        my ($id, $qId, $qBeg, $qEnd, $qSrd, $qLen, $tId, $tBeg, $tEnd, $tSrd, $tLen,
            $match, $misMatch, $baseN, $ident, $e, $score, $qLocS, $tLocS) = split "\t";
        die join("\t", $id, $qId, $qBeg, $qEnd, $qSrd, $qLen, $tId, $tBeg, $tEnd, $tSrd, $tLen,
            $match, $misMatch, $baseN, $ident, $e, $score)."\n";
        die "size unknown for $qId\n" unless exists $hq->{$qId};
        die "size unknown for $tId\n" unless exists $ht->{$tId};
        my ($qSize, $tSize) = ($hq->{$qId}, $ht->{$tId});
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
}

sub gal_tiling {
#id qId qBeg qEnd qSrd qLen tId tBeg tEnd tSrd tLen match misMatch baseN ident e score qLoc tLoc
    my ($rows, $fho) = @_;
    my $loc1 = [ map {[$_->[2], $_->[3]]} @$rows ];
    my $scores = [ map {$_->[-3]} @$rows ];
    my $ref = tiling($loc1, $scores, 2);
    my $qId = $rows->[0]->[1];
    
    for (@$ref) {
        my ($qBeg, $qEnd, $i) = @$_;
        my $row = $rows->[$i];
        my ($id, $tId, $score) = @$row[0,6,16];
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
        
        my $qLen = $qEnd - $qBeg + 1;
        my $nrqlocS = locAry2Str( [map {[$_->[0]-$rqb+1, $_->[1]-$rqb+1]} @$nrqloc] );
        my $tLen = $tEnd - $tBeg + 1;
        my $nrtlocS = locAry2Str( [map {[$_->[0]-$rtb+1, $_->[1]-$rtb+1]} @$nrtloc] );
        print $fho join("\t", $id, $qId, $qBeg, $qEnd, "+", $qLen,
            $tId, $tBeg, $tEnd, $srd, $tLen, 
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

sub gal_check_fix {
    my ($ps) = @_;
    my ($id, $qId, $qBeg, $qEnd, $qSrd, $qLen, $tId, $tBeg, $tEnd, $tSrd, $tLen,
        $match, $misMatch, $baseN, $ident, $e, $score, $qLocS, $tLocS) = @$ps;
    my ($qLoc, $tLoc) = (locStr2Ary($qLocS), locStr2Ary($tLocS));
    my $srd = ($qSrd eq $tSrd) ? "+" : "-";

    my @rqloc = $qLoc->[0];
    my @rtloc = $tLoc->[0];
    for my $i (0..@$qLoc-1) {
        next if $i == 0;
        my ($prqb, $prqe) = @{$rqloc[-1]};
        my ($prtb, $prte) = @{$rtloc[-1]};

        my ($rqb, $rqe) = @{$qLoc->[$i]};
        my ($rtb, $rte) = @{$tLoc->[$i]};
        my ($rql, $rtl) = ([$rqb, $rqe], [$rtb, $rte]);
        if($rqb <= $prqe || $rtb <= $prte) {
            my $plen = $prqe - $prqb + 1;
            my $len = $rqe - $rqb + 1;
            if($plen < $len) {
                $rqloc[-1] = $rql;
                $rtloc[-1] = $rtl;
            }
        } else {
            push @rqloc, $rql;
            push @rtloc, $rtl;
        }
    }
    $ps->[17] = locAry2Str(\@rqloc);
    $ps->[18] = locAry2Str(\@rtloc);
    return $ps;
}

sub gal_complete {
    my ($ps, $fq, $ft) = @_;
    my ($id, $qId, $qBeg, $qEnd, $qSrd, $qLen, $tId, $tBeg, $tEnd, $tSrd, $tLen,
        $match, $misMatch, $baseN, $ident, $e, $score, $qLocS, $tLocS) = @$ps;
    my ($qLoc, $tLoc) = (locStr2Ary($qLocS), locStr2Ary($tLocS));
    @$qLoc == @$tLoc || die "unequal pieces\n";
    my $nBlock = @$qLoc;

    my $seqT = seqRet([[$tBeg, $tEnd]], $tId, $tSrd, $ft);
    my $seqQ = seqRet([[$qBeg, $qEnd]], $qId, $qSrd, $fq);
    my $tSeq = getSubSeq($seqT, $tLoc);
    my $qSeq = getSubSeq($seqQ, $qLoc);
    ($match, $misMatch, $baseN) = seqCompare($qSeq, $tSeq);
    @$ps[11..13] = ($match, $misMatch, $baseN);
    return $ps;
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



1;
__END__
