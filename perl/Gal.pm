package Gal;
use strict;
use Tabix;
use Data::Dumper;
use Common;
use Location;
use Seq;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter/;
@EXPORT = qw/@HEAD_GAL @HEAD_GALL
    read_gax read_snp_cnt read_idm_cnt/;
@EXPORT_OK = qw//;

our @HEAD_GAL = qw/id tId tBeg tEnd tSrd tSize
  qId qBeg qEnd qSrd qSize
  ali mat mis qN tN ident score tLoc qLoc/;

  my $ps = [];
  my ($id, $tId, $tBeg, $tEnd, $tSrd, $tSize, 
    $qId, $qBeg, $qEnd, $qSrd, $qSize,
    $ali, $mat, $mis, $qN, $tN, $ident, $score, $tLocS, $qLocS) = @$ps;

sub read_gax {
  my ($con, $id, $beg, $end, $opt) = @_;
  my $iter = $con->query($id, $beg-1, $end);
  my @ary;
  while (my $line = $con->read($iter)) {
    my @ps = split("\t", $line);
    my ($tid, $tb, $te, $tsrd, $id, $qid, $qb, $qe, $qsrd);
    if($opt eq 't') {
      ($tid, $tb, $te, $tsrd, $id, $qid, $qb, $qe, $qsrd) = @ps;
    } else {
      ($qid, $qb, $qe, $qsrd, $id, $tid, $tb, $te, $tsrd) = @ps;
    }
    
    if($opt eq "t" && $tb < $beg) {
      $qb += $beg - $tb if $tsrd eq $qsrd;
      $qe -= $beg - $tb if $tsrd ne $qsrd;
      $tb = $beg;
    } 
    if($opt eq "q" && $qb < $beg) {
      $tb += $beg - $qb if $tsrd eq $qsrd;
      $te -= $beg - $qb if $tsrd ne $qsrd;
      $qb = $beg;
    }
    if($opt eq "t" && $te > $end) {
      $qe -= $te - $end if $tsrd eq $qsrd;
      $qb += $te - $end if $tsrd ne $qsrd;
      $te = $end;
    }
    if($opt eq "q" && $qe > $end) {
      $te -= $qe - $end if $tsrd eq $qsrd;
      $tb += $qe - $end if $tsrd ne $qsrd;
      $qe = $end;
    }
#    my $alen = $te - $tb + 1;
#    print join("\t", $qid, $qb, $qe, $tid, $tb, $te, $alen)."\n";
    push @ary, [$id, $tid, $tb, $te, $tsrd, $qid, $qb, $qe, $qsrd];
  }
  return \@ary;
}
sub read_snp_cnt {
  my ($con, $id, $beg, $end) = @_;
  my $iter = $con->query($id, $beg, $end);
  my $h;
  while (my $line = $con->read($iter)) {
    my ($tid, $tpos, $id, $qid, $qpos, $qsrd, $tbase, $qbase) = 
      split("\t", $line);
    $h->{$id} ||= 0;
    $h->{$id} ++;
  }
  return $h;
}
sub read_idm_cnt {
  my ($con, $id, $beg, $end) = @_;
  my $iter = $con->query($id, $beg, $end);
  my $h;
  while (my $line = $con->read($iter)) {
    my ($tid, $tb, $te, $id, $qid, $qb, $qe) = split("\t", $line);
    my ($tgap, $qgap) = ($te - $tb + 1, $qe - $qb + 1);
    my $tgapo = $tgap > 0 ? 1 : 0;
    my $qgapo = $qgap > 0 ? 1 : 0;
    
    $h->{$id} ||= [0, 0, 0, 0];
    $h->{$id}->[0] += $tgapo;
    $h->{$id}->[1] += $tgap;
    $h->{$id}->[2] += $qgapo;
    $h->{$id}->[3] += $qgap;
  }
  return $h;
}

1;
__END__
