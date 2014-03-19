#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  refinesv.pl - refine SV calls

=head1 SYNOPSIS
  
  refinesv.pl [-help] [-in input-file] [-out output-file]

  Options:
      -h (--help)   brief help message
      -i (--in)     input file
      -o (--out)    output file
      -1 (--1)      comp1 prefix
      -2 (--2)      comp2 prefix

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin";
use Tabix;
use Common;
use Data::Dumper;
use File::Path qw/make_path remove_tree/;
use File::Basename;
use List::Util qw/min max sum/;

my $help_flag;
my ($fi, $fo) = ("", "");
my ($pre1, $pre2) = ("", "");

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "in|i=s"  => \$fi,
  "out|o=s" => \$fo,
  "1=s"     => \$pre1,
  "2=s"     => \$pre2
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo || !$pre1 || !$pre2;

my $tgax = Tabix->new(-data=>"$pre1.gax.gz");
my $tsnp = Tabix->new(-data=>"$pre1.vnt/snp.gz");
my $tidm = Tabix->new(-data=>"$pre1.vnt/idm.gz");
my $qgax = Tabix->new(-data=>"$pre2.gax.gz");
my $qsnp = Tabix->new(-data=>"$pre2.vnt/snp.gz");
my $qidm = Tabix->new(-data=>"$pre2.vnt/idm.gz");

open(my $fhi, "<$fi") or die "cannot read $fi\n";
open(my $fho, ">$fo") or die "cannot write $fo\n";
my $cnt = 0;
while(<$fhi>) {
  chomp;
  next if /^(\#)|(tId)/;
  $cnt ++;
  next unless $cnt == 25306;
  my ($tid, $tb, $te, $id, $qid, $qb, $qe, $qlen, $tlen) = split "\t";
  print join("\t", $tid, $tb, $te, $id, $qid, $qb, $qe, $qlen, $tlen)."\n";
  my ($tlb, $tle) = ($tb - 199, $tb);
  my ($trb, $tre) = ($te, $te + 199);
  my ($tsb, $tse) = ($tb + 1, $te - 1);
  my ($qlb, $qle) = ($qb - 199, $qb);
  my ($qrb, $qre) = ($qe, $qe + 199);
  my ($qsb, $qse) = ($qb + 1, $qe - 1);

  my ($ins, $del) = (0, 0);
  if($tlen > 0) {
    my $h = read_comp($tid, $tsb, $tse, $tgax, $tsnp, $tidm, 't');
    for my $id (keys(%$h)) {
      my ($tid, $tb, $te, $tsrd, $qid, $qb, $qe, $qsrd,
        $alen, $mm, $gapo, $gap, $score) = @{$h->{$id}};
      print join(" ", "$tid:$tb-$te", "$qid:$qb-$qe", 
        "m:$alen", "mm:$mm", "gapo:$gapo", "gap:$gap", $score)."\n";
    }
  }
  if($qlen > 0) {
    my $h = read_comp($qid, $qsb, $qse, $qgax, $qsnp, $qidm, 'q');
    for my $id (keys(%$h)) {
      my ($tid, $tb, $te, $tsrd, $qid, $qb, $qe, $qsrd,
        $alen, $mm, $gapo, $gap, $score) = @{$h->{$id}};
      print join("  ", "$tid:$tb-$te", "$qid:$qb-$qe", $score)."\n";
    }
  }
  print $fho join("\t", $tid, $tb, $te, $id, $qid, $qb, $qe, 
    $qlen, $tlen)."\n";
}

sub read_comp {
  my ($id, $beg, $end, $gax, $snp, $idm, $opt) = @_;
  my $hg = read_gax($gax, $id, $beg, $end, $opt);
  my $hs = read_snp_cnt($snp, $id, $beg, $end);
  my $hi = read_idm_cnt($idm, $id, $beg, $end);

  my ($sco_m, $sco_mm, $sco_gapo, $sco_gape) = (1, -1, -2, -1);
  for my $id (keys(%$hg)) {
    my ($tid, $tb, $te, $tsrd, $qid, $qb, $qe, $qsrd, $alen)
      = @{$hg->{$id}};
    my ($mm, $tgapo, $tgap, $qgapo, $qgap) = (0) x 5;
    $mm = $hs->{$id} if exists $hs->{$id};
    ($tgapo, $tgap, $qgapo, $qgap) = @{$hi->{$id}} if exists $hi->{$id};
    my ($gapo, $gap) = ($tgapo + $qgapo, $tgap + $qgap);
    my $score = ($alen-$mm) * $sco_m + $mm * $sco_mm +
      $gapo * $sco_gapo + ($gap - $gapo) * $sco_gape;
    $hg->{$id} = [$tid, $tb, $te, $tsrd, $qid, $qb, $qe, $qsrd,
      $alen, $mm, $gapo, $gap, $score];
  }
  return $hg;
}
sub read_gax {
  my ($con, $id, $beg, $end, $opt) = @_;
  my $iter = $con->query($id, $beg-1, $end);
  my $h;
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
      $te = $beg;
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
#    print join("\t", $qid, $qb, $qe, $tid, $tb, $te)."\n";
    my $alen = $te - $tb + 1;
    $h->{$id} ||= [$tid, $tb, $te, $tsrd, $qid, $qb, $qe, $qsrd, 0];
    $h->{$id}->[1] = $tb if $tb < $h->{$id}->[1];
    $h->{$id}->[2] = $te if $te > $h->{$id}->[2];
    $h->{$id}->[5] = $qb if $qb < $h->{$id}->[5];
    $h->{$id}->[6] = $qe if $qe > $h->{$id}->[6];
    $h->{$id}->[8] += $alen;
  }
  return $h;
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

