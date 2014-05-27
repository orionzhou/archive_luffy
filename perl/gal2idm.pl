#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gal2idm.pl - Call Ins/Del/Mnp from a Gal file

=head1 SYNOPSIS
  
  gal2idm.pl [-help] [-in input-Gal] [-idm InDel-output] [-sv SV-output]

  Options:
    -h (--help)   brief help message
    -i (--in)     input Gal
    -d (--idm)    InDel output file 
    -s (--sv )    SV output file 
    -q (--qry)    query-seq file 
    -t (--tgt)    target-seq file
    -l (--len)    length threshold for Indel/SV (default: 50)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Location;
use Common;
use Seq;
use File::Path qw/make_path remove_tree/;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;

my ($fi, $fd, $fs) = ('') x 3;
my ($fq, $ft) = ('') x 2; 
my ($len) = (50); 
my ($fhi, $fhd, $fhs);
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "idm|d=s"  => \$fd,
  "sv|s=s"   => \$fs,
  "qry|q=s"  => \$fq,
  "tgt|t=s"  => \$ft,
  "len|l=i"  => \$len,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fd || !$fs;
pod2usage(2) if !$fq || !$ft;

if ($fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "Can't read $fi\n";
}

open ($fhd, ">$fd") || die "Can't write file $fd\n";
open ($fhs, ">$fs") || die "Can't write file $fs\n";

while( <$fhi> ) {
  chomp;
  next if /(^id)|(^\#)|(^\s*$)/;
  my $ps = [ split "\t" ];
  next unless @$ps == 20;
  my ($id, $tId, $tBeg, $tEnd, $tSrd, $tSize, 
    $qId, $qBeg, $qEnd, $qSrd, $qSize,
    $ali, $mat, $mis, $qN, $tN, $ident, $score, $tLocS, $qLocS) = @$ps;

  my ($rqLoc, $rtLoc) = (locStr2Ary($qLocS), locStr2Ary($tLocS));
  $rqLoc = [ sort {$a->[0] <=> $b->[0]} @$rqLoc ];
  $rtLoc = [ sort {$a->[0] <=> $b->[0]} @$rtLoc ];
  @$rqLoc == @$rtLoc || die "unequal pieces\n";
  my $nBlock = @$rqLoc;
  my @lens = map {$_->[1] - $_->[0] + 1} @$rqLoc;
 
  $nBlock > 1 || next;
  for my $i (1..$nBlock-1) {
    my ($rqb, $rqe) = ($rqLoc->[$i-1]->[1], $rqLoc->[$i]->[0]);
    my ($rtb, $rte) = ($rtLoc->[$i-1]->[1], $rtLoc->[$i]->[0]);
    my $qins = $rqe - $rqb - 1;
    my $tins = $rte - $rtb - 1;
    my ($qb, $qe) = $qSrd eq "-" ? ($qEnd - $rqe + 1, $qEnd - $rqb + 1) : 
      ($qBeg + $rqb - 1, $qBeg + $rqe - 1);
    my ($tb, $te) = $tSrd eq "-" ? ($tEnd - $rte + 1, $tEnd - $rtb + 1) : 
      ($tBeg + $rtb - 1, $tBeg + $rte - 1);

    if($qins >= $len || $tins >= $len) {
      print $fhs join("\t", $tId, $tb, $te, $id, $qId, $qb, $qe)."\n";
    } else {
      my $tBase = $tins > 0 ? 
        seqRet([[$tb+1, $te-1]], $tId, $tSrd, $ft) : '';
      my $qBase = $qins > 0 ? 
        seqRet([[$qb+1, $qe-1]], $qId, $qSrd, $fq) : '';
      print $fhd join("\t", $tId, $tb, $te, $tBase, $qBase, 
        $id, $qId, $qb, $qe)."\n";
    }
  }
}
close $fhi;
close $fhd;
close $fhs;


__END__
