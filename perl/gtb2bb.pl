#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gtb2bb.pl - convert a Gtb file to BIGBED format

=head1 SYNOPSIS
  
  gtb2bb.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (Gtb format)
    -o (--out)    output file (BigBed format)
    -s (--size)   chromosome size file

=head1 DESCRIPTION

  This program converts an input Gtb file to an output BIGBED file

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use Common;
use Location;

my ($fi, $fo, $fs) = ('') x 3;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
  "size|s=s" => \$fs
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fs;

my $base = basename($fi, ".gtb");
$fo ||= "$base.bb";


my ($fhi, $fho);
if ($fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "Can't open file $fi: $!\n";
}
open ($fho, ">$base.bed") || die "cannot write $base.bed\n";

my $hcol = {
  "default" => "0",
  "gene" => "0,0,120",
  "TE"   => "64,64,64",
  "NBS"  => "255,128,0",
  "CRP"  => "0,128,255"
};

#print $fho "#track name=gene_models itemRgb=On useScore=0\n";
while(<$fhi>) {
  chomp;
  /^(id)|(\#)/ && next;
  my ($id, $par, $chr, $beg, $end, $srd, 
    $locES, $locIS, $locCS, $loc5S, $loc3S, $phase, 
    $src, $conf, $cat1, $cat2, $cat3, $note) = split "\t";
  die "no CDS-loc for $id\n" unless $locCS;
  
  my $idStr = $id;
  $idStr .= "|$note" if $note;
  my $h = { map {$_=>1} split(/\|/, $cat3) };
  my $col = exists $h->{"NBS"} ? $hcol->{"NBS"} :
      exists $h->{"CRP"} ? $hcol->{"CRP"} :
      exists $h->{"TE"} ? $hcol->{"TE"}:
      exists $h->{"gene"} ? $hcol->{"gene"} : $hcol->{"default"};

  my @locs = sort {$a->[0] <=> $b->[0]} @{locStr2Ary($locES)};
  my $n = @locs;
  
  my $rloc = [ sort {$a->[0] <=> $b->[0]} @{locStr2Ary($locCS)} ];
  my ($rtBeg, $rtEnd) = ($rloc->[0]->[0], $rloc->[-1]->[1]);
  
  my @begs;
  my @lens;
  my ($tBeg, $tEnd);
  if($srd eq "+") {
    @begs = map {$_->[0] - 1} @locs;
    @lens = map {$_->[1] - $_->[0] + 1} @locs;
    ($tBeg, $tEnd) = ($beg+$rtBeg-1, $beg+$rtEnd-1);
  } else {
    @begs = reverse map {$end-$beg+1 - $_->[1]} @locs;
    @lens = reverse map {$_->[1] - $_->[0] + 1} @locs;
    ($tBeg, $tEnd) = ($end-$rtEnd+1, $end-$rtBeg+1);
  }
  print $fho join("\t", $chr, $beg - 1, $end, $idStr, 0, $srd, 
    $tBeg - 1, $tEnd, $col,
    $n, join(",", @lens), join(",", @begs) )."\n";
}
close $fhi;
close $fho;


runCmd("bedSort $base.bed $base.bed", 1);
runCmd("bedToBigBed -tab $base.bed $fs $fo", 1);
runCmd("rm $base.bed", 1);



__END__
