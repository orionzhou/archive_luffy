#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gtb2bed.pl - convert a Gtb file to BED format

=head1 SYNOPSIS
  
  gtb2bed.pl [-help] [-in input-file] [-out output-file]

  Options:
      -h (--help)   brief help message
      -i (--in)     input Gtb
      -o (--out)    output Bed

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;
use Location;
use Gtb;

my ($fi, $fo) = ('') x 2;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

my ($fhi, $fho);
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

my $t = readTable(-inh=>$fhi, -header=>1);
close $fhi;

my $hcol = {
  "default" => "0",
  "gene" => "0,0,120",
  "TE"   => "64,64,64",
  "NBS"  => "255,128,0",
  "CRP"  => "0,128,255"
};

#print $fho "#track name=gene_models itemRgb=On useScore=0\n";
for my $i (0..$t->nofRow-1) {
  my ($id, $par, $chr, $beg, $end, $srd, $locE, $locI, $locC, $loc5, $loc3, $phase, $src, $conf, $cat1, $cat2, $cat3, $note) = $t->row($i);
  next if $cat2 ne "mRNA";
  die "no CDS-loc for $id\n" unless $locC;
  
  my $idStr = $id;
  $idStr .= "|$note" if $note;
  my $h = { map {$_=>1} split(/\|/, $cat3) };
  my $col = exists $h->{"NBS"} ? $hcol->{"NBS"} :
      exists $h->{"CRP"} ? $hcol->{"CRP"} :
      exists $h->{"TE"} ? $hcol->{"TE"}:
      exists $h->{"gene"} ? $hcol->{"gene"} : $hcol->{"default"};

  my @locs = sort {$a->[0] <=> $b->[0]} @{locStr2Ary($locE)};
  my $n = @locs;
  
  my $rloc = [ sort {$a->[0] <=> $b->[0]} @{locStr2Ary($locC)} ];
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
close $fho;



__END__
