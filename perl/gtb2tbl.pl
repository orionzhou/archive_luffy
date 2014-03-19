#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gtb2tbl.pl - convert a Gtb file to TBL (Bioconductor-Gviz) format

=head1 SYNOPSIS
  
  gtb2tbl.pl [-help] [-in input-file] [-out output-file]

  Options:
    -help   brief help message
    -in     input file
    -out    output file

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
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;

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

print $fho join("\t", qw/chromosome start end width strand feature gene exon transcript symbol/)."\n";
while( <$fhi> ) {
  chomp;
  next if /(^id)|(^\#)|(^\s*$)/;
  my $ps = [ split "\t" ];
  next unless @$ps >= 18;
  my ($id, $par, $chr, $beg, $end, $srd, $locES, $locIS, $locCS, $loc5S, $loc3S, $phaseS, $src, $conf, $cat1, $cat2, $cat3, $note) = @$ps;
  if($cat2 eq "mRNA") {
    $locCS || die "no CDS for $id\n";
    $locES || die "no exon for $id\n";
    my $reloc = locStr2Ary($locES);
    my @rlocs = map {[@$_, 'protein_coding']} @{locStr2Ary($locCS)};
    push @rlocs, map {[@$_, 'utr5']} @{locStr2Ary($loc5S)} if $loc5S;
    push @rlocs, map {[@$_, 'utr3']} @{locStr2Ary($loc3S)} if $loc3S;
    for my $i (1..@rlocs) {
      my ($rb, $re, $type) = @{$rlocs[$i-1]};
      if($loc5S && $type eq "utr3" && @rlocs==3) {
        $rb += 1;
      }
      my $idx = first_index {$_->[0] <= $rb && $re <= $_->[1]} @$reloc;
      $idx > -1 || die "$rb-$re($type) not in exon\n";
      my ($b, $e) = $srd eq "+" ? ($beg+$rb-1, $beg+$re-1) : ($end-$re+1, $end-$rb+1);
      print $fho join("\t", $chr, $b, $e, $e-$b+1, $srd, $type, $par, $id.($idx+1), $id, $note)."\n";
    }
  } else {
    $locES || die "no exon for $id\n";
    my $reloc = locStr2Ary($locES);
    for my $i (0..@$reloc-1) {
      my ($rb, $re) = @{$reloc->[$i]};
      my ($b, $e) = $srd eq "+" ? ($beg+$rb-1, $beg+$re-1) : ($end-$re+1, $end-$rb+1);
      print $fho join("\t", $chr, $b, $e, $e-$b+1, $srd, $cat2, $par, $id.$i, $id, $note)."\n";
    }
  }
}
close $fho;



__END__
