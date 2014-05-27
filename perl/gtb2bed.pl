#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gtb2bed.pl - convert a Gtb file to BED file

=head1 SYNOPSIS
  
  gtb2bed.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (Gtb format)
    -o (--out)    output prefix

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use File::Path qw/make_path remove_tree/;
use File::Basename;
use Common;
use Location;

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

my $fhi;
if ($fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "Can't open file $fi: $!\n";
}

my $dir = $fo;
-d $dir || make_path($dir);

my $hh;
for my $pre (qw/mrna cds intron utr5 utr3/) {
  open (my $fho, ">$dir/$pre.bed") || die "cannot write $dir/$pre.bed\n";
  $hh->{$pre} = $fho;
}

while(<$fhi>) {
  chomp;
  /^(id)|(\#)/ && next;
  my ($id, $par, $chr, $beg, $end, $srd, 
    $locES, $locIS, $locCS, $loc5S, $loc3S, $phase, 
    $src, $conf, $cat1, $cat2, $cat3, $note) = split "\t";
  next if $cat2 ne "mRNA";
  die "no CDS-loc for $id\n" unless $locCS;
  my $len = $end - $beg + 1;
  my $hl = {
    'mrna' => "1-$len",
    'cds' => $locCS,
    'intron' => $locIS,
    'utr5' => $loc5S,
    'utr3' => $loc3S,
  };
  for my $suf (keys(%$hl)) {
    $hl->{$suf} ne '' || next;
    my $rloc = locStr2Ary($hl->{$suf});
    my @loc;
    if($srd eq "-") {
      @loc = map {[$end - $_->[1] + 1, $end - $_->[0] + 1]} @$rloc;
    } else {
      @loc = map {[$beg + $_->[0] - 1, $beg + $_->[1] - 1]} @$rloc;
    }
    my $fho = $hh->{$suf};
    for (@loc) {
      my ($beg, $end) = @$_;
      print $fho join("\t", $chr, $beg - 1, $end, $id, "", $srd)."\n";
    }
  }
}
close $fhi;
for my $pre (keys(%$hh)) {
  my $fho = $hh->{$pre};
  close $fho;
  runCmd("bedSort $fo/$pre.bed $fo/$pre.bed");
}

__END__
