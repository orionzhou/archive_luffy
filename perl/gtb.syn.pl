#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gtb.syn.pl - extract synonymous positions from a Gtb file

=head1 SYNOPSIS
  
  gtb.syn.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (Gtb)
    -o (--out)    output file (Tbl)

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

my ($fhi, $fho);
if ($fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "cannot read file $fi\n";
}

if ($fo eq '' || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

while(<$fhi>) {
  chomp;
  /(^id)|(^\#)|(^\s*$)/i && next;
  my $ps = [ split("\t", $_, -1) ];
  @$ps >= 18 || die "not 19 fileds:\n$_\n";
  my ($id, $par, $chr, $beg, $end, $srd, 
    $locES, $locIS, $locCS, $loc5S, $loc3S, $phaseS,
    $src, $conf, $cat1, $cat2, $cat3, $note) = @$ps;
  $cat1 eq "mRNA" || next;
  die "no CDS-loc for $id\n" unless $locCS;
  my $rloc = locStr2Ary($locCS);
  my @phases = split(",", $phaseS);
  @$rloc == @phases || die "not enough phases: $locCS, $phaseS\n";
  my $spos = $rloc->[-1]->[1];
  for my $i (0..$#phases) {
    my $phase = $phases[$i];
    my ($rb, $re) = @{$rloc->[$i]};
    for my $j ($rb..$re) {
      my $rpos = (3-$phase)%3 + $j-$rb+1;
      $rpos % 3 == 0 || next;
      $j == $spos && next;
      my $pos = $srd eq "-" ? $end - $j + 1 : $beg + $j - 1;
      print $fho join("\t", $chr, $pos)."\n";
    }
  }
}
close $fhi;
close $fho;

__END__
