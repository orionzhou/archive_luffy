#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gtb2tbl.pl - convert a Gtb file to Tbl (location) file

=head1 SYNOPSIS
  
  gtb2tbl.pl [-help] [-in input-file] [-out output-file]

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
  open ($fhi, $fi) || die "Can't open file $fi: $!\n";
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
    $locES, $locIS, $locCS, $loc5S, $loc3S, $phase, 
    $src, $conf, $cat1, $cat2, $cat3, $note) = @$ps;
  $cat1 eq "mRNA" || next;
  die "no CDS-loc for $id\n" unless $locCS;
  my $len = $end - $beg + 1;
  my $hl = {
    'mrna' => "1-$len",
    'cds'  => $locCS,
    'utr5' => $loc5S,
    'utr3' => $loc3S,
    'intron' => $locIS,
  };
  for my $type (keys(%$hl)) {
    my $rloc = locStr2Ary($hl->{$type});
    my $loc = $srd eq "-" ?
      [ map {[$end - $_->[1] + 1, $end - $_->[0] + 1]} @$rloc ] :
      [ map {[$beg + $_->[0] - 1, $beg + $_->[1] - 1]} @$rloc ];
    for (@$loc) {
      print $fho join("\t", $chr, @$_, $srd, $id, $type, $cat2)."\n";
    }
  }
}
close $fhi;
close $fho;
runCmd("sort -k1,1 -k2,2n -k3,3n $fo -o $fo");
runCmd("bgzip -c $fo > $fo.gz");
runCmd("tabix -s 1 -b 2 -e 3 -f $fo.gz");

__END__
