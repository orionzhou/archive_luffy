#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gtbcatte.pl - fix TE categorization for a Gtb file

=head1 SYNOPSIS
  
  gtbcatte.pl [-help] [-in input-Gtb] [-out output-Gtb]

  Options:
      -help   brief help message
      -in     input Gtb
      -out    output Gtb

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;

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

for my $i (0..$t->lastRow) {
  my ($id, $par, $chr, $beg, $end, $srd, $locE, $locI, $locC, $loc5, $loc3, $phase, $src, $conf, $cat1, $cat2, $cat3, $note) = $t->row($i);
  if($cat1 =~ /^gene$/i) {
      $t->setElm($i, "cat3", "gene");
  } elsif($cat1 =~ /transposable_element/i) {
      $t->setElm($i, "cat3", "TE");
  }
}
print $fho $t->tsv(1);
close $fho;



__END__
