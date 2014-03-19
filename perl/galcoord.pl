#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  galcoord.pl - transform the coordinates of a Gal file

=head1 SYNOPSIS
  
  galcoord.pl [-help] [-in input-file] [-opt option] [-out output-file]

  Options:
    -help   brief help message
    -in     input file
    -out    output file
    -opt    option (qry / tgt)
    -qry    query-size  file [required if opt='qry']
    -tgt    target-size file [required if opt='tgt']

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;
use Seq;

my ($fi, $fo, $opt) = ('') x 3;
my ($fq, $ft) = ('') x 2;
my ($fhi, $fho);
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
  "opt|p=s"  => \$opt,
  "qry|q=s"  => \$fq,
  "tgt|t=s"  => \$ft
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo || !$opt;

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

my ($hq, $ht);
if( -s $fq ) {
  my $t = readTable(-in=>$fq, -header=>1);
  $hq = { map {$t->elm($_, "id") => $t->elm($_, "length")} (0..$t->nofRow-1) };
}
if( -s $ft ) {
  my $t = readTable(-in=>$ft, -header=>1);
  $ht = { map {$t->elm($_, "id") => $t->elm($_, "length")} (0..$t->nofRow-1) };
}

print $fho join("\t", @HEAD_GAL)."\n";

while(<$fhi>) {
  chomp;
  next if /(^id)|(^\#)|(^\s*$)/;
  my @ps = split "\t";
  if($opt =~ /^qry$/i) {
    my ($qId, $qBeg, $qEnd) = @ps[1..3];
    if($qId =~ /^(\w+)\-([0-9e\+]+)\-([0-9e\+]+)$/) {
      my ($qi, $beg, $end) = ($1, $2, $3);
      $ps[1] = $qi;
      $ps[2] = $beg + $qBeg - 1;
      $ps[3] = $beg + $qEnd - 1;
      die "no size for $qi\n" unless exists $hq->{$qi};
      $ps[5] = $hq->{$qi};
    }
  } elsif($opt =~ /^tgt$/i) {
    my ($tId, $tBeg, $tEnd) = @ps[6..8];
    if($tId =~ /^(\w+)\-([0-9e\+]+)\-([0-9e\+]+)$/) {
      my ($ti, $beg, $end) = ($1, $2, $3);
      $ps[6] = $ti;
      $ps[7] = $beg + $tBeg - 1;
      $ps[8] = $beg + $tEnd - 1;
      die "no size for $ti\n" unless exists $ht->{$ti};
      $ps[10] = $ht->{$ti};
    }
  } else {
    die "unknown opt: $opt\n";
  }
  print $fho join("\t", @ps)."\n";
}
close $fhi;
close $fho;


__END__
