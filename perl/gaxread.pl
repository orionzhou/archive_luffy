#!/usr/bin/perl -w
use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Tabix;
use Common;
use Gal;
use Location;
use Data::Dumper;
use List::Util qw/min max sum/;

my ($fi, $fo) = ('') x 2;
my $locS = "";
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
  "string|s=s"  => \$locS,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi;
pod2usage(2) if !$locS;

my ($chr, $beg, $end);
if($locS =~ /^([\w\-\_]+)\:([0-9e\.,]+)\-([0-9e\.,]+)$/i) {
  ($chr, $beg, $end) = ($1, $2, $3);
  $beg =~ s/,//g;
  $end =~ s/,//g;
} else {
  die "unknown string: $locS\n";
}

my $fho;
if ($fo eq "" || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

my $gax = Tabix->new(-data=>$fi);
my $ary = read_gax($gax, $chr, $beg, $end);
$ary = [ sort {$a->[0] <=> $b->[0]} @$ary ];

my @ids = map {$_->[0]} @$ary;
my $ref = group(\@ids);

for my $id (keys(%$ref)) {
  my ($idxb, $cnt) = @{$ref->{$id}};
  my @rows = map {$ary->[$idxb + $_]} (0..$cnt-1);
  my ($tid, $qid) = ($rows[0]->[1], $rows[0]->[5]);
  my ($tsrd, $qsrd) = ($rows[0]->[4], $rows[0]->[8]);
  my $tb = min( map {$_->[2]} @rows );
  my $te = max( map {$_->[3]} @rows );
  my $qb = min( map {$_->[6]} @rows );
  my $qe = max( map {$_->[7]} @rows );

  my @tl = map {[$_->[2], $_->[3]]} @rows;
  my @ql = map {[$_->[6], $_->[7]]} @rows;
  my @rtl = $tsrd eq "-" ? map {[$te - $_->[1] + 1, $te - $_->[0] + 1]} @tl
    : map {[$tb + $_->[0] - 1, $tb + $_->[1] - 1]} @tl;
  my @rql = $qsrd eq "-" ? map {[$qe - $_->[1] + 1, $qe - $_->[0] + 1]} @ql
    : map {[$qb + $_->[0] - 1, $qb + $_->[1] - 1]} @ql;
  my $ali = locAryLen(\@rtl);
  my $qlen = $qe - $qb + 1;
  my $tlen = $te - $tb + 1;
  print $fho "$tid:$tb-$te|$tsrd|$tlen $qid:$qb-$qe|$qsrd|$qlen $ali\n";
}

