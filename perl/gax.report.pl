#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gax.report.pl - qry loc => tgt loc

=head1 SYNOPSIS
  
  gax.report.pl [-help] [-qry qry-genome] [-tgt tgt-genome] [-loc loc-str]

  Options:
    -h (--help)   brief help message
    -l (--loc)    location string
    -q (--qry)    qry genome (default: HM034)
    -t (--tgt)    tgt genome (default: HM101)
    -o (--out)    output file (default: stdout)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Bio::DB::Fasta;
use Tabix;
use Common;
use Gal;
use Location;
use Data::Dumper;
use List::Util qw/min max sum/;

my ($qry, $tgt) = ('HM034', 'HM101');
my ($fo) = ('');
my $locS;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "qry|q=s"  => \$qry,
  "tgt|t=s"  => \$tgt,
  "out|o=s"  => \$fo,
  "loc|l=s"  => \$locS,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$locS;

my $data = $ENV{'data'};
my $qry_fas = "$data/genome/$qry/11_genome.fas";
my $tgt_fas = "$data/genome/$tgt/11_genome.fas";
my $qdb = Bio::DB::Fasta->new($qry_fas);
my $tdb = Bio::DB::Fasta->new($tgt_fas);

my $dir = "$data/misc3/$qry\_$tgt/23_blat";
chdir $dir || die "cannot chdir to $dir\n";

my ($chr, $beg, $end);
if($locS =~ /^([\w\-\_]+)\:([0-9e\.,]+)\-([0-9e\.,]+)$/i) {
  ($chr, $beg, $end) = ($1, $2, $3);
  $beg =~ s/,//g;
  $end =~ s/,//g;
} elsif($locS =~ /^([\w\-]+)\:(\d+)$/) {
  ($chr, $beg) = ($1, $2);
} elsif($locS =~ /^([\w\-]+)$/) {
  $chr = $1;
  $beg = 1;
} else {
  die "unknown string: $locS\n";
}

my $rev;
if(defined $qdb->length($chr)) {
  $rev = 1;
} elsif(defined $tdb->length($chr)) {
  $rev = 0;
} else {
  die "no db size for $chr\n";
}

($qdb, $tdb) = ($tdb, $qdb) if $rev;
$end ||= $tdb->length($chr);

my $fi = $rev ? "41.9/gax.gz" : "31.9/gax.gz";
#$fi = $rev ? "41.5/gax.gz" : "31.5/gax.gz";

my $fho;
if ($fo eq "" || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

my $gax = Tabix->new(-data => $fi);
my $ary = read_gax($gax, $chr, $beg, $end);
$ary = [ sort {$a->[0] <=> $b->[0]} @$ary ];

my @ids = map {$_->[0]} @$ary;
my $ref = group(\@ids);

my $tsize = $tdb->length($chr);
print $fho "==========\n";
print $fho "$chr:$beg-$end|$tsize\n";
print $fho "==========\n";

my $h;
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
  my $qsize = $qdb->length($qid);
  my $str = "$tb-$te|$tsrd $qid:$qb-$qe|$qsrd|$qsize $tlen|$qlen|$ali";
  $h->{$id} = [$ali, $str];
}
for my $id (sort {$h->{$b}->[0] <=> $h->{$a}->[0]} keys(%$h)) {
  print $fho $h->{$id}->[1]."\n";
}


