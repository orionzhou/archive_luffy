#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gax2bed.pl - convert a Gax (long) file to BED file

=head1 SYNOPSIS
  
  gax2bed.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file
    -o (--out)    output file
    -p (--opt)    option ('qry' or 'tgt', default: 'tgt')

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Location;

my ($fi, $fo) = ('') x 2;
my $opt = "tgt";
my ($fhi, $fho);
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
  "opt|p=s"  => \$opt,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$opt;

if ($fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "Can't open file $fi: $!\n";
}

if ($fo eq "" || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

$opt = lc($opt);
while( <$fhi> ) {
  chomp;
  next if /(^id)|(^\#)|(^\s*$)/;
  my ($tid, $tb, $te, $tsrd, $id, $qid, $qb, $qe, $qsrd) = split "\t";
  if($opt eq "qry") {
    print $fho join("\t", $qid, $qb - 1, $qe)."\n";
  } elsif($opt eq "tgt") {
    print $fho join("\t", $tid, $tb - 1, $te)."\n";
  } else {
    die "unknow optiotn: $opt\n";
  }
}
close $fhi;
close $fho;


__END__
