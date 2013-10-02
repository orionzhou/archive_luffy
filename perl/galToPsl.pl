#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  galToPsl.pl - convert GAL file to PSL format

=head1 SYNOPSIS
  
  galToPsl.pl [-help] [-in input-file] [-out output-file]

  Options:
      -help   brief help message
      -in     input file
      -out    output file
      -qry    query size file 
      -tgt    target size file

=head1 DESCRIPTION

  This program converts a GAL file to an output PSL file

=cut
  
#### END of POD documentation.
#-----------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;
use Gal;

my ($fi, $fo) = ('', '');
my ($fhi, $fho);
my ($fq, $ft) = ('') x 2; 
my $help_flag;

#----------------------------------- MAIN -----------------------------------#
GetOptions(
    "help|h"   => \$help_flag,
    "in|i=s"   => \$fi,
    "out|o=s"  => \$fo,
    "qry|q=s"  => \$fq,
    "tgt|t=s"  => \$ft,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;
pod2usage(2) if !$fq || !$ft;

if ($fi eq "stdin") {
    $fhi = \*STDIN;
} else {
    open ($fhi, $fi) || die "Can't open file $fi: $!\n";
}

if ($fo eq "stdout") {
    $fho = \*STDOUT;
} else {
    open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

my $tq = readTable(-in=>$fq, -header=>1);
my $hq = { map {$tq->elm($_, "id") => $tq->elm($_, "length")} (0..$tq->nofRow-1) };
my $tt = readTable(-in=>$ft, -header=>1);
my $ht = { map {$tt->elm($_, "id") => $tt->elm($_, "length")} (0..$tt->nofRow-1) };

gal2Psl($fhi, $fho, $hq, $ht);


__END__
