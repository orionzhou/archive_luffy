#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  galTiling.pl - coordinate tiling of a Gal file

=head1 SYNOPSIS
  
  galTiling.pl [-help] [-in input-file] [-out output-file]

  Options:
      -help   brief help message
      -in     input file - needs to be sorted by qId
      -out    output file

=cut
  
#### END of POD documentation.
#-----------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Location;
use Gal;

my ($fi, $fo) = ('') x 2;
my ($fhi, $fho);
my $help_flag;

#----------------------------------- MAIN -----------------------------------#
GetOptions(
    "help|h"   => \$help_flag,
    "in|i=s"   => \$fi,
    "out|o=s"  => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

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

print $fho join("\t", qw/id qId qBeg qEnd qSrd qLen tId tBeg tEnd tSrd tLen
    match misMatch baseN ident e score qLoc tLoc/)."\n";

my @rows;
my $tag = '';
while(<$fhi>) {
    chomp;
    next if /(^id)|(^\#)|(^\s*$)/;
    my @ps = split "\t";
    $ps[-2] = locStr2Ary($ps[-2]);
    $ps[-1] = locStr2Ary($ps[-1]);
    if($ps[1] ne $tag && $tag ne "") {
        gal_tiling(\@rows, $fho);
        @rows = ();
    }
    push @rows, \@ps;
    $tag = $ps[1];
}
gal_tiling(\@rows, $fho) if @rows > 0;
close $fhi;
close $fho;


__END__
