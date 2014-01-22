#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  galexpand.pl - expand a file of Gal[wide] format to Gal[long] format

=head1 SYNOPSIS
  
  galexpand.pl [-help] [-in input-file] [-out output-file]

  Options:
      -help   brief help message
      -in     input file
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

print $fho join("\t", @HEAD_GALL)."\n";

while( <$fhi> ) {
    chomp;
    next if /(^id)|(^\#)|(^\s*$)/;
    my $ps = [ split "\t" ];
    next unless @$ps == 19;
    my ($id, $qId, $qBeg, $qEnd, $qSrd, $qSize, $tId, $tBeg, $tEnd, $tSrd, $tSize,
        $match, $misMatch, $baseN, $ident, $e, $score, $qLocS, $tLocS) = @$ps;
    my ($rqloc, $rtloc) = (locStr2Ary($qLocS), locStr2Ary($tLocS));
    for my $i (0..@$rqloc-1) {
        my ($rqb, $rqe) = @{$rqloc->[$i]};
        my ($rtb, $rte) = @{$rtloc->[$i]};
        my ($qb, $qe) = $qSrd eq "-" ? ($qEnd-$rqe+1, $qEnd-$rqb+1)
            : ($qBeg+$rqb-1, $qBeg+$rqe-1);
        my ($tb, $te) = $tSrd eq "-" ? ($tEnd-$rte+1, $tEnd-$rtb+1)
            : ($tBeg+$rtb-1, $tBeg+$rte-1);
        print $fho join("\t", $id, $qId, $qb, $qe, $qSrd, $tId, $tb, $te, $tSrd)."\n";
    }
}
close $fhi;
close $fho;


__END__
