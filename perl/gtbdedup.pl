#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gtbdedup.pl - remove redundant gene models from a Gtb file (based on CDS locations)

=head1 SYNOPSIS
  
  gtbdedup.pl [-help] [-in input-Gtb] [-out output-Gtb]

  Options:
      -help   brief help message
      -in     input Gtb
      -out    output Gtb

=cut
  
#### END of POD documentation.
#-----------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;
use Location;

my ($fi, $fo) = ('') x 2;
my $help_flag;

#----------------------------------- MAIN -----------------------------------#
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
my @idxs;
my $h;
for my $i (0..$t->lastRow) {
    my ($id, $par, $chr, $beg, $end, $srd, $locE, $locI, $locC, $loc5, $loc3, $phase, $src, $conf, $cat1, $cat2, $cat3, $note) = $t->row($i);
    my $rloc = [sort {$a->[0] <=> $b->[0]} @{locStr2Ary($locC)}];
    my $loc = $srd eq "-" ? [map {[$end-$_->[1]+1, $end-$_->[0]+1]} @$rloc] : 
        [map {[$beg+$_->[0]-1, $beg+$_->[1]-1]} @$rloc];
    my $str = join("|", $chr, locAry2Str($loc));
    my $len = $end - $beg + 1;
    if(exists $h->{$str}) {
        my ($pidx, $plen) = @{$h->{$str}};
        if($plen <= $len) {
            push @idxs, $i;
        } else {
            push @idxs, $pidx;
            $h->{$str} = [$i, $len];
        }
    } else {
        $h->{$str} = [$i, $len];
    }
}
printf "%5d | %5d removed\n", scalar(@idxs), $t->nofRow;

$t->delRows(\@idxs);
print $fho $t->tsv(1);
close $fho;



__END__
