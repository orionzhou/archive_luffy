#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gtbsetdiff.pl - subtract one Gtb file (query) from another Gtb file (target)

=head1 SYNOPSIS
  
  gtb_compare.pl [-help] [-query query-Gtb] [-target target-Gtb] [-out subtracted-Gtb]

  Options:
      -help    brief help message
      -query   query Gtb file
      -target  target Gtb file
      -out     output (subtracted) Gtb file

=cut
  
#### END of POD documentation.
#-----------------------------------------------------------------------------
use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;

my ($fhq, $fht, $fho);
my ($fq, $ft, $fo) = ('') x 4;
my $help_flag;

#----------------------------------- MAIN -----------------------------------#
GetOptions(
    "help|h"  => \$help_flag,
    "query|qry|q=s"  => \$fq,
    "target|tgt|t=s" => \$ft,
    "out|o=s"  => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fq || !$ft || !$fo;

open ($fhq, "<$fq") || die "cannot open file $fq for reading: $!\n";
open ($fht, "<$ft") || die "cannot open file $ft for reading: $!\n";
open ($fho, ">$fo") || die "cannot open file $fo for writing: $!\n";

my $tq = readTable(-inh=>$fhq, -header=>1);
my $tt = readTable(-inh=>$fht, -header=>1);
close $fhq;
close $fht;

my $hq;
for my $i (0..$tq->nofRow-1) {
    my $str = join("#", map {$tq->elm($i, $_)} qw/chr srd locC/);
    $hq->{$str} = $i;
}

my $ht;
my @idxs_rm;
for my $i (0..$tt->nofRow-1) {
    my $str = join("#", map {$tt->elm($i, $_)} qw/chr srd locC/);
    if(exists($ht->{$str})) {
        push @idxs_rm, $ht->{$str};
    }
    $ht->{$str} = $i;
    if(exists($hq->{$str})) {
        push @idxs_rm, $i;
    }
}

my $no = $tt->nofRow;
$tt->delRows(\@idxs_rm);
print $fho $tt->tsv(1);
close $fho;

printf STDOUT "%5d / %5d substracted\n", scalar(@idxs_rm), $no;


__END__
