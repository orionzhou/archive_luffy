#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  qsub.blat.pl - split a fasta file and build qsub commands

=head1 SYNOPSIS
  
  qsub.seq.pl [-help] [-in input-file] [-out output-directory]
                      [-n number-batches] [-db blat-db]

  Options:
    -h (--help)   brief help message
    -i (--in)     input (fasta) file
    -o (--out)    output directory
    -n (--num)    number of qsub batches (def: 1)
    -t (--tgt)    blat target (def: "HM101")
    -g (--tag)    qsub job tag (def: "pz")

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin";
use File::Path qw/make_path remove_tree/;
use File::Spec;
use Common;
use Data::Dumper;

my ($fi, $dir) = ('') x 2;
my ($n, $tgt) = (1, "HM101");
my $tag = "pz";
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "in|i=s"  => \$fi,
  "out|o=s" => \$dir,
  "num|n=i" => \$n,
  "tgt|t=s" => \$tgt,
  "tag|g=s" => \$tag,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$dir;

$fi = File::Spec->rel2abs($fi);
$dir = File::Spec->rel2abs($dir);

runCmd("rm -rf $dir") if -d $dir;
-d $dir || make_path($dir);
chdir $dir || die "cannot chdir to $dir\n";

runCmd("ln -sf $fi part.fas");

my $part = $n * 16;
my $digits = getDigits($part);
runCmd("pyfasta split -n $part part.fas");
runCmd("rm part.fas part.fas.*");

my $ps = runCmd("du -h part.*.fas", 2);
my @sizes;
for (@$ps) { push @sizes, [ split " ", $_ ]->[0]; }
@sizes = sort @sizes;
print "\n##### stats begin #####\n";
printf "range: %s  -  %s\n", $sizes[0], $sizes[$#sizes];
print "##### stats end   #####\n\n";

print "\n##### qsub command begins #####\n";
for my $i (0..$n-1) {
  my $beg = $i * 16;
  print "qsub blat -N blat.$tag.$i -v PRE=$dir/part,SUF=fas,BEG=$beg,DIG=$digits,TGT=$tgt -l qos=weightlessqos\n";
}
print "##### qsub command ends   #####\n\n";



exit 0;
