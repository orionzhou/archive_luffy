#!/usr/bin/perl -w
use strict;
use FindBin;
use lib "$FindBin::Bin";
use Common;
use File::Path qw/make_path remove_tree/;
use Time::HiRes qw/gettimeofday tv_interval/;
use POSIX qw/ceil floor/;
use List::Util qw/min max sum/;
use Bio::SeqIO;

my $fs = "/home/youngn/zhoup/Data/genome/pan3/11_genome.fa";
my $dir = "/home/youngn/zhoup/Data/genome/pan3/18_mappability";

make_path($dir) unless -d $dir;
my $f01 = "$dir/01.fa";
# winslide.pl -i 15.sizes -step 5 -size 60 | seqdbret.pl -d 11_genome.fa -o $dir/01.fa


