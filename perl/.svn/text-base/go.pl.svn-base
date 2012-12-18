#!/usr/bin/perl -w
use strict; use Init; use Common; use Localdb; use Readfile; use Writefile;
use Path::Class; use Time::HiRes qw/gettimeofday tv_interval/; use Data::Dumper;
use List::Util qw/min max sum/; use POSIX qw/ceil floor/; use Math::Round qw/round/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;

my $dirW = dir($DIR_Misc1, "cat");
my $f01 = file($dirW, "01_mt_30.txt");
my $f02 = file($dirW, "02_mt_30_filtered.txt");
#filter($f01, $f02);

