#!/usr/bin/perl -w
#统计输入文件中不同单倍组的频率
use strict;

open(IN, "<in") or die("cannot open infile\n");
#open(OUT, ">out") or die("cannot create outfile\n");

my $buffer = " ";
my %hg;
while(<IN>) {
	chomp;
	if($_) {
		my @ele_arr = split(" ",$_);
		if($buffer !~ / $ele_arr[0] /) {
			$buffer .= $ele_arr[0]." ";
			$hg{$ele_arr[0]} = 1;
		} else {
			$hg{$ele_arr[0]} += 1;
		}
	}
}

while(my ($key,$value) = each(%hg)) {
	print $key."\t".$value."\n";
}