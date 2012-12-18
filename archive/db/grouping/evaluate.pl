#!/usr/bin/perl -w
#输入文件中的SNP进行评估
unshift (@INC, "/home/orion/Scripts/DB/INC");
require "get_input.pl";
require "function.pl";
use strict;

our $dbh;
my $query;
my $sth;

my %snp_arr;
my $snp_total = " ";
my @snp_total_arr;
my @dis_arr;
my @freq_arr;
my %value_snp;

open(IN, "<in") or die("Cannot read from in file\n");
open(OUT, ">out") or die("Cannot create output file\n");

our $range_include_string = &get_include_range();
our $range_exclude_string = "303-315";
our $ignore_chr_string = "N,R,Y,W";
while(<IN>) {
	chomp;
	$_ =~ s/\W+//g;
	if($_ && !exists($snp_arr{$_})) {
		my $snp_str = &get_snp_by_id($_);
		my $snp_str_reconstructed = &reconstruct_snp_string($snp_str);
		$snp_arr{$_} = $snp_str_reconstructed;
		#print $_."\t".$snp_str_reconstructed."\n";
		my @tmp = split(" ",$snp_str_reconstructed);
		foreach my $snp (@tmp) {
			if($snp) {
				if($snp_total !~ / $snp /) {
					$snp_total .= $snp." ";
					push (@snp_total_arr, $snp);
				}
			}
		}
	}
}
#print $snp_total."\n";
while(my ($snp_id,$snp_str) = each(%snp_arr)) {
		for(my $i=0; $i<@snp_total_arr; $i++) {
			if(!$dis_arr[$i]) {
				$dis_arr[$i] = " ";
			}
			if($snp_str =~ /$snp_total_arr[$i]/) {
				$dis_arr[$i] .= "1 ";
			} else {
				$dis_arr[$i] .= "0 ";
			}
		}
}

my $sample_count = scalar(keys(%snp_arr));
print OUT join("\t","snp","freq","value","distribution")."\n";
for(my $i=0; $i<@snp_total_arr; $i++) {
	my $dis = $dis_arr[$i];
	my $count = 0;
	while($dis =~ /1/g) {
		$count ++;
	}
	push (@freq_arr, $count);
	my $pic = ($count<$sample_count/2) ? $count : $sample_count-$count;
	if(!exists($value_snp{$pic})) {
		$value_snp{$pic} = "";
	}
	$value_snp{$pic} .= $i." ";
}
foreach my $value (sort {$b<=>$a} keys %value_snp) {
	my @snp_keys_arr = split(" ",$value_snp{$value});
	my @snp_pair_arr;
	for(my $i=0; $i<@snp_keys_arr; $i++) {
		for(my $j=$i+1; $j<@snp_keys_arr; $j++) {
			if(&judge_redundant($snp_keys_arr[$i],$snp_keys_arr[$j]) == 1) {
				push (@snp_pair_arr, $snp_keys_arr[$i]."=".$snp_keys_arr[$j]);
			}
		}
	}
	#print grouped snp keys array
	my @snp_grouped_keys_arr = &cluster_id_pair(\@snp_pair_arr);
	foreach my $group_key_string (@snp_grouped_keys_arr) {
		my @group_keys = split(" ",$group_key_string);
		my $tmp = "";
		foreach my $snp_key (@group_keys) {
			if($snp_key) {
				$tmp .= $snp_total_arr[$snp_key]." ";
			}
		}
		chop($tmp);
		print OUT join("\t",$tmp,$freq_arr[$group_keys[0]]."/".$sample_count,$value,$dis_arr[$group_keys[0]])."\n";
	}
	#print ungrouped snp keys array
	my $all_grouped_snp_keys = join(" ",@snp_grouped_keys_arr);
	foreach my $snp_key (@snp_keys_arr) {
		if($all_grouped_snp_keys !~ / $snp_key /) {
			print OUT join("\t",$snp_total_arr[$snp_key],$freq_arr[$snp_key]."/".$sample_count,$value,$dis_arr[$snp_key])."\n";
		}
	}
}

sub judge_redundant() {
	my ($key1, $key2) = @_;
	my $dis1 = $dis_arr[$key1];
	my $dis2 = $dis_arr[$key2];
	if($dis1 eq $dis2) {
		return 1;
	} else {
		$dis1 =~ s/0/2/g;
		$dis1 =~ s/1/0/g;
		$dis1 =~ s/2/1/g;
		if($dis1 eq $dis2) {
			return 1;
		} else {
			return 0;
		}
	}
}
