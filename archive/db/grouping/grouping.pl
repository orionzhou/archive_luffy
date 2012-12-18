#!/usr/bin/perl -w
#根据给定的一组SNP对in中的snp_id进行dissection
unshift (@INC, "/home/orion/Scripts/DB/INC");
require "get_input.pl";
require "function.pl";
use strict;

our $dbh;
my $query;
my $sth;

my %snp_arr;

open(IN, "<in") or die("Cannot read from in file\n");
open(SNP, "<snp") or die("Cannot read snp file\n");
open(OUT, ">out") or die("cannot create out file\n");

our $range_include_string = &get_include_range();
our $range_exclude_string = "303-315";
our $ignore_chr_string = "N,R,Y,W";
while(<IN>) {
	chomp;
	$_ =~ s/\W+//g;
	if($_ && !exists($snp_arr{$_})) {
		my $snp_str = &get_snp_by_id($_);
		$snp_arr{$_} = &reconstruct_snp_string($snp_str);
	}
}

my @group_arr = (join("-",keys(%snp_arr)));
print OUT join("\n",@group_arr)."\n";
print OUT scalar(keys(%snp_arr))."\n";
while(<SNP>) {
	chomp;
	$_ =~ s/\s/\t/g;
	$_ =~ s/[^.\d\tA-Zd]//g;
	if($_) {
		my @snp_in_line = split("\t",$_);
		#for(my $i=0;$i<scalar(@snp_in_line);$i++) {
		for(my $i=0;$i<1;$i++) {
			print OUT $snp_in_line[$i]."\n";
			my @group_arr_2;
			foreach my $group_string (@group_arr) {
				my @ids = split("-",$group_string);
				my ($group_string_1, $group_string_2) = ("", "");
				foreach my $snp_id (@ids) {
					if($snp_arr{$snp_id} =~ /$snp_in_line[$i]/) {
						$group_string_1 .= $snp_id."-";
					} else {
						$group_string_2 .= $snp_id."-";
					}
				}
				if($group_string_1 ne "") {
					chop($group_string_1);
					push (@group_arr_2,$group_string_1);
				}
				if($group_string_2 ne "") {
					chop($group_string_2);
					push (@group_arr_2,$group_string_2);
				}
			}
			@group_arr = @group_arr_2;
			my @scalar_arr;
			foreach my $group_string (@group_arr) {
				my @tmp = split("-",$group_string);
				push (@scalar_arr, scalar(@tmp));
			}
			print OUT join("\t",@group_arr)."\n";
			print OUT join("\t",@scalar_arr)."\n";
		}
	}
}