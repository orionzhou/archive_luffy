#!/usr/bin/perl -w
#将identical sequence的id进行扩展，成为相差任意个数snp的mtDNA数据库id list
unshift (@INC, "/home/orion/Scripts/DB/INC");
require "get_input.pl";
require "function.pl";
use strict;

our $dbh;
our $group_ids = &get_group_ids();
our $output_mode = &get_output_mode();
our $range_include_string = &get_include_range();
our $range_exclude_string = &get_exclude_range();
our $ignore_chr_string = &get_ignore_chrs();
our $diff_chr_string = &get_diff_chrs();

my @group_id_arr = split(",",$group_ids);
my @snp_id_arr = &group_to_snp(@group_id_arr);
my @variant_arr;
my @diff_selected = split(",",$diff_chr_string);
my $query;
my $sth;

foreach my $snp_id (@snp_id_arr) {
	$query = "SELECT snp FROM mt_snp WHERE id=$snp_id";
	$sth = $dbh->prepare($query);
	$sth->execute();
	my $rs = $sth->fetchrow_array();
	push (@variant_arr, &reconstruct_snp_string($rs));
}

open(IN, "<in_for_expand") or die("cannot read infile for expand\n");
open(OUT, ">in") or die("cannot write to outfile\n");

while(<IN>) {
	chomp;
	my @ele_arr = split(/[\*]/,$_);
	my $id = $ele_arr[1];
	$query = "SELECT snp FROM mt_snp WHERE id=$id";
	$sth = $dbh->prepare($query);
	$sth->execute();
	my $rs = $sth->fetchrow_array();
	my $ht_common = &reconstruct_snp_string($rs);
	for(my $a=0; $a<@snp_id_arr; $a++) {
		my $diff = &compare($ht_common,$variant_arr[$a]);
		foreach my $tmp (@diff_selected) {
			if($diff eq $tmp) {
				print OUT "*".sprintf("%03.0f",$snp_id_arr[$a]);
			}
		}
	}
	print OUT "*\n";
}