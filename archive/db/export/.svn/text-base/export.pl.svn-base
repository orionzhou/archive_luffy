#!/usr/bin/perl -w
#根据指定参数输出某种格式的mtSNP数据文件
unshift (@INC, "/home/orion/Scripts/DB/INC");
require "get_input.pl";
require "function.pl";
use strict;

our $dbh;
our $input_mode = &get_input_mode();
our $group_ids;
our $output_mode;
if($input_mode == 1) {
	$group_ids = &get_group_ids();
	$output_mode = &get_output_mode();
}
our $range_include_string = &get_include_range();
our $range_exclude_string = &get_exclude_range();
our $ignore_chr_string = &get_ignore_chrs();

my $query;
my $sth;

open(OUT, ">out") or die("Cannot create output file\n");
if($input_mode == 1) {
	my @group_id_arr = split(",",$group_ids);
	my @snp_id_arr = &group_to_snp(@group_id_arr);
	&output_normal(\@snp_id_arr);
} elsif($input_mode == 2) {
	open(IN, "<in") or die("cannot read input file\n");
	while(<IN>) {
		chomp;
		my @snp_id_arr;
		my @ary_tmp = split(/[\*]/,$_);
		foreach my $snp_id_str (@ary_tmp) {
			if($snp_id_str) {
				push (@snp_id_arr,int($snp_id_str));
			}
		}
		&output_extra(\@snp_id_arr);
		print OUT "\n";
	}
}

sub output_normal() {
	my ($id_arr_ref) = @_;
	our $range_include_string;
	my @range_arr = split("-",$range_include_string);
	my $i = 0;
	print OUT join("\t", "id", "start", "stop", "SNPs")."\n";
	while(my $snp_id = $id_arr_ref->[$i]) {
		$i ++;
		my $new_id = sprintf("%03.0f",$i)."_".$snp_id;
		my $query = "SELECT snp FROM mt_snp WHERE id=$snp_id";
		my $sth = $dbh->prepare($query);
		$sth->execute();
		my $rs = $sth->fetchrow_array();
		print OUT join("\t", $new_id, $range_arr[0], $range_arr[scalar(@range_arr)-1], &reconstruct_snp_string($rs))."\n";
	}
}

sub output_extra() {
	my ($id_arr_ref) = @_;
	our $range_include_string;
	my @range_arr = split("-",$range_include_string);
	my $i = 0;
	print OUT join("\t","id","SNP","start","stop","hg_exp","hg_est")."\n";
	while(my $snp_id = $id_arr_ref->[$i]) {
		$i ++;
		my $new_id = sprintf("%03.0f",$i)."_".$snp_id;
		$query = "SELECT hg_exp,hg_est,snp FROM mt_snp WHERE id=$snp_id";
		$sth = $dbh->prepare($query);
		$sth->execute();
		my @rs = $sth->fetchrow_array();
		print OUT join("\t", $new_id, &reconstruct_snp_string($rs[2]), $range_arr[0], $range_arr[scalar(@range_arr)-1], $rs[0], $rs[1])."\n";
	}
}