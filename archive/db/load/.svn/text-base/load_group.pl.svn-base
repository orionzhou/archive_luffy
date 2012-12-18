#!/usr/bin/perl -w
#根据group文件的内容，对mt_group表进行选择性插入，对mt_main的group_id字段进行更新
unshift (@INC, "INC") && require "db_connect.pl";
use strict;
our $dbh;
my $working_dir = "./";
my $group_file = "in_group";

open(GROUP, $working_dir.$group_file) or die("Couldn't read group file");
while(<GROUP>) {
	chomp;
	my $line = $_;
	my @ele_arr = split("\t",$line);
	my $region = shift(@ele_arr);
	my $name = shift(@ele_arr);
	my $acc_str = (@ele_arr) ? shift(@ele_arr):"";
	my $pmid = (@ele_arr) ? shift(@ele_arr):"";
	my $title = (@ele_arr) ? shift(@ele_arr):"";
	my $query = "select count(*) from mt_group where pmid=\"$pmid\"";
	my $sth = $dbh->prepare($query);
	$sth->execute();
	my $count = $sth->fetchrow_array();
	$sth->finish();
	if(!$count) {
		#do insert
		my @acc_arr = &acc_construct($acc_str);
		my $query = "insert into mt_group (region,name,pmid,title) VALUES ($region, \"$name\", \"$pmid\", \"$title\")";
		#print $query."\n".@acc_arr."\n";
		my $rows = $dbh->do($query);
		if(!defined($rows)) {
			die("Error in inserting $acc_str\n");
		} else {
			$query = "select LAST_INSERT_ID()";
			$sth = $dbh->prepare($query);
			$sth->execute();
			my $insert_id = $sth->fetchrow_array();
			print "$insert_id-$region-$name-$pmid\t".substr($title,0,20)."...\tinserted\n";
			foreach my $acc (@acc_arr) {
				my $query = "update mt_main set group_id=".$insert_id." where acc=\"$acc\"";
				my $rs = $dbh->do($query);
				die("error while updating mt_main with $acc in $pmid\n") if !$rs;
			}
		}
	} else {
		print "already exists $pmid\t".substr($title,0,20)."...\n";
	}

}

sub acc_construct() {
	my ($acc) = @_;
	my @acc_str_arr = split(",",$acc);
	my @acc_arr;
	foreach my $acc_str (@acc_str_arr) {
		my @pair = split("-",$acc_str);
		if(@pair==1) {
			push (@acc_arr,$acc_str);
		} elsif(@pair==2) {
			$pair[0] =~ /\b([A-Z]{1,3})(\d{4,6}\b)/i;
			my $prefix = $1;
			my $start = $2;
			my $surfix_len = length($start);
			$pair[1] =~ /\b([A-Z]{1,3})(\d{4,6}\b)/i;
			if($1 ne $prefix) {
				die("Error in processing ".$acc_str."\n");
			}
			my $stop = $2;
			for(my $i=$start;$i<=$stop;$i++) {
				push (@acc_arr, $prefix.sprintf("%0".$surfix_len."d",$i));
			}
		} else {
			die("Error in processing ".$acc_str."\n");
		}
	}
	return @acc_arr;
}