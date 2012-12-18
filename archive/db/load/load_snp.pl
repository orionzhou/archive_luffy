#!/usr/bin/perl -w
#根据mt_group中某(几)条记录核查并查找其对应的mt_main表中的序列，对snp_id字段为空的记录进行bl2seq，并将结果存储到mt_snp表中，同时更新mt_group表中snp_id字段
unshift (@INC, "/home/orion/Scripts/DB/INC");
require "db_connect.pl";
require "/home/orion/Scripts/align/al2seq_db.pl";
use strict;
use Bio::Seq;
use Bio::SeqIO;

my @group_id_arr = (1);
our $dbh;

my $i=0;
foreach my $group_id (@group_id_arr) {
	my @id_arr = &check_id($group_id);
	foreach my $id (@id_arr) {
		my $query = "select * from mt_main where id=$id";
		my $sth = $dbh->prepare($query);
		$sth->execute();
		my $hashref = $sth->fetchrow_hashref();
		my $acc = $hashref->{"acc"};
		my $snp_id = $hashref->{"snp_id"};
		$i ++;
		if($snp_id == 0) {
			#do bl2seq
			my $seqobj = Bio::Seq->new( -display_id => sprintf("%03.0f",$i)."_$acc", 
				-seq => $hashref->{"seq"} );
			my $rCRS_path = "../../align/rCRS.gb";
			my $rCRS_file = Bio::SeqIO->new(-file => $rCRS_path, -format=>'genbank');
			my $rCRS = $rCRS_file->next_seq();
			my @rs_arr = split("\t",&al2seq($seqobj, $rCRS));
			#insert result to mt_snp and update mt_main.snp_id
			$query = "insert into mt_snp (start,stop,snp) VALUES (\"$rs_arr[0]\", \"$rs_arr[1]\", \"$rs_arr[2]\")";
			my $rs = $dbh->do($query);
			if(!defined($rs)) {
				die("error in inserting into mt_snp\n");
			} else {
				$query = "select LAST_INSERT_ID()";
				$sth = $dbh->prepare($query);
				$sth->execute();
				my $insert_snp_id = $sth->fetchrow_array();
				$query = "update mt_main set snp_id=$insert_snp_id where id=$id";
				$rs = $dbh->do($query);
				die("error in updating mt_main.snp_id\n") if !$rs;
				print "seq_id($id) acc($acc) -> snp_id($insert_snp_id)\n";
			}
		} else {
			print ("seq_id $id - $acc has already been aligned\n");
		}
	}
}

sub check_id() {
	my ($group_id) = @_;
	my @id_arr;
	my $query = "select id from mt_main where group_id=$group_id";
	my $sth = $dbh->prepare($query);
	$sth->execute();
	while(my $id = $sth->fetchrow_array()) {
		push (@id_arr,$id);
	}
	$query = "select seq_id from mt_group where id=$group_id";
	$sth = $dbh->prepare($query);
	$sth->execute();
	my $count = 0;
	my $id_str = $sth->fetchrow_array();
	my @id_range_arr = split(",",$id_str);
	foreach my $id_range (@id_range_arr) {
		my @id_pair = split("-",$id_range);
		if(@id_pair == 1) {
			$count += 1;
		} elsif(@id_pair ==2) {
			$count += ($id_pair[1]-$id_pair[0]+1);
		} else {
			die("group(id:$group_id) - seq_id string format error\n");
		}
	}
	if(@id_arr == $count) {
		return @id_arr;
	} else {
		die("group(id:$group_id) seq_id check error:\t".@id_arr."(db) != $count\n");
		return 1;
	}
}