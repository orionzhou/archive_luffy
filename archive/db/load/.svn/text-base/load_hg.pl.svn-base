#!/usr/bin/perl -w
#根据mtDNAmanager对单倍组划分的结果(直接利用snp_id或者提取id字段中acc号作为唯一标识符)更新mt_snp的hg_exp和hg_est
unshift (@INC, "/home/orion/Scripts/DB/INC");
require "db_connect.pl";
use strict;

our $dbh;
my $working_dir = "./";
my $in_file = "in_hg";

open(IN, $working_dir.$in_file) or die("cannot open in_hg file\n");

print "Which identifier provided?\n"
	."\t1:\tmt_snp.id\n"
	."\t2:\tAccession Number\n"
	."(default 1)";
my $id_type = <STDIN>;
chomp($id_type);
$id_type = ($id_type) ? $id_type:1;
die("illegal input\n") if $id_type !~ /^[12]$/;

while(<IN>) {
	chomp;
	if($_) {
		my ($id,$hg_exp,$hg_est) = split("\t",$_);
		$hg_exp = $hg_exp ? $hg_exp : "ND";
		$hg_est = $hg_est ? $hg_est : "ND";
		my $snp_id;
		if($id_type == 1) {
			my @ary = split("_",$id);
			$snp_id = $ary[1];
		} elsif($id_type == 2) {
			$id =~ /\_([.A-Z0-9]+)$/;
			my $acc = $1;
			die("cannot determine acc no. from $id\n") if !$acc;
			my $query = "SELECT b.id FROM mt_main AS a, mt_snp AS b WHERE a.snp_id=b.id AND a.acc=\"$acc\"";
			my $sth = $dbh->prepare($query);
			$sth->execute();
			my ($matrix_ref);
			$matrix_ref = $sth->fetchall_arrayref();
			my ($rows) = ( !defined($matrix_ref) ? 0 : scalar(@{$matrix_ref}) );
			if($rows != 1) {
				die("$id reflected to 0 or >1 records in db\n");
			} else {
				$snp_id = $matrix_ref->[0][0];
			}
		}
		my $query = "UPDATE mt_snp SET hg_exp=\"$hg_exp\", hg_est=\"$hg_est\" WHERE id=$snp_id";
		my $rs = $dbh->do($query);
		if(!$rs){
			die("Update mt_snp failed at $id\n"); 
		} else {
			print $id." - updated\n";
		}
	}
}