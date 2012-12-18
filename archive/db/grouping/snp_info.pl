#!/usr/bin/perl -w
#输入文件中的SNP进行评估
unshift (@INC, "/home/orion/Scripts/DB/INC");
require "db_connect.pl";
use strict;

our $dbh;
open(IN, "<snp") or die("Cannot read from in file\n");
open(OUT, ">out") or die("cannot create out file\n");
while(<IN>) {
	chomp;
	$_ =~ s/\s/\t/g;
	my @snp_arr = split("\t",$_);
	foreach my $snp (@snp_arr) {
		print OUT "===========================\n"
			.$snp;
		if($snp =~ /^([.0-9]+)([A-Zd])$/) {
			my ($pos,$value) = ($1,$2);
			my $query = "SELECT * FROM mt_snp_info WHERE pos=".int($pos);
			my $sth = $dbh->prepare($query);
			$sth->execute();
			my $flag = 0;
			while(my @rs = $sth->fetchrow_array()) {
				my $flag2 = 0;
				if($pos =~ /[\.]/) {
					if($rs[3] =~ /^([ATCG])\-\1[ATCG]*$value/) {
						$flag = 1;
						$flag2 = 1;
					}
				} elsif($value eq "d") {
					if($rs[3] =~ /^[ATCG]\-del/) {
						$flag = 1;
						$flag2 = 1;
					}
				} else {
					if($rs[3] =~ /^[ATCG]\-$value/) {
						$flag = 1;
						$flag2 = 1;
					}
				}
				if($flag2 == 1) {
					print OUT "\t".join("|",$rs[3],$rs[4],$rs[2])."\n";
				}
			}
			if($flag == 0) {
				print OUT "\t->N/A, located in:\n";
				
				my $query = "SELECT * FROM mt_locus_info WHERE (start<=".int($pos)." AND stop>=".int($pos).") OR (start>stop AND (stop>=".int($pos)." OR start<=".int($pos)."))";
				my $sth = $dbh->prepare($query);
				$sth->execute();
				while(my @rs = $sth->fetchrow_array()) {
					print OUT "\t".$rs[1]."(".$rs[2]."-".$rs[3].")"." : ".$rs[5]."\n";
				}
				
				print OUT "\t->Two nearest SNPs:\n";
				$query = "SELECT * FROM mt_snp_info WHERE pos<".int($pos)." ORDER BY pos DESC LIMIT 1";
				$sth = $dbh->prepare($query);
				$sth->execute();
				my @rs = $sth->fetchrow_array();
				print OUT "\t".join("\t",$rs[1],join("|",$rs[3],$rs[4],$rs[2]))."\n" if @rs;
				
				$query = "SELECT * FROM mt_snp_info WHERE pos>".int($pos)." ORDER BY pos ASC LIMIT 1";
				$sth = $dbh->prepare($query);
				$sth->execute();
				@rs = $sth->fetchrow_array();
				print OUT "\t".join("\t",$rs[1],join("|",$rs[3],$rs[4],$rs[2]))."\n" if @rs;
			}
		} else {
			print OUT "\terror snp string value\n";
			exit;
		}
	}
}