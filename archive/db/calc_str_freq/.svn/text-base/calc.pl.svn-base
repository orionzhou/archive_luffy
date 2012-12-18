#!/usr/bin/perl -w
require "../db_connect.pl";
use strict;
use Switch;

our $dbh;
my $working_dir = "./";
my $in_file = "in";
my $out_file = "out";

open(IN, $working_dir.$in_file) or die("Couldn't read in file");
open(OUT, ">".$working_dir.$out_file) or die("Couldn't create out file");

my @option = ("1-双纯合相同","2-双杂合相同","3-杂纯共一","4-杂杂共一","5-其他");  

print OUT join("\t","option:",@option,"\n");
print OUT join("\t","id","indi_pair","option|result","\n");
while(<IN>) {
	chomp;
	if($_) {
		my @line_arr;
		push (@line_arr,$_);
		for(my $i=2;$i<=4;$i++) {
			my $tmp = <IN>;
			chomp($tmp);
			push (@line_arr,$tmp);
		}
		die("Incomplete dataset at ".join("\n",@line_arr)."\n") if @line_arr != 4;
		my @ele_arr1 = split("\t",$line_arr[0]);
		my @ele_arr2 = split("\t",$line_arr[1]);
		my @ele_arr3 = split("\t",$line_arr[2]);
		my @ele_arr4 = split("\t",$line_arr[3]);
		if(@ele_arr1!=17 || @ele_arr2!=17 || @ele_arr3!=17 ||@ele_arr4!=17) {
			die("Incomplete dataset at ".join("\n",@line_arr)."\n");
		}
		print OUT join("\t",$ele_arr1[0],$ele_arr1[1]."#".$ele_arr3[1],"\t");
		for(my $i=1;$i<=15;$i++) {
			my $j = $i+1;
			my @allele_arr = ($ele_arr1[$j],$ele_arr2[$j],$ele_arr3[$j],$ele_arr4[$j]);
			print OUT &calc($i,@allele_arr)."\t";
		}
		print OUT "\n";
	}
}

sub calc() {
	my ($locus_id,@allele_arr) = @_;
	my @freq_arr;
	for(my $i=0;$i<=3;$i++) {
		my $query = "select freq from str_freq where locus_id=\"".$locus_id."\" AND allele=\"".$allele_arr[$i]."\"";
		my $sth = $dbh->prepare($query);
		$sth->execute();
		my $freq = $sth->fetchrow_array();
		$sth->finish();
		die("unrecorded allele ".$allele_arr[$i]." at locus $locus_id\n") if !defined($freq);
		push(@freq_arr,$freq);
	}
	my $option;
	my $result;
	if($allele_arr[0]==$allele_arr[1] && $allele_arr[2]==$allele_arr[3] && $allele_arr[0]==$allele_arr[3]) {
		$option = 1;
		$result = &formula_1($freq_arr[0]);
	} elsif( 
				($allele_arr[0]!=$allele_arr[1] && $allele_arr[2]!=$allele_arr[3]) 
				&& 
				(
					($allele_arr[0]==$allele_arr[2] && $allele_arr[1]==$allele_arr[3]) 
					||
					($allele_arr[0]==$allele_arr[3] && $allele_arr[1]==$allele_arr[2])
				)
			) {
		$option = 2;
		$result = &formula_2($freq_arr[0],$freq_arr[1]);
	} elsif(
				(
					($allele_arr[0]==$allele_arr[1] && $allele_arr[2]!=$allele_arr[3]) 
					&&
					($allele_arr[0]==$allele_arr[2] || $allele_arr[0]==$allele_arr[3])
				)
				||
				(
					($allele_arr[0]!=$allele_arr[1] && $allele_arr[2]==$allele_arr[3])
					&&
					($allele_arr[0]==$allele_arr[2] || $allele_arr[1]==$allele_arr[2])
				)
			) {
		$option = 3;
		if($allele_arr[0]==$allele_arr[1] && $allele_arr[2]!=$allele_arr[3]) {
			$result = &formula_3($freq_arr[0]);
		} else {
			$result = &formula_3($freq_arr[2]);
		}
	} elsif(
				($allele_arr[0]!=$allele_arr[1] || $allele_arr[2]!=$allele_arr[3])
				&&
				(
					($allele_arr[0]==$allele_arr[2] && $allele_arr[1]!=$allele_arr[3])
					||
					($allele_arr[0]==$allele_arr[3] && $allele_arr[1]!=$allele_arr[2])
					||
					($allele_arr[0]!=$allele_arr[2] && $allele_arr[1]==$allele_arr[3])
					||
					($allele_arr[0]!=$allele_arr[3] && $allele_arr[1]==$allele_arr[2])
				)
			) {
		$option = 4;
		if($allele_arr[0]==$allele_arr[2] && $allele_arr[1]!=$allele_arr[3]) {
			$result = &formula_3($freq_arr[0]);
		} elsif($allele_arr[0]==$allele_arr[3] && $allele_arr[1]!=$allele_arr[2]) {
			$result = &formula_3($freq_arr[0]);
		} elsif($allele_arr[0]!=$allele_arr[2] && $allele_arr[1]==$allele_arr[3]) {
			$result = &formula_3($freq_arr[1]);
		}else {
			$result = &formula_3($freq_arr[1]);
		}
	} else {
		$option = 5;
		$result = &formula_5();
	}
	return $option."|".$result;
}

sub formula_1() {
	my ($pi) = @_;
	if($pi==0) {
		return "Undefined";
	} else {
		return sprintf("%5f",(($pi+1)*($pi+1))/(4*$pi*$pi));
	}
}

sub formula_2() {
	my ($pi,$pj) = @_;
	if($pi==0 || $pj==0) {
		return "Undefined";
	} else {
		return sprintf("%5f",(2*$pi*$pj+$pi+$pj+1)/(8*$pi*$pj));
	}
}

sub formula_3() {
	my ($pi) = @_;
	if($pi==0) {
		return "Undefined";
	} else {
		return sprintf("%5f",($pi+1)/(4*$pi));
	}
}

sub formula_4() {
	my ($pi) = @_;
	if($pi==0) {
		return "Undefined";
	} else {
		return sprintf("%5f",(2*$pi+1)/(8*$pi));
	}
}

sub formula_5() {
	return sprintf("%5f",0.25);
}