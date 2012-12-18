#!/usr/bin/perl -w
#根据指定参数某个group_id对应所有序列的snp统计信息
unshift (@INC, "/home/orion/Scripts/DB/INC");
require "get_input.pl";
require "function.pl";
use strict;
use Bio::Seq;
use Bio::SeqIO;

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
my @snp_id_arr;

open(OUT, ">out") or die("Cannot create output file\n");
if($input_mode == 1) {
	my @group_id_arr = split(",",$group_ids);
	@snp_id_arr = &group_to_snp(@group_id_arr);
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
	}
}

my $rCRS_path = "../../align/rCRS.gb";
my $rCRS_file = Bio::SeqIO->new(-file=>$rCRS_path, -format=>'genbank')
	or die("cannot open rCRS file\n");
my $rCRS = $rCRS_file->next_seq();

my %pos_snp_arr;
my @variant_arr;
my $count_snp_unique=0;
my $count_sites_unique=0;
my %count_variant;

for(my $i=0; $i<@snp_id_arr; $i++) {
	$query = "SELECT snp FROM mt_snp WHERE id=".$snp_id_arr[$i];
	$sth = $dbh->prepare($query);
	$sth->execute();
	my $rs = $sth->fetchrow_array();
	my $tmp = &reconstruct_snp_string($rs);
	&analyze($tmp);
	push (@variant_arr, $tmp);
}

foreach (sort {$a <=> $b} keys(%pos_snp_arr)) {
	my $pos = $_;
	my @snp_value_arr = split(" ",$pos_snp_arr{$_});
	my %tmp;
	foreach my $snp_value (@snp_value_arr) {
		if($snp_value) {
			if(exists($tmp{$snp_value})) {
				$tmp{$snp_value} ++;
			} else {
				$tmp{$snp_value} = 1;
			}
		}
	}
	$count_sites_unique ++;
	while(my ($a,$b) = each(%tmp)) {
		print OUT join("\t",$pos.$a,$b)."\t";
		$count_snp_unique ++;
		if(exists($count_variant{$b})) {
			$count_variant{$b} .= " ".$pos.$a;
		} else {
			$count_variant{$b} = $pos.$a;
		}
	}
	print OUT "\n";
}

print OUT "\n$count_snp_unique Unique SNPs at $count_sites_unique Unique Sites\n";
print OUT "\n";

foreach (sort {$a<=>$b} keys(%count_variant)) {
	my @tmp = split(" ",$count_variant{$_});
	print OUT $_."\t".@tmp."\t".$count_variant{$_}."\n";
}

sub analyze() {
	my ($snp_str) = @_;
	my @snp_arr = split(" ",$snp_str);
	foreach my $snp (@snp_arr) {
		if($snp =~ /^([.0-9]+)([A-Zd])$/) {
			my $pos = $1;
			my $value = $2;
			if(!exists($pos_snp_arr{$pos})) {
				$pos_snp_arr{$pos} = $value;
			} else {
				$pos_snp_arr{$pos} .= " ".$value;
			}
		}
	}
}