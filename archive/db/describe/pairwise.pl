#!/usr/bin/perl -w
#根据指定参数某个group_id对应所有序列的snp统计信息
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
our $diff_chr_string = &get_diff_chrs();

my @snp_id_arr;
my @diff_selected = split(",",$diff_chr_string);

open(OUT, ">out") or die("Cannot create output file\n");
if($input_mode == 1) {
	my @group_id_arr = split(",",$group_ids);
	@snp_id_arr = &group_to_snp(@group_id_arr);
} elsif($input_mode == 2) {
	open(IN, "<in") or die("cannot read input file\n");
	while(<IN>) {
		chomp;
		my @ary_tmp = split(/[\W]/,$_);
		foreach my $snp_id_str (@ary_tmp) {
			if($snp_id_str) {
				push (@snp_id_arr,int($snp_id_str));
			}
		}
	}
}

my @variant_arr;
my %diff_count;
my @seq_pair_identical;

foreach my $snp_id (@snp_id_arr) {
	my $snp_str = &get_snp_by_id($snp_id);
	push (@variant_arr, &reconstruct_snp_string($snp_str));
}

for(my $a=0; $a<@snp_id_arr; $a++) {
	for(my $b=0; $b<=$a; $b++) {
		#print OUT "\t";
	}
	for(my $b=$a+1; $b<@snp_id_arr; $b++) {
		my $diff = &compare($variant_arr[$a],$variant_arr[$b]);
		if(exists($diff_count{$diff})) {
			$diff_count{$diff} ++;
		} else {
			$diff_count{$diff} = 1;
		}
		#print OUT $diff."\t";
		for(my $i=0; $i<@diff_selected; $i++) {
			if($diff == $diff_selected[$i]) {
				push (@seq_pair_identical, sprintf("%03.0f",$snp_id_arr[$a])."=".sprintf("%03.0f",$snp_id_arr[$b]));
			}
		}
	}
	#print OUT "\n";
}

#print OUT "\n";
print OUT join("\t","Pairwise Difference(s)","Count")."\n";
foreach my $diff (sort {$a<=>$b} keys %diff_count) {
	print OUT join("\t",$diff,$diff_count{$diff})."\n";
}

my $buffer = " ";
foreach my $seq_pair (@seq_pair_identical) {
	my ($g,$h) = split("=",$seq_pair);
	if($buffer=~/ $g / || $buffer=~/ $h /) {
		if($buffer!~/ $g /) {
			$buffer .= "$g ";
		}
		if($buffer!~/ $h /) {
			$buffer .= "$h ";
		}
	} else {
		$buffer .= "$g $h ";
	}
}
my @id_array = split(" ",$buffer);
print OUT "\n".@id_array." Sequences are undistinguishable due to (".$diff_chr_string.") differences\n";

if($diff_chr_string eq "0") {
	&identical_seq_distribution(@seq_pair_identical);
}

sub identical_seq_distribution() {
	#注意split内的分隔符支持正则表达式
	my @seq_pair_arr = @_;
	my @cluster;
	print OUT "\n";
	print OUT "Distribution of Identical Sequences:\n";
	foreach my $seq_str (@seq_pair_arr) {
		my ($seq1,$seq2) = split("=",$seq_str);
		my $flag = 0;
		for(my $i=0; $i<scalar(@cluster); $i++) {
			if($cluster[$i] =~ /\*$seq1\*/ || $cluster[$i] =~ /\*$seq2\*/) {
				if($cluster[$i] !~ /\*$seq1\*/) {
					$cluster[$i] .= $seq1."*";
				}
				if($cluster[$i] !~ /\*$seq2\*/) {
					$cluster[$i] .= $seq2."*";
				}
				$flag = 1;
			}
		}
		if($flag == 0) {
			push (@cluster, "*".$seq1."*".$seq2."*");
		}
	}
	my %buff;
	foreach my $tmp (@cluster) {
		my @ary = split(/[\*]/,$tmp);
		if(exists($buff{@ary-1})) {
			$buff{@ary-1} ++;
		} else {
			$buff{@ary-1} = 1;
		}
		print OUT (@ary-1)."\t".$tmp."\n";
	}
	print OUT "\n";
	
	print OUT join("\t","Size of Cluster","Count")."\n";
	foreach my $tmp (sort {$a<=>$b} keys %buff) {
		print OUT join("\t",$tmp,$buff{$tmp})."\n";
	}
	print OUT "\n";
}