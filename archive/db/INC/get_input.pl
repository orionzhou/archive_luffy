#!/usr/bin/perl -w
#根据用户输入选项作为其他脚本的运行参数
require "db_connect.pl";
use strict;

our $dbh;

sub get_input_mode() {
	my %input_mode_arr = (
		1 => "Manually Select group id(s) existing in current db",
		2 => "From input file (./in)"
		);
	print "=============================================\n";
	foreach my $key (sort {$a<=>$b} keys %input_mode_arr) {
		print "\t".join("\t",$key.":",$input_mode_arr{$key})."\n";
	}
	print "=============================================\n"
		."how do you want to do the export(default 1):";
	my $input_mode = <STDIN>;
	chomp($input_mode);
	$input_mode = ($input_mode eq "") ? 1:$input_mode;
	if(!$input_mode_arr{$input_mode}) {
		die("illegal input mode\n");
	}
	return $input_mode;
}

sub get_group_ids() {
	my $group_id_str = "*";
	my $query = "SELECT id,name,pmid,title FROM mt_group";
	my $sth = $dbh->prepare($query);
	$sth->execute();
	print "=============================================\n"
		."\t".join("\t","id","pubmed_id","title")."\n";
	while(my $hashref = $sth->fetchrow_hashref()) {
		my $tmp1 = $hashref->{"pmid"} ? $hashref->{"pmid"} : "";
		my $tmp2 = $hashref->{"title"} ? substr($hashref->{"title"},0,25)."..." : "";
		print "\t".join("\t",$hashref->{"id"},$hashref->{"name"},$tmp1,$tmp2,"\n");
		$group_id_str .= $hashref->{"id"}."*";
	}
	print "=============================================\n"
		."enter the group id(s) you wish to export\n"
		."(divided by \",\"):";
	my $id_str = <STDIN>;
	chomp($id_str);
	$id_str =~ s/[^\d\,]//g;
	die("error input\n") if $id_str eq "";
	my @group_id_arr = split(",",$id_str);

	#check group_id(s)
	foreach my $group_id (@group_id_arr) {
		if($group_id_str !~ /\*$group_id\*/) {
			die("group_id($group_id) do not exist in db\n");
		}
	}
	return $id_str;
}

sub get_output_mode() {
	my %output_mode_arr = (
		1 => "from sequence data (table mt_main)",
		2 => "from mt_group.snp_id"
		);
	print "=============================================\n";
	foreach my $key (sort {$a<=>$b} keys %output_mode_arr) {
		print "\t".join("\t",$key.":",$output_mode_arr{$key})."\n";
	}
	print "=============================================\n"
		."select an output mode(default 1):";
	my $output_mode = <STDIN>;
	chomp($output_mode);
	$output_mode = ($output_mode eq "") ? 1:$output_mode;
	if(!$output_mode_arr{$output_mode}) {
		die("illegal output mode\n");
	}
	return $output_mode;
}

sub get_include_range() {
	my %range_include = (
		1, ["mtGenome","1-16569"], 
		2, ["control region","16024-16569,1-437,438-576"],
		3, ["HVRI+HVRII","16024-16365,73-340"],
		4, ["NonHVRI/HVRII","1-72,341-16023,16366-16569"],
		5, ["NADH dehydrogenase subunit 1","3307-4262"],
		6, ["NADH dehydrogenase subunit 2","4470-5511"],
		7, ["Cytochrome C oxidase subunit I","5904-7445"],
		8, ["Cytochrome C oxidase subunit II","7586-8269"],
		9, ["ATP synthase F0 subunit 8","8366-8572"],
		10, ["ATP synthase F0 subunit 6","8527-9207"],
		11, ["Cytochrome C oxidase subunit III","9207-9990"],
		12, ["NADH dehydrogenase subunit 3","10059-10404"],
		13, ["NADH dehydrogenase subunit 4L","10470-10766"],
		14, ["NADH dehydrogenase subunit 4","10760-12137"],
		15, ["NADH dehydrogenase subunit 5","12337-14148"],
		16, ["NADH dehydrogenase subunit 6","14149-14673"],
		17, ["Cytochrome b","14147-15887"],
		18, ["12s rRNA","648-1601"],
		19, ["16s rRNA","1671-3229"],
		20, ["All tRNA's combined","577-647,1602-1670,3230-3304,4263-4331,4329-4400,"
			."4402-4469,5512-5579,5587-5655,5657-5729,5761-5826,5826-5891,7446-7514,"
			."7518-7585,8295-8364,9991-10058,10405-10469,12138-12206,12207-12265,12266-12336,"
			."14674-14742,15888-15953,15956-16023"],
		21, ["Control Region ouside of HVRI/HVRII","16366-16569,1-72,341-576"],
		22, ["Noncoding regrion ouside of Control Region","3305-3306,4401-4401,5580-5586,5656-5656,"
			."5730-5760,5892-5903,7515-7517,8270-8294,8365-8365,14743-14746,15954-15955"]
		);
	my $range_include_ref = \%range_include;

	print "=============================================\n";
	foreach my $key (sort {$a<=>$b} keys %$range_include_ref) {
		print "\t".join("\t",$key.":",$range_include_ref->{$key}->[1]."(".$range_include_ref->{$key}->[0].")")."\n";
	}
	print "=============================================\n"
		."select the range you want to include in output(default 1):";
	my $range_include_option = <STDIN>;
	chomp($range_include_option);
	$range_include_option = ($range_include_option eq "") ? 1:$range_include_option;
	if(!exists($range_include_ref->{$range_include_option})) {
		die("illegal input option\n");
	}
	my $range_include_string = $range_include_ref->{$range_include_option}->[1];
	return $range_include_string;
}

sub get_exclude_range() {
	my %range_exclude = (
		1, ["poly-C region","303-315"],
		2, ["none",""]
		);
	my $range_exclude_ref = \%range_exclude;
	print "=============================================\n";
	foreach my $key (sort {$a<=>$b} keys %$range_exclude_ref) {
		print "\t".join("\t",$key.":",$range_exclude_ref->{$key}->[1]."(".$range_exclude_ref->{$key}->[0].")")."\n";
	}
	print "=============================================\n"
		."select a range you want to exclude in output(default 1):";
	my $range_exclude_option = <STDIN>;
	chomp($range_exclude_option);
	$range_exclude_option = ($range_exclude_option eq "") ? 1:$range_exclude_option;
	if(!exists($range_exclude_ref->{$range_exclude_option})) {
		die("illegal input option\n");
	}
	my $range_exclude_string = $range_exclude_ref->{$range_exclude_option}->[1];
	return $range_exclude_string;
}

sub get_ignore_chrs() {
	my %ignore_chr_arr = (
		1 => "N,R,Y,W",
		2 => "N",
		3 => ""
		);
	print "=============================================\n";
	foreach my $key (sort {$a<=>$b} keys %ignore_chr_arr) {
		print "\t".join("\t",$key.":",$ignore_chr_arr{$key})."\n";
	}
	print "=============================================\n"
		."Select the ambiguous nucleotide characters you want to exclude(default 1):";
	my $ignore_option = <STDIN>;
	chomp ($ignore_option);
	$ignore_option = ($ignore_option eq "") ? "1":$ignore_option;
	if(!$ignore_chr_arr{$ignore_option}) {
		die("illegal input option\n");
	}
	return $ignore_chr_arr{$ignore_option};
}

sub get_diff_chrs() {
	my %diff_to_select = (
		1 => "0,1",
		2 => "0"
	);
	print "=============================================\n";
	foreach my $k (sort {$a<=>$b} keys %diff_to_select) {
		print "\t".$k.":\t".$diff_to_select{$k}."\n";
	}
	print "=============================================\n"
		."Select the definition of identical sequences\n"
		."Number of Pairwise Differences(Default option 1):";
	my $diff_select_option = <STDIN>;
	chomp ($diff_select_option);
	$diff_select_option = ($diff_select_option eq "") ? "1":$diff_select_option;
	if(!exists($diff_to_select{$diff_select_option})) {
		die("illegal input option\n");
	}
	return $diff_to_select{$diff_select_option};
}