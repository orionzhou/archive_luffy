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
		my @ary_tmp = split(/\W+/,$_);
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

my %range = (
		1, ["mtGenome","1-16569"], 
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
my $range_ref = \%range;

my %pos_snp_arr;
my @variant_arr;
my $count_snp_unique=0;
my $count_sites_unique=0;
my %count_variant;

&write_head();
for(my $i=0; $i<@snp_id_arr; $i++) {
	$query = "SELECT snp FROM mt_snp WHERE id=".$snp_id_arr[$i];
	$sth = $dbh->prepare($query);
	$sth->execute();
	my $rs = $sth->fetchrow_array();

	print OUT sprintf("%03.0f",$i+1)."_".$snp_id_arr[$i]."\t";
	foreach my $key (sort {$a<=>$b} keys %$range_ref) {
		our $range_include_string = $range_ref->{$key}->[1];
		print OUT &snp_summary(&reconstruct_snp_string($rs))."\t";
	}
	print OUT "\n";
}

sub write_head() {
	print OUT "\t";
	foreach my $key (sort {$a<=>$b} keys %$range_ref) {
		my $range_string = $range_ref->{$key}->[1];
		my @range_pair_arr = split(",",$range_string);
		my $length = 0;
		foreach my $range_pair (@range_pair_arr) {
			my @tmp = split("-",$range_pair);
			if(@tmp == 1) {
				$length += 1;
			} elsif(@tmp == 2) {
				$length += ($tmp[1]-$tmp[0]+1);
			}
		}
		print OUT $range_ref->{$key}->[0]."\t\t\t\t".$length."\t";
	}
	print OUT "\n";
	print OUT "\t".join("\t","transition","transversion","deletion","insertion","total")."\n";
}


sub snp_summary() {
	my ($snp_str) = @_;
	my @snp_arr = split(" ",$snp_str);
	my $buffer = "";
	
	my ($count_transition, $count_transversion, $count_deletion, $count_insertion) = (0, 0, 0, 0);
	foreach my $snp (@snp_arr) {
		if($snp =~ /^([.0-9]+)([A-Zd])$/) {
			my $pos = $1;
			my $value = $2;
			if($pos =~ /\./) {
				$count_insertion ++;
			} elsif($value =~ /d/) {
				$count_deletion ++;
			} else {
				#my $seqobj = Bio::Seq->new(-id=>"Ref_nt_$pos", -seq=>$rCRS->subseq($pos,$pos));
				#my $nt_ref_revcom = $seqobj->revcom()->seq();
				#$buffer .= $seqobj->seq().$pos.$value." ";
				my $nt = $rCRS->subseq($pos,$pos);
				my $nt_transitioned = &base_transition($nt);
				if($nt_transitioned eq $value) {
					$count_transition ++;
				} elsif($nt_transitioned=~/[ATCG]/i && $value=~/[ATCG]/i) {
					$count_transversion ++;
				} else {
					print "Unknown mutation character at $pos-$value\n";
					exit;
				}
			}
		} elsif($snp ne "") {
			$buffer = "error snp string value\n";
			exit;
		}
	}
	$buffer .= join("\t", $count_transition, $count_transversion, $count_deletion,
		$count_insertion, $count_transition+$count_transversion+$count_deletion+$count_insertion);
	return $buffer; 
}