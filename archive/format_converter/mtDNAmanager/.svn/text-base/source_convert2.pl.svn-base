#!/usr/bin/perl -w
use strict;
use Bio::Seq;
use Bio::SeqIO;

my $dir = "/home/orion/Scripts/test/align/";
my $file = "source_japan";

my $rCRS_path = "rCRS.gb";
my $rCRS_file = Bio::SeqIO->new(-file => $rCRS_path, -format=>'genbank');
my $rCRS = $rCRS_file->next_seq();

my @range = ("16024-16365","40-302","316-370");
open(IN, $dir.$file) or die("Cannnot open source file");

my $out = "tmp";
open(STDOUT, ">".$dir.$out) or die("Cannot write result file");
while(<IN>) {
	chomp;
	if($_) {
		my @ele_arr = split("\t",$_);
		my $id = shift(@ele_arr);
		my $group = shift(@ele_arr);
		shift(@ele_arr);
		my %mut_arr;
		for(my $i=1;$i<=2;$i++) {
			my $tmp = shift(@ele_arr);
			if($tmp) {
				my @ele_arr = split(" ",$tmp);
				foreach my $ele (@ele_arr) {
					if($ele=~/N/ || $ele=~/p/) {
						#do nothing
					} elsif($ele =~ /^((\d+)[A-Z])$/) {
						$mut_arr{$2} = $1;
					} elsif($ele=~/^((\d+\.\d+)[A-Z])$/) {
						$mut_arr{$2} = $1;
					} elsif($ele=~/^(\d+)d$/) {
						$mut_arr{$1} = $1."d";
					} elsif($ele=~/^(\d+)$/) {
						my $tmpseq = Bio::Seq->new(-seq=>$rCRS->subseq($1,$1));
						$mut_arr{$1} = $1.$tmpseq->revcom()->seq();
					}else {
						die("Unrecognized variant");
					}
				}
			}
		}
		print $id.".".$group."\t";
		foreach (sort {$a <=> $b } keys(%mut_arr)) {
			print &printw($_,$mut_arr{$_},@range);
		}
		print "\n";
	}
}

sub printw() {
	my ($pos,$value,@range) = @_;
	my @begin;
    my @end;
    if(@range) {
      foreach my $ele (@range) {
        my @pair = split("-",$ele);
        push (@begin,$pair[0]);
        push (@end,$pair[1]);
      }
    } else {
      push (@begin,1);
      push (@end,16569);
    }
	my $flag = 0;
	my $i = 0;
    foreach (@begin) {
    	if($pos>=$_ && $pos<=$end[$i]) {
        	$flag = 1;
        }
        $i++;
    }
    if($flag == 1) {
    	return $value."\t";
    } else {
    	return "";
    }
}