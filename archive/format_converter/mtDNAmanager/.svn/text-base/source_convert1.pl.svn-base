#!/usr/bin/perl -w
use strict;

my $dir = "/home/orion/Scripts/test/align/";
my $file = "source_hongkong";

my @range = ("16024-16365","40-302","316-370");
open(IN, $dir.$file) or die("Cannnot open source file");

my $out = "tmp";
open(STDOUT, ">".$dir.$out) or die("Cannot write result file");
while(<IN>) {
	chomp;
	if($_) {
		my @ele_arr = split("\t",$_);
		my $number = shift(@ele_arr);
		my $group = shift(@ele_arr);
		for(my $i=1;$i<=$number;$i++) {
			my %mut_arr;
			foreach my $ele (@ele_arr) {
				if($ele =~ /^((\d+)[A-Z])$/) {
					$mut_arr{$2} = $1;
				} elsif($ele=~/^((\d+\.\d+)[A-Z])$/) {
					$mut_arr{$2} = $1;
				} elsif($ele=~/^\-(\d+)$/) {
					$mut_arr{$1} = $1."d";
				} else {
					die("Unrecognized variant");
				}
			}
			my $name = ($i==1)?$group:$group.".".$i;
			print $name."\t";
			foreach (sort {$a <=> $b } keys(%mut_arr)) {
				print &printw($_,$mut_arr{$_},@range);
			}
			print "\n";
		}
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