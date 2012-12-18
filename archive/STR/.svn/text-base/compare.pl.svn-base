#!/usr/bin/perl -w
use strict;
my $found_path = "F:/Genius/STR/perl+php+fasta/result/result.txt" ;
my $reported_path = "F:/Genius/STR/perl+php+fasta/reported_STRs/result.txt" ;
my $comparison_path = "F:/Genius/STR/perl+php+fasta/result/comparison.txt" ;
open( FOUND, $found_path ) or die("could not open finding file");
open( REPORTED, $reported_path ) or die("could not open reported file");
open( COMPARISON, ">".$comparison_path ) or die("could not write into comparison file");
my @reported_arr;
while(<REPORTED>)
{
	chomp;
        my @ele_reported_arr = split(/\t/,$_);
        if(scalar(@ele_reported_arr)>3)
        {
        	push @reported_arr, $_ ;
                #print $ele_reported_arr[0]," will be considered.\n";
        }
        else
        {
        	#print $ele_reported_arr[0]," is discared.\n";
        }
}
@reported_arr = sort(@reported_arr);
my $hit = 0;
#print "A total of ",scalar(@reported_arr)," STR was identified within the 416 reported.\n";
while(<FOUND>)
{
	chomp;
        my @ele_found_arr = split(/\t/,$_);
        my $start_found = $ele_found_arr[7];
        my $end_found = $ele_found_arr[8];
        my $motif_found = $ele_found_arr[3];
        my $repeat_found = $ele_found_arr[4];
        print join("\t",$motif_found,$repeat_found,$start_found,$end_found),"\t";
        print COMPARISON join("\t",@ele_found_arr),"\t";
        foreach my $reported (@reported_arr)
        {
        	my @ele_reported = split(/\t/,$reported);
                my $i=3;
                while($ele_reported[$i])
                {
                	if( abs($ele_reported[$i+2]-$start_found)<=10
                         && abs($ele_reported[$i+3]-$end_found)<=10 )
                        {
                        	print join("\t",$ele_reported[0],$ele_reported[1],
                                $ele_reported[$i],$ele_reported[$i+1],$ele_reported[$i+2],
                                $ele_reported[$i+3]),"\t";
                                print COMPARISON join("\t",$ele_reported[0],$ele_reported[1],
                                $ele_reported[$i],$ele_reported[$i+1],$ele_reported[$i+2],
                                $ele_reported[$i+3]),"\t";
                                $hit ++ ;
                        }
                        $i +=4 ;
                }
        }
        print "\n";
        print COMPARISON "\n";
        #print join("==",@ele_found_arr),"\n";
}
print "\nA total of $hit hits!";
