#!/usr/bin/perl
open( FILE, "./blast_1.txt") or die("Could not open file!");
$i=0;
while(<FILE>)
{
	chomp;
	@ele_arr = split(/\t/,$_);
	print ($i+1,"---",$ele_arr[0],":",scalar(@ele_arr)," elements---");
	if(scalar(@ele_arr)==1)
	{print ("this line is discarded.\n");}
	else
	{
		my $fname = sprintf("%03d",$i+1) . "_" . $ele_arr[0] . ".txt" ;
		print STDOUT (">>>> ",$fname,"\n");
		open ( WR , ">" . $fname) or die ("could not write file.");
		print WR (">",$ele_arr[0],"\n");
		print WR ($ele_arr[1]," ",$ele_arr[2]);
		close(WR);
	}
	$i++;
}
close(FILE);
