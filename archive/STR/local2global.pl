#!/usr/bin/perl -w
use strict;
open ( PRE, "F:/Genius/STR/perl+php+fasta/result/result1.txt" ) or die("could not open file");
open ( RESULT, ">>F:/Genius/STR/perl+php+fasta/result/result.txt" ) or die("could not open file");
my $j=1;
while( <PRE> )
{
	chomp;
        my @ele_arr = split( /\t/, $_ );
        print join("--",@ele_arr),"\n";
	my ($contig) = $ele_arr[0] =~ /(NT_[0-9]+\.[0-9]{1,2})/;
	my ($part) = $ele_arr[0] =~ /part([0-9]{4})/;
	my $start_contig = ($part-1) * 400 * 70 + 1 ;
	my $start_level2 = $ele_arr[5] + $start_contig - 1 ;
        my $end_level2 = $ele_arr[6] + $start_contig - 1 ;
	my %loc_contig = ( 'NT_113967.1'=>1, 'NT_113968.1'=>84822, 'NT_113969.1'=>201385,
	'NT_113970.1'=>1017558, 'NT_113971.1'=>1104114, 'NT_113972.1'=>1274235,
	'NT_113973.1'=>2128239, 'NT_011896.9'=>2709521, 'NT_086998.1'=>9024956,
	'NT_011878.9'=>9901323, 'NT_087001.1'=>11214554, 'NT_113819.1'=>11653955,
	'NT_011875.11'=>12308579, 'NT_011903.12'=>22360817, 'NT_025975.2'=>57228750,
	'NT_091573.1'=>57377045, 'NT_113974.1'=>57443438);
	my $start_level3 = $start_level2 + $loc_contig{$contig} - 1 ;
        my $end_level3 = $end_level2 + $loc_contig{$contig} - 1 ;
        print RESULT join("\t",@ele_arr,$start_level3,$end_level3), "\n";
        $j++;
}
close( PRE );
close( RESULT );
