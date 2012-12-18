#!/usr/bin/perl -w
use strict;
open( GENE, "F:/Genius/STR/perl+php+fasta/result/Y_gene.txt" ) or die("What's wrong");
open ( STR, "F:/Genius/STR/perl+php+fasta/result/result.txt" ) or die("What's wrong");
open ( GENE_STR, ">F:/Genius/STR/perl+php+fasta/result/gene_str.txt" ) or die("What's wrong");
my @str_arr;
my @gene_arr;
my $i=0;
while (<STR>)
{
	chomp;
        push @str_arr, $_;
}
my $j=0,$i=0;
while (<GENE>)
{
	chomp;
        my @gene_info_arr = split(/\t/,$_);
        if( $gene_info_arr[0] =~ /^\d+$/ )
        {
        	push @gene_arr, $_;
        }
}
close(STR);
close(GENE);
while( $str_arr[$j] && $str_arr[$j+1] )
{
	print $str_arr[$j],"\t";
        print GENE_STR $str_arr[$j],"\t";
        my @this_str = split(/\t/,$str_arr[$j]);
        my @next_str = split(/\t/,$str_arr[$j+1]);
        my @gene = split(/\t/,$gene_arr[$i]);
        while( $gene[0] > $this_str[8] &&
         ($gene[0] < $next_str[8] || $gene[0] < $next_str[7]) )
        {
        	print $gene[0]."-".$gene[2],"\t";
                print GENE_STR $gene[0]."-".$gene[2],"\t";
                $i++;
                @gene = split(/\t/,$gene_arr[$i]);
        }
        print "\n";
        print GENE_STR "\n";
        $j++;
}
print "\n",$i,"genes has been located.";