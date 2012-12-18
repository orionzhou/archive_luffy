#!/usr/bin/perl -w
use LWP::Simple;
my $utils = "http://www.ncbi.nlm.nih.gov/entrez/eutils";
my $db     = "unists";
my $report = "native";

my $file_prefix = "E:/Scripts/STR/";
open( STRin, $file_prefix."STRs.txt" ) or die("Couldn't open file!");
open( STRout, ">".$file_prefix."STR_id.txt" ) or die("Couldn't create file!");
my @unfound_arr;
while( <STRin> )
{
    chomp;
    my $row = $_;
    $row =~ s/\s//i ;
    my @ele_arr = split("\t",$row);
    my $query  = $ele_arr[0] . " AND human[ORGN]" ;
    #print join(" ",@ele_arr[0,1,3]),"\n";
    my $esearch = "$utils/esearch.fcgi?" . "db=$db&retmax=10&usehistory=y&term=";
    my $esearch_result = get($esearch . $query);
    $esearch_result =~ m|<Count>(\d+)</Count>.*<IdList>(.+)</IdList>|s;
    my $count = $1;
    my $ids = $2;
    if($count>0 && $count<8)
    {
      #print "$count\n";
      #print SNPs "$count\n";
      while( $ids =~ m|<Id>(\d+)</Id>|sg )
      {
        print STRout $1,"\t";
        print $1,"\t";
        #print SNPs $1,"\t";
      }
      #&esearch_result_handler($esearch_result);
      #print $esearch_result;
    }
    else
    {
      print STRout "none";
      push (@unfound_arr,$query);
      #print SNPs "0";
    }
    print STRout "\n";
    print "\n";
    #print SNPs "\n\n";
}
print "\n\n",scalar(@unfound_arr)," records unfound:\n";
print join("\n",@unfound_arr);
