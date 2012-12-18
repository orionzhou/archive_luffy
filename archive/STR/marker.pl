#!/usr/bin/perl -w
use strict;
my $file_prefix = "D:/seq/marker_info/";
my @chr_arr = (1..23);
my @ele_arr;
my $temp;
#open(WR,">".$file_prefix."marker.txt");
foreach my $chr (@chr_arr)
{
  open(MAR,$file_prefix."info".$chr.".txt");
  $temp = <MAR>;
  $temp = <MAR>;
  while(<MAR>)
  {
    chomp;
    @ele_arr = split(" ",$_);
    if($chr==23)
    {
      print "X\t";
    }
    else
    {
      print $chr."\t";
    }
    if(scalar(@ele_arr) == 11)
    {
      if($ele_arr[0] eq "*")
      {
        print "*";
      }
      elsif($ele_arr[0] eq "x")
      {
        print "x";
      }
    }
    print "\t";
    print join("\t",@ele_arr[scalar(@ele_arr)-10..scalar(@ele_arr)-1]);
    print "\n";
  }
}