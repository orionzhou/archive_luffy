#!/usr/bin/perl -w
use strict;
my $file_pre = "E:/Scripts/STR/";
open ( IDEO, $file_pre."ideogram.txt" );
my @ele_arr;
while(<IDEO>)
{
  chomp;
  @ele_arr = split("\t",$_);
  if($ele_arr[2])
  {
    if($ele_arr[2] =~ /^(\d)+$/ || $ele_arr[2] eq 'ter')
    {
      print join ("\t",@ele_arr[0..2,5,6,9]),"\n";
    }
  }
  else
  {
    print join ("\t",@ele_arr[0..2,5,6,9]),"\n";
  }
}