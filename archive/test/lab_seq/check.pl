#!/usr/bin/perl -w
use strict;
use Bio::Seq;
use Bio::SeqIO;
my $dir = "E:/Scripts/test/lab_seq/";
my $in = Bio::SeqIO->new(-file => $dir.'rCRS.gb', -format => 'genbank');
my @file_name_arr = ("01_part2", "03", "04", "06");

my $rCRS = $in->next_seq();
open ( REV, ">".$dir."recovered.txt" ) or die("Cannot Write Result File");

foreach my $file_name (@file_name_arr)
{
  open( XLS, $dir.$file_name.".txt" ) or die("Cannot open XLS-file\n");
  print $file_name."\n";

  my $buffer = <XLS>;
  chomp($buffer);
  my @loc_arr = split(/\t/, $buffer);
  $buffer = <XLS>;
  chomp($buffer);
  my @ele_arr = split(/\t/, $buffer);
  my $i = 3;
  while($loc_arr[$i])
  {
    my $rCRS_res = $rCRS->subseq($loc_arr[$i],$loc_arr[$i]);
    my $provided_res = $ele_arr[$i];
    if($rCRS_res ne $provided_res)
    {
      print $loc_arr[$i]."[".$rCRS_res."/".$provided_res."]\t";
    }
    $i++;
  }
  print "\n";

}