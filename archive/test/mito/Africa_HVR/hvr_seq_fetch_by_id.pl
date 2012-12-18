#!/usr/bin/perl -w
use LWP::Simple;
my $eutils = "http://www.ncbi.nlm.nih.gov/entrez/eutils/";
my $path = "E:/Scripts/test/Africa_HVR/";

open( ID, $path."idlist.txt" ) or die("couldn't read idlist");

my $count = 1;
my $records_per_round = 5;
while( !eof(ID) )
{
  my @id_arr = ();
  for(my $i=0; $i<$records_per_round; $i++ )
  {
    my $line = <ID>;
    chomp($line);
    #my @ele_arr = split("\t",$line);
    push @id_arr, $line;
  }
  &fetch_several_seq($count,$records_per_round,@id_arr);
  $count += $records_per_round;
}

sub fetch_several_seq
{
  my ($cnt,$records_per_round,@id_arr) = @_;
  my $efetch = "efetch.fcgi?db=nuccore"
  	."&retmode=text&rettype=gb&id=".join(",",@id_arr);
  my $file_name = sprintf( "%04.0f", $cnt );
  print "fetching ".$cnt."\t-\t".($cnt+$records_per_round-1)
  	."\t(id:".join(",",@id_arr).")...\t";
  my $efetch_result = get($eutils . $efetch);
  #print $efetch_result;
  my @result_arr = split("//\n\n",$efetch_result);
  for(my $j=$cnt; $j<$cnt+$records_per_round; $j++)
  {
    open( SEQ, ">".$path."seq/".sprintf("%04.0f", $j).".txt")
    	or die("Could not create seq-file\n");
    print SEQ $result_arr[($j-1)%5]."//\n";
  }
  print "complete.\n";
}