#!/usr/bin/perl -w
use LWP::Simple;
my $eutils = "http://www.ncbi.nlm.nih.gov/entrez/eutils/";
my $path = "E:/Scripts/test/mito/";

my $path_dir = $path."idlist/";
my $path_seq = $path."mito_seq/";
opendir( IDLIST, $path_dir ) or die("couldn't read idlist\n");
while( my $list = readdir(IDLIST) ) {
     if( -f $path_dir.$list ) {
     	print $list."\n";
        open( ONELIST, $path_dir.$list ) or die("could not read $list\n");
        my $count = 1;
	my $records_per_round = 5;
	while( !eof(ONELIST) )
	{
	  my @id_arr = ();
	  my @acc_arr = ();
          print "\t";
	  for(my $i=0; $i<$records_per_round; $i++ )
	  {
	    my $line = <ONELIST>;
            if( $line ) {
	    chomp($line);
            my @ele_arr = split("\t",$line);
	    push @id_arr, $ele_arr[1];
	    push @acc_arr, $ele_arr[0];
            print $ele_arr[0]." ";
            }
	  }
	  &fetch_several_seq(substr($list,0,2),$count,join(",",@id_arr),@acc_arr);
	  $count += $records_per_round;
	}
     }
}

sub fetch_several_seq
{
  my ($listid,$cnt,$ids,@acc_arr) = @_;
  my $efetch = "efetch.fcgi?db=nuccore"
  	."&retmode=text&rettype=gb&id=".$ids;
  my $efetch_result = get($eutils . $efetch);
  #print $efetch_result;
  my @result_arr = split("//\n\n",$efetch_result);
  for(my $j=$cnt; $j<$cnt+scalar(@acc_arr); $j++)
  {
    open( SEQ, ">".$path_seq.substr($listid,0,2)."_".sprintf("%02.0f", $j)
    	."_".$acc_arr[($j-1)%5].".gb") or die("Could not create seq-file\n");
    print SEQ $result_arr[($j-1)%5]."//\n";
  }
  print "complete\n";
}
