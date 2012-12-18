#!/usr/bin/perl -w
use LWP::Simple;
my $eutils = "http://www.ncbi.nlm.nih.gov/entrez/eutils/";
my $dir = "E:/Scripts/test/hvr/";

&id_fetch();
#&seq_fetch();
sub id_fetch
{
  for(my $i=57; $i<=76; $i++)
  {
    my $ACC = "EU2574" . $i;
    my $esearch = "esearch.fcgi?db=nuccore&term=$ACC";
    my $query1 = get( $eutils . $esearch );
    $query1 =~ m|<Count>(\d+)</Count>|s;
    my $count = $1;
    if($count == 1)
    {
      $query1 =~ m|<Id>(\d+)</Id>|s;
      my $id = $1;
      print $ACC."\t".$id."\n";
    }
    else
    {
      print "More than 1 records found!\n";
    }
  }
}

sub seq_fetch
{
  open( IDLIST, $dir."ACC-id.txt" ) or die("Cannot open ACC-id file\n");
  while( !eof(IDLIST) )
  {
    my $line = <IDLIST>;
    chomp($line);
    my @ele_arr = split("\t",$line);
    my $efetch = "efetch.fcgi?db=nuccore"
  	."&retmode=text&rettype=gb&id=".$ele_arr[1];
    print $ele_arr[0]."...\t\t";
    my $efetch_result = get($eutils . $efetch);
    open( SEQ, ">".$dir."seq/".$ele_arr[0].".txt" )
    	or die("Cannot create seq file");
    print SEQ $efetch_result;
    print "complete.\n";
  }
}