#!/usr/bin/perl -w
use LWP::Simple;
my $eutils = "http://www.ncbi.nlm.nih.gov/entrez/eutils/";

&fetch_mito_seqs();

sub fetch_mito_seqs
{
  my $esearch = "esearch.fcgi?db=nuccore"
	."&term=Homo[Organism]+AND+mitochondrion[All+Fields]+AND+15000:17000[SLEN]"
        ."+NOT+pseudogene[All+Fields]&usehistory=y";
  my $count_query = get($eutils . $esearch ."&retmax=10");
  $count_query =~
    m|<Count>(\d+)</Count>.*<QueryKey>(\d+)</QueryKey>.*<WebEnv>(\S+)</WebEnv>|s;
  my $total_count = $1;
  my $querykey = $2;
  my $webenv = $3;
  print "Totally ".$total_count." Mito_Seqs\n\n";
  my $buffer = "";
  for(my $i=0; $i<$total_count; $i++)
  {
    open( Seq, "E:/Scripts/test/mito/mito_seq/".sprintf("%04.0f",$i+1).".txt" )
    	or die("could not open file");
    $buffer = <Seq>;
    if(eof)
    {
      printf("%04.0f is missing\n",$i+1);
      &fetch_one_seq($i,$querykey,$webenv);
    }
  }
}

sub fetch_one_seq
{
  my ($id,$querykey,$webenv) = @_;
  my $efetch = "efetch.fcgi?db=nuccore"
  	."&retmode=text&rettype=gb&query_key=$querykey&WebEnv=$webenv"
        ."&retstart=".$id."&retmax=1";
  my $file_name = sprintf("%04.0f",$id+1);
  print "\tfetching ".$file_name."...\t\t";
  my $efetch_result = get($eutils . $efetch);
  open( MITO, ">E:/Scripts/test/mito/mito_seq/".$file_name.".txt" )
  	or die("Cannot create file");
  print MITO $efetch_result;
  print "complete.\n";
}