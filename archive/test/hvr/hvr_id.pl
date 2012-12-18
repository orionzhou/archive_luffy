#!/usr/bin/perl -w
use LWP::Simple;
my $eutils = "http://www.ncbi.nlm.nih.gov/entrez/eutils/";

&fetch_hvr_id();

sub fetch_mito_id
{
  my $esearch = "esearch.fcgi?db=Nucleotide"
	."&term=hypervariable+segment+AND+Chinese+AND+mitochondrial";
  my $count_query = get($eutils . $esearch ."&retmax=10000");
  $count_query =~
    m|<Count>(\d+)</Count>.*<QueryKey>(\d+)</QueryKey>.*<WebEnv>(\S+)</WebEnv>|s;
  my $total_count = $1;
  my $querykey = $2;
  my $webenv = $3;
  print "Totally ".$total_count." hvr_Seqs\n\n";
  #print $count_query;
  open( ID, ">E:/Scripts/test/hvr/idlist.txt" ) or die("Cannot create file");
  while( $count_query =~ m|<Id>(\d+)</Id>|g )
  {
    print ID $1,"\n";
  }
}