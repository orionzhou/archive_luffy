#!/usr/bin/perl -w
use LWP::Simple;
my $eutils = "http://www.ncbi.nlm.nih.gov/entrez/eutils/";

&fetch_hvr_id();

sub fetch_hvr_id
{
  my $esearch = "esearch.fcgi?db=Nucleotide"
	."&term=mitochondrial+AND+100:1500[SLEN]+AND+African+AND+human[ORGN]";
  my $count_query = get($eutils . $esearch ."&retmax=10000");
  #print $count_query,"\n";
  $count_query =~
    m|<Count>(\d+)</Count>|s;
  my $total_count = $1;
  print "Totally ".$total_count." African_hvr_Seqs\n\n";
  #print $count_query;
  open( ID, ">E:/Scripts/test/Africa_HVR/idlist.txt" ) or die("Cannot create file");
  while( $count_query =~ m|<Id>(\d+)</Id>|g )
  {
    print ID $1,"\n";
  }
}