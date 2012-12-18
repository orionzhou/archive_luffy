#!/usr/bin/perl -w
require 'chr_band.pl';
use strict;
use DBI;
my $dbh = DBI->connect( "DBI:mysql:database=seq;host=localhost,
	"genius","prodigy", {'RaiseError' => 1});
my $file_prefix = "E:/Genius/papers/CNV/";
#open( CNV, ">".$file_prefix."CNV.txt" );

my @chr_arr = (1..22,"X","Y");
my $total = 0;
foreach my $chr (@chr_arr)
{
  my $query  = "select distinct(LocusID),LocusStart,LocusEnd from DGV "
  	."where LocusChr='chr".$chr."' and VariationType="
  	."'CopyNumber'";
  my $sqr = $dbh->prepare($query);
  $sqr->execute();
  my @locus_arr;
  while(my $ref = $sqr->fetchrow_hashref())
  {
    push(@locus_arr,join("*",$ref->{'LocusID'},$ref->{'LocusStart'},$ref->{'LocusEnd'}));
  }
  my $count = scalar(@locus_arr);
  my $interval;
  my $max;
  if($count>20)
  {
    $interval = ($count-1)/(20-1);
    $max = 20;
  }
  else
  {
    $interval = 1;
    $max = $count;
  }
  for(my $i=0;$i<$max;$i++)
  {
    my $pos = int(0 + ($i*$interval));
    my @ele_arr = split("\\*",$locus_arr[$pos]);
    my $band = ideogram_localize($chr,$ele_arr[2]);
    print join("\t",$chr,@ele_arr,$band),"\n";
  }
}