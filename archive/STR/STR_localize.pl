#!/usr/bin/perl -w
use strict;
use DBI;
my $dbh = DBI->connect( "DBI:mysql:database=seq;host=202.117.161.53",
	"genius","prodigy", {'RaiseError' => 1});
my $file_pre = "E:/Scripts/STR/";
open( STR, $file_pre."STRs_2.txt" );
my @ele_arr;
my $count_localized = 0;
my $count_undetermined = 0;
while(<STR>)
{
  chomp;
  @ele_arr = split("\t",$_);
  if($ele_arr[9] || !$ele_arr[6])
  {
    print "\n";
    $count_localized ++;
  }
  else
  {
    my $string_query = "select * from ideogram where chr='".$ele_arr[4]
    	."' and bp_start<".$ele_arr[6]." and bp_stop>".$ele_arr[6]
        ." and band regexp '^[0-9]{2}\$'" ;
    my $sqr = $dbh->prepare($string_query);
    $sqr->execute();
    while(my $ref = $sqr->fetchrow_hashref())
    {
      print join("\t",$ref->{'chr'}.$ref->{'arm'}.$ref->{'band'},
      	$ref->{'bp_start'},$ref->{'bp_stop'});
    }
    print "\n";
    $count_undetermined ++;
  }
}
$dbh->disconnect();