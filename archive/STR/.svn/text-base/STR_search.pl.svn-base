#!/usr/bin/perl -w
use DBI;
my $dbh = DBI->connect( "DBI:mysql:database=seq;host=202.117.161.53",
	"genius","prodigy", {'RaiseError' => 1});


my $file_prefix = "E:/Scripts/STR/";
open( STRin, $file_prefix."STRs.txt" ) or die("Couldn't open file!");
#open( STRout, ">".$file_prefix."STR_id.txt" ) or die("Couldn't create file!");
my @unfound_arr;
while( <STRin> )
{
    chomp;
    my @ele_arr = split("\t",$_);
    #print join(" ",@ele_arr[0,1,3]),"\n";
    my $sqr = $dbh->prepare("SELECT * from sts where feature_name='".$ele_arr[0]
    	."' and group_label='reference'");
    $sqr->execute();
    my $count = 0;
    while(my $ref = $sqr->fetchrow_hashref())
    {
      print join("\t", $ref->{'feature_id'},$ref->{'chr_start'},$ref->{'chr_stop'},
      	),"\t";
      $count ++;
    }
    if($count == 0)
    {
      print "none";
      push(@unfound_arr,$ele_arr[0]);
    }
    print "\n";
    #print SNPs "\n\n";
}
print "\n\n",scalar(@unfound_arr)," records unfound:\n";
print join("\n",@unfound_arr);
