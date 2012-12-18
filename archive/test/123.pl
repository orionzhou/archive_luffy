#!/usr/bin/perl -w
use strict;
#print "Please enter your name:\n\n";
#$ggg = <STDIN>;
#print ("Your name is ",$ggg,"\n");
my @temp = ("kkk",4,"bjjj",64.34);
$" = "===";
my $array = join($",@temp);
print $array,"\n";
system "echo @temp";