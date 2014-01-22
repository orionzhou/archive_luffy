#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;
use Net::FTP;

my $ftp = Net::FTP->new("data.ncgr.org", Debug=>0) || die "failed connect\n";
$ftp->login("pipestone", "smRRkn++") || die "login failed ", $ftp->message;
$ftp->cwd("old/variants/20131223") || die $ftp->message;

my $dir = "/home/youngn/zhoup/Data/misc3/hapmap_mt40/30_vnt";
for my $fn ($ftp->ls()) {
    $ftp->get($fn, "$dir/$fn") || die $ftp->message;
}
