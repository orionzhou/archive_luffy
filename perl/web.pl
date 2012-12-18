#!/usr/bin/perl
#use strict; use Init; use Common; use Localdb; use Run; 
use Net::SFTP; use Net::SSH::Perl;
use Time::HiRes qw/gettimeofday tv_interval/; use Data::Dumper;
use List::Util qw/min max sum/; use POSIX qw/ceil floor/; use Math::Round qw/round/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;

my @hosts = qw/l4 l5 l6 l13 l14 l15 l16/;
my ($user, $pw) = qw/zhoup Abc123!!/;
for (@hosts) {
    my $host = "$_.msi.umn.edu";
    my $ssh = Net::SSH::Perl->new($host);
    $ssh->login($user, $pw);
    my $cmd = "top -b -n 1";
    my ($stdout, $stderr, $exit) = $ssh->cmd($cmd);
    my @ary = split("\n", $stdout);
    print $host."\n";
    print "\t".join("\t\n", @ary[1..4])."\n";
}  
