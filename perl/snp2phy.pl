#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  snp2phy.pl - convert a snp file to phylip file

=head1 SYNOPSIS
  
  snp2phy.pl [-help] [options] [-in input] [-out output] 

  Options:
      -h, --help       brief help message
      -i, --in         input
      -o, --out        output

=head1 VERSION
  
  0.1
  
=cut
  
#### END of POD documentation.
#-----------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib $FindBin::Bin;
use Data::Dumper;
use List::Util qw/min max sum/; use POSIX qw/ceil floor/;

my ($fi, $fo) = ('') x 2;
my $fhi;
my $fho;
my $help_flag;

#----------------------------------- MAIN -----------------------------------#
GetOptions(
    "help|h"  => \$help_flag,
    "in|i=s"  => \$fi,
    "out|o=s" => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;

$fhi = \*STDIN;
$fho = \*STDOUT;
unless ($fi eq "stdin" || $fi eq "-" || $fi eq "") {
    open ($fhi, "<$fi") || die "Can't open file $fi: $!\n";
}
unless ($fo eq "stdout" || $fo eq "-" || $fo eq "") {
    open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

my (@names, @poss, @data);
my ($nind, $npos) = (0, 0);

my $first = 1;
while( <$fhi> ) {
    chomp;
    my $line = $_;
    my @ps = split("\t", $line);
    if($first) {
        $nind = @ps - 4;
        @data = map {""} (0..$nind-1);
    } else {
        die "not $nind + 4 cols\n" unless $nind + 4 == @ps;
    }
    my ($ref, $alt) = ($ps[2], $ps[3]);
    die "ref[$ref] alt[$alt] not SNP\n" if length($ref) != 1 || length($alt) != 1;
    push @poss, $ps[1];
    $npos ++;

    for my $i (0..$nind-1) {
        if($ps[$i+4] =~ /([A-Za-z0-9\-_]+)\=([01\.])\/([01\.])/) {
            push @names, $1 if $first;
            my $nt = ($2 eq "0" && $3 eq "0") ? $ref : ($2 eq "1" && $3 eq "1") ? $alt : "N";
            $data[$i] .= $nt;
        } else {
            die "unknown allele: $ps[$i+4]\n";
        }
    }
    $first = 0 if $first;
}
close $fhi;

print $fho "$nind $npos\n";
for my $r (0..ceil($npos/250)-1) {
    for my $i (0..$nind-1) {
        my $str = substr($data[$i], $r*250, 250);
        if($r == 0) {
            printf $fho "%-20s%s\n", $names[$i], $str;
        } else {
            printf $fho "%s\n", $str;
        }
    }
    print $fho "\n";
}
close $fho;



