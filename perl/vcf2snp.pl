#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  vcf2snp.pl - convert a VCF (human readable) file to SNP file

=head1 SYNOPSIS
  
  vcf2snp.pl [-help] [-in input] [-out output] 

  Options:
      -h, --help       brief help message
      -i, --in         input
      -o, --out        output

=cut
  
#### END of POD documentation.
#-----------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib $FindBin::Bin;
use Data::Dumper;
use List::Util qw/min max sum/;

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

while( <$fhi> ) {
    chomp;
    my $line = $_;
    my @ps = split("\t", $line);
    @ps == 5 || die "not 5 lines:\n$line\n";
    my ($chr, $beg, $ref, $alts, $str) = @ps;
    my $end = $beg + length($ref) - 1;
    my @alts = split(",", $alts);
    
    my @types;
    for my $alt (@alts) {
        my $type;
        if(length($ref) == 1 && length($alt) == 1) {
            $type = 'snp';
        } elsif(length($ref) > length($alt)) {
            $type = ($ref =~ /^$alt/i) ? 'del' : 'mnp';
        } elsif(length($ref) < length($alt)) {
            $type = ($alt =~ /^$ref/i) ? "ins" : "mnp"; 
        } else {
            print $line."\n";
            $type = 'mnp';
        }
        push @types, $type;
    }
    if($str =~ /([A-Za-z0-9\-_]+)\=([\d\.])\/([\d\.])/) {
        if($2 eq '.' || $3 eq '.') {
        } elsif($2 eq $3) {
            print $fho join("\t", $chr, $beg, $end, $ref, $alts[$2-1], $types[$2-1])."\n";
        } elsif($2 ne $3 && $2 == 0) {
            print $fho join("\t", $chr, $beg, $end, $ref, '', 'het')."\n";
        }
    } else {
        die "unknown allele: $str\n";
    }
}
close $fhi;
close $fho;


