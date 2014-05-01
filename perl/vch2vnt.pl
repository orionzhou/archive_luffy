#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  vcf2vnt.pl - convert a VCF (human readable) file to SNP/IDM file

=head1 SYNOPSIS
  
  vcf2snp.pl [-help] [-in input] [-out output dir] 

  Options:
      -h (--help)       brief help message
      -i (--in)         input
      -o (--out)        output dir

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib $FindBin::Bin;
use Data::Dumper;
use File::Path qw/make_path remove_tree/;
use List::Util qw/min max sum/;

my ($fi, $do) = ('') x 2;
my $fhi;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "in|i=s"  => \$fi,
  "out|o=s" => \$do,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$do;

$fhi = \*STDIN;
unless ($fi eq "stdin" || $fi eq "-" || $fi eq "") {
  open ($fhi, "<$fi") || die "Can't open file $fi: $!\n";
}
-d $do || make_path($do);

my $hf = {
  1  => "$do/snp",
  0  => "$do/het", 
  9  => "$do/idm",
};
for my $type (keys(%$hf)) {
  my $f = $hf->{$type};
  open(my $fh, ">$f") || die "cannot write file $f\n";
  $hf->{$type} = $fh;
}

while( <$fhi> ) {
  chomp;
  my $line = $_;
  my @ps = split("\t", $line);
  @ps == 5 || die "not 5 lines:\n$line\n";
  my ($chr, $beg, $ref, $alts, $str) = @ps;
  my @alts = split(",", $alts);
  
  my @types;
  for my $alt (@alts) {
    my $type = (length($ref) == 1 && length($alt) == 1) ? 1 : 9;
    push @types, $type;
  }
  if($str =~ /([A-Za-z0-9\-_]+)\=([\d\.])\/([\d\.])/) {
    if($2 eq '.' || $3 eq '.') {
    } elsif($2 eq $3) {
      my ($alt, $type) = ($alts[$2-1], $types[$2-1]);
      my $fho = $hf->{$type};
      
      if($type == 1) {
        print $fho join("\t", $chr, $beg, $ref, $alt)."\n";
      } else {
        my $end = $beg + length($ref);
        ($ref, $alt) = (substr($ref, 1), substr($alt, 1));
        print $fho join("\t", $chr, $beg, $end, $ref, $alt)."\n";
      }
    } elsif($2 ne $3 && $2 == 0) {
      my ($alt, $type) = ($alts[$3-1], $types[$3-1]);
      $type = 0;
      my $fho = $hf->{$type};
      print $fho join("\t", $chr, $beg, $ref, $alt)."\n";
    }
  } else {
      die "unknown allele: $str\n";
  }
}
close $fhi;
for (keys(%$hf)) { close $hf->{$_}; }


