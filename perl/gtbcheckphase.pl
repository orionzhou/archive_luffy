#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gtbcheckphase.pl - check and fix a Gtb file

=head1 SYNOPSIS
  
  gtbcheckphase.pl [-help] [-in input-file] [-seq refseq-file] [-out output-file]

  Options:
      -help   brief help message
      -in     input file
      -out    output file
      -seq    refseq file

=cut
  
#### END of POD documentation.
#-----------------------------------------------------------------------------
use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;
use Gtb;

my ($fi, $fo, $fs) = ('') x 3;
my $help_flag;

#----------------------------------- MAIN -----------------------------------#
GetOptions(
    "help|h"  => \$help_flag,
    "in|i=s"  => \$fi,
    "out|o=s" => \$fo,
    "seq|s=s" => \$fs,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo || !$fs;

my ($fhi, $fho);
if ($fi eq "stdin" || $fi eq "-") {
    $fhi = \*STDIN;
} else {
    open ($fhi, $fi) || die "Can't open file $fi: $!\n";
}

if ($fo eq "stdout" || $fo eq "-") {
    $fho = \*STDOUT;
} else {
    open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

print $fho join("\t", qw/id parent chr beg end strand locE locI locC loc5 loc3 phase source conf cat1 cat2 cat3 note/)."\n";
my $n_fixed = 0;
while( <$fhi> ) {
    chomp;
    next if /(^id)|(^\#)|(^\s*$)/;
    my $ps = [ split "\t" ];
    next unless @$ps >= 18;
    my @phases = gtbcheckphase($ps, $fs);
    if(@phases) {
        $n_fixed ++;
        $ps->[11] = join(",", @phases);
    }
    print $fho join("\t", @$ps)."\n";
}
print "  $n_fixed non-0 frames fixed\n";
close $fhi;
close $fho;

__END__
