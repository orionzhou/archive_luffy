#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  galbest.pl - Filter a GAL file keeping only best hit(s)

=head1 SYNOPSIS
  
  galbest.pl [-help] [-in input-file] [-out output-file]

  Options:
      -help   brief help message
      -in     input file
      -out    output file

=cut
  
#### END of POD documentation.
#-----------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Location;
use Gal;

my ($fi, $fo) = ('') x 2;
my ($fhi, $fho);
my $help_flag;

#----------------------------------- MAIN -----------------------------------#
GetOptions(
    "help|h"   => \$help_flag,
    "in|i=s"   => \$fi,
    "out|o=s"  => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

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

print $fho join("\t", @HEAD_GAL)."\n";
my $h;
while( <$fhi> ) {
    chomp;
    next if /(^id)|(^\#)|(^\s*$)/;
    my $ps = [ split "\t" ];
    next unless @$ps == 19;
    my ($qId, $score) = @$ps[1,16];
    if(!exists $h->{$qId}) {
        $h->{$qId} = [$score, $ps];
    } elsif($score > $h->{$qId}->[0]) {
        $h->{$qId} = [$score, $ps];
    } elsif($score == $h->{$qId}->[0]) {
        push @{$h->{$qId}}, $ps;
    }
}

for my $qId (sort(keys(%$h))) {
    my ($score, @pss) = @{$h->{$qId}};
    for my $ps (@pss) {
        print $fho join("\t", @$ps)."\n";
    }
}
close $fhi;
close $fho;


__END__
