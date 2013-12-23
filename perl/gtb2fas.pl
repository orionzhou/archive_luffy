#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gtb2fas.pl - convert a Gtb file to FASTA format

=head1 SYNOPSIS
  
  gtb2fas.pl [-help] [-seq refseq-fasta] [-opt option] [-in input-file] [-out output-file]

  Options:
      -help   brief help message
      -in     input file
      -out    output file
      -seq    reference sequence file
      -opt    ouput option (default: protein)

=cut
  
#### END of POD documentation.
#-----------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Location;
use Seq;

my ($fi, $fo, $fs) = ('') x 3;
my $opt = "protein";
my $help_flag;

#----------------------------------- MAIN -----------------------------------#
GetOptions(
    "help|h"   => \$help_flag,
    "in|i=s"   => \$fi,
    "out|o=s"  => \$fo,
    "seq|s=s"  => \$fs,
    "opt|p=s"  => \$opt,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fs;

my ($fhi, $fho);
if ($fi eq '' || $fi eq "stdin" || $fi eq "-") {
    $fhi = \*STDIN;
} else {
    open ($fhi, $fi) || die "Can't open file $fi: $!\n";
}

if ($fo eq '' || $fo eq "stdout" || $fo eq "-") {
    $fho = \*STDOUT;
} else {
    open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

my $seqH = Bio::SeqIO->new(-fh=>$fho, -format=>"fasta");
while( <$fhi> ) {
    chomp;
    next if /(^id)|(^\#)|(^\s*$)/;
    my $ps = [ split "\t" ];
    next unless @$ps >= 18;
    my ($id, $par, $chr, $beg, $end, $srd, $locES, $locIS, $locCS, $loc5S, $loc3S, $phaseS, $src, $conf, $cat1, $cat2, $cat3, $note) = @$ps;
    $cat2 eq "mRNA" || next;
    $locCS || die "no CDS for $id\n";
    my $rloc = locStr2Ary($locCS);
    my $loc = $srd eq "-" ? [map {[$end-$_->[1]+1, $end-$_->[0]+1]} @$rloc] : 
        [map {[$beg+$_->[0]-1, $beg+$_->[1]-1]} @$rloc];

    my $seq;
    if($opt =~ /^cds$/i) {
        my $seqstr = seqRet($loc, $chr, $srd, $fs);
        $seq = Bio::Seq->new(-id=>$id, -seq=>$seqstr);
    } elsif($opt =~ /^pro/i ) {
        my $seqstr = seqRet($loc, $chr, $srd, $fs);
        my @phases = split(",", $phaseS);
        $seq = Bio::Seq->new(-id=>$id, -seq=>$seqstr)->translate(-frame=>$phases[0]);
    } elsif($opt =~ /^mrna$/) {
        $loc = [[$beg, $end]];
        my $seqstr = seqRet($loc, $chr, $srd, $fs);
        $seq = Bio::Seq->new(-id=>$id, -seq=>$seqstr);
    } elsif($opt =~ /^mrna\+$/) {
        $loc = [[$beg-1000, $end+1000]];
        my $seqstr = seqRet($loc, $chr, $srd, $fs);
        $seq = Bio::Seq->new(-id=>$id, -seq=>$seqstr);
    } else {
        die "unknown opt: $opt\n";
    }
    $seqH->write_seq($seq);
}
close $fhi;
$seqH->close();

__END__
