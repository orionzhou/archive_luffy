#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  seqextract.pl - extract given sequences (with specified ranges) from an input fasta file

=head1 SYNOPSIS
  
  seqextract.pl [-help] [-in input-fasta] [-out output-fasta] [-name name-file] <IDs>

  Options:
      -help   brief help message
      -in     input file (can be 'stdin')
      -out    output file (can be 'stdout')
      -name   a file containing sequence names (and/or ranges)

=head1 DESCRIPTION

  This program extracts sequences with specified ranges from the input fasta file

=head1 OPTIONS

=over 6
  
=item B<-help>
  
  Print a usage summary.

=item B<IDs>
  
  Fasta IDs

=item B<output-file>

  To write to stdout, the user could either specify 'stdout' or simply leave this
  augument empty.

=back
  
=head1 BUGS
  
=head1 REFERENCES
  
=head1 VERSION
  
  0.1
  
=cut
  
#### END of POD documentation.
#-----------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::Seq;
use Bio::SeqIO;

my ($fi, $fo, $fn) = ('') x 3;
my ($fhi, $fho);
my $help_flag;

#----------------------------------- MAIN -----------------------------------#
GetOptions(
    "help|h"   => \$help_flag,
    "in|i=s"   => \$fi,
    "out|o=s"  => \$fo,
    "name|n=s" => \$fn,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

my @ids = @ARGV;
if ($fi eq "stdin") {
    $fhi = \*STDIN;
} else {
    open ($fhi, $fi) || die "Can't open file $fi: $!\n";
}

if ($fo eq "stdout") {
    $fho = \*STDOUT;
} else {
    open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

my $h = { map {$_=>[['', '']]} @ids };
if($fn && -s $fn) {
    open(FHN, "<$fn") or die "cannot read $fn\n";
    my $flag_first = 1;
    while(<FHN>) {
        chomp;
        next if /^\#/;
        if($flag_first && /^id/i) {
            $flag_first = 0;
            next;
        }
        my @ps = split "\t";
        next if @ps == 0;
        my $id = $ps[0];
        my ($beg, $end) = ('', '');
        ($beg, $end) = @ps[1..2] if @ps >= 3;
        $beg <= $end || die "$id: $beg - $end\n";
        push @{$h->{$id}}, [$beg, $end];
    }
    close FHN;
}

my $seqHI = Bio::SeqIO->new(-fh=>$fhi, -format=>'fasta');
my $seqHO = Bio::SeqIO->new(-fh=>$fho, -format=>'fasta');
my $cnt = 0;
while(my $seqO = $seqHI->next_seq()) {
    my ($id, $len, $seq) = ($seqO->id, $seqO->length, $seqO->seq);
    if(exists $h->{$id}) {
        for (@{$h->{$id}}) {
            my ($beg, $end) = @$_;
            if($beg && $end) {
                $beg <= $len || die "$id\[len=$len] but you want $beg-$end\n";
                $end <= $len || die "$id\[len=$len] but you want $beg-$end\n";
                my $idN = join("-", $id, $beg, $end);
                my $seqN = substr($seq, $beg-1, $end-$beg+1);
                $seqHO->write_seq( Bio::Seq->new(-id=>$idN, -seq=>$seqN) );
            } else {
                $seqHO->write_seq($seqO);
            }
            $cnt ++;
        }
    }
}
$seqHI->close();
$seqHO->close();
printf "  %4d sequences extracted\n", $cnt;

exit 0;



