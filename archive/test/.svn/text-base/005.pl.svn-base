#!/usr/bin/perl -w
use strict;
use Bio::Tools::Run::StandAloneBlast;
use Bio::SeqIO;
my $seq1_in = Bio::SeqIO -> new (-file => 'D:/software/blast/for_blast/chrY.txt',
	-format => 'fasta' );
#my $seq2_in = Bio::SeqIO -> new (-file => 'D:/software/blast/db/hs_ref_chrY.fa',
#	-format => 'fasta' );
my $seq2_in = Bio::SeqIO -> new (-file => 'D:/software/blast/for_blast/chrY1.txt',
	-format => 'fasta' );
my $seq1 = $seq1_in -> next_seq();
my $seq2 = $seq2_in -> next_seq();
my @params = ( program => 'blastn',
	'outfile' => 'D:/software/blast/blast_result/1v10.txt' );
my $factory = Bio::Tools::Run::StandAloneBlast->new (@params);
my $m = 'PAM100' ;
$factory->M($m);
$factory->bl2seq( $seq1, $seq2 );
print "Okay! Now check your result!\n";
