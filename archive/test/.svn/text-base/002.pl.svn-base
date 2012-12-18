#!/usr/bin/perl
use Bio::Seq;
use Bio::SeqIO;
use strict;

#use Bio::DB::GenBank;
#my $database = new Bio::DB::GenBank;
#my $seq = $database->get_Seq_by_gi('405830');

#my $out = Bio::SeqIO->newFh ( -file => '>002.gbk', -format=>'genbank' );
my $out = Bio::SeqIO->new( -file => '>test/002.fasta', -format => 'fasta' );
my $seq1 = Bio::Seq->new( -display_id => 'AY88888', -seq => 'AAAAATTTTT',
	-desc => 'Homo Sapiens Mitochondrion DNA' );
print join("\t",$seq1->id,":",$seq1->length,$seq1->desc),"\n";
#print $out $seq;
$out->write_seq($seq1);