#!/usr/bin/perl
use Bio::DB::GenBank;
use Bio::DB::Query::GenBank;
use Bio::Seq;
use Bio::SeqIO;
use strict;
my $query = "homo sapiens[Organism] AND dystrophin[Title] AND RefSeq[Property]";
my $result = Bio::DB::Query::GenBank->new (-db=>'nucleotide',
	-query=>$query, -mindate =>'2005') or die("Query Failed");
my $count = $result->count;
print "Server returned you ",$count, " items.\n";
my @id_arr = $result->ids;
print "Their ids are :",join($",@id_arr),"\n";
my $gb = new Bio::DB::GenBank;
my $stream = $gb->get_Stream_by_id([sprintf(join(",",@id_arr))]);
my $out = Bio::SeqIO->newFh(-file=>'>>test/003.out',-format=>'genbank');

while(my $seq = $stream->next_seq)
{
	#print $stream->length;
	print $seq->id,"\t",$seq->length;
	print $out $seq;
        print "\n";
}
