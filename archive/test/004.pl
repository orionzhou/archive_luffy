#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlast;
my @params = ('database' => 'drosoph.aa',
	'outfile' => 'D:/software/blast/blast_result/376.out',
        'program' => 'blastp',
        _READMETHOD => "Blast");
my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);
my $seq_in = Bio::SeqIO->new( -file=>'D:/software/blast/for_blast/kkk.txt',
	-format=>'fasta' );
my $query = $seq_in->next_seq();
$factory->blastall($query);
#my $exe = $factory->program_dir();
#print $exe,"\n";