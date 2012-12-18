#!/usr/bin/perl -w
use strict;
use Bio::Seq;
use Bio::SeqIO;
my $dir_in = "E:/Scripts/test/mito/mito_seq/";
my $dir_out = "E:/Scripts/test/mito/mito_seq/Seq_fasta/";

for(my $i=580; $i<=641; $i++)
{
  my $in = Bio::SeqIO->new( -file => $dir_in."EF184".sprintf("%03.0f",$i).'.gb',
   -format=>'genbank');
  my $seq = $in->next_seq();
  my $out = Bio::SeqIO->new( -file => '>'.$dir_out."EF184".sprintf("%03.0f",$i).'.fasta',
   -format => 'fasta' );
  $out->write_seq($seq);
  print "EF184".$i.".gb complete\n";
}