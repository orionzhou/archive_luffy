#!/usr/bin/perl -w
use strict;
use Bio::Seq;
use Bio::SeqIO;
my $in = Bio::SeqIO->new(-file => '../align/rCRS.gb',
				-format=>'genbank');
#my $out = Bio::SeqIO->new(-file=> '>001.fasta', -format=>'fasta');
while (my $seq = $in->next_seq())
{
	print $seq->primary_id()."\t".$seq->display_id()."\n";
        print "first 10 bases:\t", $seq->subseq(1,10), "\n";
        print "length:\t", $seq->length, "\n";
        if( defined $seq->species )
        {
          my $species = $seq->species();
          print "Species:\t",$species->binomial," [",$species->common_name,"]\n";
        }
	my @features = $seq->get_SeqFeatures();
        #my @features = $seq->get_all_SeqFeatures();
        print "Features:\n";
	foreach my $feat (@features)
        {
	  print join("\t",$feat->primary_tag.":",$feat->start,$feat->end,
          	$feat->length,$feat->strand),"\n";
          my $feat_str = $feat->gff_string;
          print "\t",$feat_str,"\n";
          my @sub_feat_arr = $feat->get_SeqFeatures();
          foreach my $sub_feat (@sub_feat_arr)
          {
            print "\t\t",$sub_feat->primary_tag,"\n";
          }
          #print "\t",$feat->seq->seq(),"\n";
        }
        #$out->write_seq($seq);
}
