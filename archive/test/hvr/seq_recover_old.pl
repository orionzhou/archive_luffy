#!/usr/bin/perl -w
use strict;
use Bio::Seq;
use Bio::SeqIO;
my $in = Bio::SeqIO->new(-file => 'E:/Scripts/test/hvr/CRS.gb',
				-format => 'genbank');
my $dir = "E:/Scripts/test/hvr/";
my $dir_out = "E:/Genius/Groups/Pop_Genetics/Mito_HVR/hvr1/";
my @pop_arr = ("Tu", "Yugu", "Sala", "Baoan", "Dongxiang", "Han_Xian");

my $CRS = $in->next_seq();
my $HVR1 = $CRS->subseq(16090, 16400);

for(my $i=6; $i<=6; $i++)
{
  open( XLS, $dir."sup".$i.".txt" ) or die("Cannot open XLS-file\n");
  my $out = Bio::SeqIO->new(-file => ">".$dir_out.$pop_arr[$i-1].".fasta",
  	    			-format => 'fasta');
  my $seq_after = Bio::Seq->new( -display_id => "Cambridge Refrence Sequence",
    	-seq => $HVR1,
	-desc => 'Homo Sapiens Mitochondrial control region hypervariable region I' );
  $out->write_seq( $seq_after );
  my $buffer = <XLS>;
  chomp( $buffer );
  my @loc_arr = split( "\t", $buffer );
#  print scalar( @ele_arr ),"\n";
#  foreach my $ele (@ele_arr)
#  {
#    print &hvr_seq($ele-16090+1),"-";
#  }
#  print "\n";

  $buffer = <XLS>;
#  chomp( $buffer );
#  @ele_arr = split( "\t", $buffer );
#  print scalar( @ele_arr ),"\n";
#  foreach my $ele (@ele_arr)
#  {
#    if($ele eq '')
#    {
#      $ele = " ";
#    }
#    print $ele,"-";
#  }
#  print "\n";

  my $sample_count = 0;
  while(<XLS>)
  {
    chomp;
    my @ele_arr = split("\t",$_,-1);
    #print scalar( @ele_arr ),"\n";
    my $j = 0;
    my $local_hvr = $HVR1;
    my @local_hvr_arr = split("",$local_hvr);
    foreach my $ele (@ele_arr)
    {
      if($ele ne "")
      {
        #print $loc_arr[$j].":".&hvr_seq($loc_arr[$j]-16090+1)."->".$ele."\t";
        my $op = 0;
        my $ele_after = "";
        if($ele =~ /[ATCGN]/)
        {
          $op = 1;
          $ele_after = $ele;
          if($loc_arr[$j] !~ /^\d+$/)
          {
            $op = 3;
            $ele_after = $local_hvr_arr[int($loc_arr[$j]-16090)].$ele;
            print $loc_arr[$j].":".$ele_after."\t";
          }
        }
        elsif($ele eq "d")
        {
          $op = 2;
          $ele_after = "";
        }
        if($op == 0)
        {
          print "\nInvalid input $ele at ".$loc_arr[$j]."!";
          exit;
        }
        $local_hvr_arr[int($loc_arr[$j]-16090)] = $ele_after;
      }
      $j ++;
    }
    $sample_count ++;
    $seq_after = Bio::Seq->new( -display_id => $pop_arr[$i-1]."_".$sample_count,
    	-seq => join("",@local_hvr_arr),
	-desc => 'Homo Sapiens Mitochondrial control region hypervariable region I' );
    $out->write_seq( $seq_after );
    print "\t$sample_count done\n";
  }
  print $pop_arr[$i-1]." - ".$sample_count." samples\n";
}

sub hvr_seq()
{
  my ($loc) = @_;
  return substr($HVR1, $loc-1, 1);
}