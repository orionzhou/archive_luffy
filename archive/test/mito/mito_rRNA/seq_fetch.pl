#!/usr/bin/perl -w
use strict;
use Bio::Seq;
use Bio::SeqIO;

my $dir = "E:/Scripts/test/mito_tRNA/seq";
my $unfound = 0;
&traverse_dir( $dir, 0 );

sub seq_read
{
  my ( $dir, $file ) = @_;
  my $path = $dir."/".$file;
  if( $file =~ /([A-Z]{2}\d{5,6}).*\.txt/ )
  {
    open( SEQ, $path ) or die("could not open seq_file:\n\t$path\n");
    my $AC = $1;
    print "\t-\t$AC\n";
    my $buffer = "";
    read( SEQ, $buffer, 20000);
    #$buffer =~ s/\s//;
    if( $buffer =~
  	/(AATAGGTTT)([A-Z]+)(GGACGAAC).+(GCTAAACCT)([A-Z]+)(AACAGGGTT)/ )
    {
      #print "Pattern found\n";
      my $out = Bio::SeqIO->new( -file => ">".$dir."/".substr($file,0,length($file)-4)
      .".fasta", -format => 'fasta' );
      my $seq1 = Bio::Seq->new( -display_id => $AC, -seq => $1.$2.$3,
      	-desc => 'Homo Sapiens Mitochondrion 12s rRNA' );
      my $seq2 = Bio::Seq->new( -display_id => $AC, -seq => $4.$5.$6,
      	-desc => 'Homo Sapiens Mitochondrion 16s rRNA' );
      $out->write_seq($seq1);
      $out->write_seq($seq2);
      return 0;
    }
    else
    {
      print substr($path,30,length($path)-30).": pattern not found!\n";
      return 1;
    }
  }
  else
  {
    print "\n";
    return 0;
  }
}

sub traverse_dir()
{
  my ($dir,$level) = @_;
  opendir( SEQ, $dir ) or die("Could not open seq_dir:\n\t$dir\n");
  my $entity;
  my $dircount = 0;
  my @dir_arr = ();
  while( $entity = readdir(SEQ) )
  {
    if( -d $dir."/".$entity )
    {
      if( $entity ne "." && $entity ne ".." )
      {
        push @dir_arr, $entity;
      }
      $dircount ++;
    }
    elsif( -f $dir."/".$entity )
    {
      print "\t" x $level,"File\t-\t",$entity;
      if(seq_read( $dir, $entity ) == 1)
      {
        $unfound ++;
      }
    }
  }
  if( $dircount == 2 )
  {
    return;
  }
  else
  {
    foreach my $dir_to_access (@dir_arr)
    {
      print "\t" x $level,"Directory\t-\t",$dir_to_access,"\n";
      &traverse_dir( $dir."/".$dir_to_access, $level+1 );
    }
  }
}

print "\nTotally $unfound seqs unfound\n";