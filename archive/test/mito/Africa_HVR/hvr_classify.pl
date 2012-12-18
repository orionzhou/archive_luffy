#!/usr/bin/perl -w
my $dir = "E:/Scripts/test/Africa_HVR/Seq";
opendir( SEQDIR, $dir ) or die("Could not open seq_dir\n");

while(my $entity = readdir(SEQDIR))
{
  if( -f $dir."/".$entity )
  {
    #print $entity,"\n";
    open( SEQ, $dir."/".$entity ) or die("could not read file");
    my $buffer = "";
    read( SEQ, $buffer, 2000 );
    $buffer =~ /LOCUS\s+([A-Z0-9]+)/;
    my $ACC = $1;
    $buffer =~ /GI:(\d+)\n/;
    my $gi = $1;
    my %item = ("haplotype",0,"country",1,"isolate",2,"isolation_source",3,
  	  "note",4);
    my @records = ("","","","","");
    while($buffer =~ m|(\s{10,}\/(.+))\n|g)
    {
      my $line = $2;
      if($2 =~ /(country|isolate|note|haplotype|isolation_source)=\"(.*)\"/)
      {
        if( substr($2,0,6) ne "codons" )
        {
          #print "$1:$2\t";
          $records[$item{$1}] = $2;
        }
      }
    }
    print join("\t",$ACC,$gi,@records),"\n";
  }
}