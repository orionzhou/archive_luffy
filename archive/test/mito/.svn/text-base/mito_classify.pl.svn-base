#!/usr/bin/perl -w
my $dir = "E:/Scripts/test/mito/mito_seq/";
print join("\t","Locus","GI","haplotype","country","isolate",
	"isolation_source","note"),"\n";
for(my $i=1; $i<5; $i++)
{
  my $file_name = sprintf( "%04.0f", $i ).".txt";
  open( MITO, $dir.$file_name ) or die("could not read file");
  my $buffer = "";
  read( MITO, $buffer, 2000 );
  $buffer =~ /LOCUS\s+([A-Z0-9]+)/;
  my $locus = $1;
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
  print join("\t",$locus,$gi,@records),"\n";
}
