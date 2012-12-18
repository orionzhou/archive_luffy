#!/usr/bin/perl -w
use strict;
use Bio::Seq;
my $search_dir = "F:/Genius/STR/perl+php+fasta/handled_data";
opendir ( MYDIR , $search_dir ) or die("could not open search DIR");
my @dir_array = readdir(MYDIR);
my @dirs = ();
foreach my $dir (@dir_array)
{
	if($dir =~ /^[0-9]{3}_NT_[0-9]+\.[0-9]{1,2}_part[0-9]{4}\([0-9\-]+\)\.txt/ )
	{
		#print ($dir,"\n");
		push @dirs,$dir;
	}
}
closedir(MYDIR);
print ("There are totally ",scalar(@dirs)," files to handle, so, just be patient...\n\n");

my $file_path = "F:/Genius/STR/perl+php+fasta/for_blast/blast_1.txt";
open ( MYFILE , $file_path ) or die("could not open the file to search");
my $j=1;
while ( <MYFILE> )
{
   chomp;
   my @line_arr = split( /\n/ , $_);
   open ( RESULT, '>>F:/Genius/STR/perl+php+fasta/result/blast_result.txt');
   foreach my $line (@line_arr)
   {
        my @ele_arr = split( /\t/ , $line );
        if(scalar(@ele_arr) < 3)
        {
        	print "Line $j is abandoned, cause there's no infomation actually.\n";
                print RESULT join("\t",$ele_arr[0],"abandoned"),"\n";
        }
        elsif( scalar(@ele_arr) ==3 )
        {
        	my $seq1_str = $ele_arr[1], my $seq2_str = $ele_arr[2];
                my $seq1 = Bio::Seq->new ( -display_id => 'Primer1', -seq => $seq1_str );
                my $seq2 = Bio::Seq->new ( -display_id => 'Primer2', -seq => $seq2_str );
                my $seq1_c = $seq1->revcom(), my $seq2_c = $seq2->revcom();
                $ele_arr[3] = $seq1_c->seq(), $ele_arr[4] = $seq2_c->seq();
        	print $j,": ";
                print join("  ",@ele_arr),"\n";
                print RESULT $ele_arr[0],"\t";
                foreach my $fname (@dirs)
                {
                	open ( SEARCH_FILE, $search_dir."/".$fname ) or die("Cannot open file");
                        my $seq = '' ;
                        while( <SEARCH_FILE> )
                        {
                        	chomp;
                        	$seq .= $_;
                        }
                        $seq =~ s/\W//ig;
                        study($seq);
                	if ( $seq =~ /(($ele_arr[1])[atcg]+($ele_arr[4]))/ig
                             || $seq =~ /(($ele_arr[4])[atcg]+($ele_arr[1]))/ig
                             || $seq =~ /(($ele_arr[3])[atcg]+($ele_arr[2]))/ig
                             || $seq =~ /(($ele_arr[2])[atcg]+($ele_arr[3]))/ig )
                        {
                        	#print join("%%",$2,$3),"\n";
                                my $Target = $1;
                                my $end = pos($seq);
                                my $length = length($Target);
                                my $start = pos($seq) - $length + 1;
                                open ( DYS, ">>F:/Genius/STR/perl+php+fasta/reported_STRs/"
                                	.sprintf("%03.0f",$j)."_".$ele_arr[0].".txt" )
               				or die("Could not create result file");
                                &local2global( $fname, $start, $length, $Target );
                        }
                }
                print "\n";
                print RESULT "\n";
        }
   $j++;
   }
}
sub local2global
{
	my ( $fname, $start_level1, $length, $target ) = @_;
	my ($contig) = $fname =~ /(NT_[0-9]+\.[0-9]{1,2})/;
        my ($part) = $fname =~ /part([0-9]{4})/;
        my $start_contig = ($part-1) * 400 * 70 + 1 ;
        my $start_level2 = $start_level1 + $start_contig - 1 ;
        my %loc_contig = ( 'NT_113967.1'=>1, 'NT_113968.1'=>84822, 'NT_113969.1'=>201385,
        'NT_113970.1'=>1017558, 'NT_113971.1'=>1104114, 'NT_113972.1'=>1274235,
        'NT_113973.1'=>2128239, 'NT_011896.9'=>2709521, 'NT_086998.1'=>9024956,
        'NT_011878.9'=>9901323, 'NT_087001.1'=>11214554, 'NT_113819.1'=>11653955,
        'NT_011875.11'=>12308579, 'NT_011903.12'=>22360817, 'NT_025975.2'=>57228750,
        'NT_091573.1'=>57377045, 'NT_113974.1'=>57443438);
        my $start_level3 = $start_level2 + $loc_contig{$contig} - 1 ;
        my $end_level3 = $start_level3 + $length - 1 ;
        print "\t Gocha in $contig !\n";
        print "\t Start_level1:$start_level1, Length:$length","\n";
        print "\t Start_level2:$start_level2, Start_contig:$start_contig\n";
        print "\t Start_level3:$start_level3 (Absolute Position)\n";
        print RESULT join("\t","succeeded",$contig,$part,$start_level1,$start_level2,$start_level3,$end_level3),"\t";
        print DYS ">",$contig,"|part",$part,"|",$start_level3,"-",$end_level3,"\n";
        print DYS $target,"\n","\n";
        close( DYS );
}