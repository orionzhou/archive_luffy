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
        if ( $seq =~ /((AGCCACCTGGGTA[atcg]{2}TGAGG)[atcg]+(CAGA[atcg]{2}GAAAAGCTGCAACA))/ig
           	|| $seq =~ /((TGTTGCAGCTTTTC[atcg]{2}TCTG)[atcg]+(CCTCA[atcg]{2}TACCCAGGTGGCT))/ig
                || $seq =~ /((CCTCA[atcg]{2}TACCCAGGTGGCT)[atcg]+(TGTTGCAGCTTTTC[atcg]{2}TCTG))/ig
                || $seq =~ /((CAGA[atcg]{2}GAAAAGCTGCAACA)[atcg]+(AGCCACCTGGGTA[atcg]{2}TGAG))/ig
                )
        {
		my $Target = $1;
                my $end = pos($seq);
                my $length = length($Target);
                my $start = pos($seq) - $length + 1;
                &local2global( $fname, $start, $length );
        }
}
sub local2global
{
	my ( $fname, $start_level1, $length ) = @_;
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
        print "\t Gocha in $contig !\n";
        print "\t Start_level1:$start_level1, Length:$length","\n";
        print "\t Start_level2:$start_level2, Start_contig:$start_contig\n";
        print "\t Start_level3:$start_level3 (Absolute Position)\n";
}
