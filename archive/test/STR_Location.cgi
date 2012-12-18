#!/usr/local/bin/perl -w
use CGI;
use Bio::Seq;
use Bio::SeqIO;
use strict;
my $obj = new CGI;
print $obj->header(-type=>"text/html; charset=gb2312") ;
print $obj->start_html( -title=>"STR_Location", -BGCOLOR=>"#DFFCFD");
&print_form($obj);
&print_result($obj) if($obj->param);
print $obj->end_html();

sub print_form
{
	my ($obj) = @_;
	print $obj->h2({-align=>"center"},"For retrieving STR and it's flanking sequence");
	print $obj->startform(-method=>'get');
	print "<table align='center' border=1 cellpadding=2 cellspacing=2 width='80%' bgcolor='#ECF9FF'>",
		"<tr align=center style='font-size:18px; color:#029EE1'><td width=50%>";
	print "Start_Position: ";
	print $obj->textfield(-name=>'pos_start',-default=>'',-size=>8,-maxlength=>8),$obj->br();
	print "</td><td>End_Position: ";
	print $obj->textfield(-name=>'pos_end',-default=>'',-size=>8,-maxlength=>8),$obj->br();
	print "</td></tr><tr align=center style='font-size:15px'>",
		"<td colspan=2 height=35 valign=middle>Show flanking sequences of&nbsp;&nbsp;";
	my %labels = ('300'=>'300 bp','350'=>'350 bp','400'=>'400 bp','450'=>'450 bp','500'=>'500 bp');
	print $obj->popup_menu( -name=>'flanking_len', -value=>['300','350','400','450','500'],
		-default=>'300', -labels=>\%labels );
	print "</td></tr><tr align=center><td colspan=2>";
	print $obj->submit('action','Submit');
	print "</td></tr></table>";
        print $obj->hr();
}

sub print_result
{
	my ($obj) = @_;
        my $pos_start = $obj->param("pos_start");
        my $pos_end = $obj->param("pos_end");
        my $flanking_len = $obj->param("flanking_len");
        if( $pos_start !~ /^[0-9]{1,8}$/ || $pos_end !~ /^[0-9]{1,8}$/ )
        {
        	print $obj->h3({-align=>"center"},"Are you kidding me!");
                exit;
        }
        elsif( $pos_start>=$pos_end || abs($pos_start-$pos_end)<15 || abs($pos_start-$pos_end)>400 )
        {
        	print $obj->h3({-align=>"center"},"Please input a valid number!");
                exit;
        }
        my @contig_name = ( 'NT_113967.1', 'NT_113968.1', 'NT_113969.1', 'NT_113970.1',
        	'NT_113971.1', 'NT_113972.1', 'NT_113973.1', 'NT_011896.9', 'NT_086998.1',
                'NT_011878.9', 'NT_087001.1', 'NT_113819.1', 'NT_011875.11', 'NT_011903.12',
                'NT_025975.2', 'NT_091573.1', 'NT_113974.1' );
        my @contig_start = ( 1, 84822, 201385, 1017558, 1104114, 1274235, 2128239, 2709521,
        	9024956, 9901323, 11214554, 11653955, 12308579,
                22360817, 57228750, 57377045, 57443438);
        my @contig_length = ( 34821, 86563, 766173, 36556, 80121, 754004, 581282, 6265435,
        	276367, 813231, 39401, 554624, 10002238, 4867933, 98295, 66393, 329517 );
        my $tag = -1 ;
        my $i = 0 ;
        foreach my $start (@contig_start)
        {
        	my $length = $contig_length[$i];
                if( $pos_start >= $start && $pos_end <= ($start+$length-1) )
                {
                	$tag = $i;
                }
                $i++;
        }
        if($tag == -1)
        {
        	print $obj->h3({-align=>"center"},"已测序的Y染色体序列尚不包括您输入的范围!");
                exit;
        }
        my $pos_start_1 = $pos_start - $contig_start[$tag] + 1 ;
        my $length = $pos_end - $pos_start + 1 ;
        my $part = int($pos_start_1 / (400*70)) + 1 ;
        my $pos_start_2 = $pos_start_1 % (400*70);
        my $file_prefix = sprintf("%03.0f",$tag+1)."_".$contig_name[$tag]."_part".sprintf("%04.0f",$part);
        print $obj->h2({-align=>"center"},"Result");
        print "<table align='center' border=1 cellpadding=2 cellspacing=2 width='80%' bgcolor='#FDF2FC'>",
		"<tr align=center style='font-size:14px; color:#029EE1'><td colspan=2>";
        print "Located in ",$file_prefix,"\n";
        my $directory = "/var/www/html/webroot/personal/genius/Y_STR/handled_data";
        my $file_name = "";
	$file_name = $file_prefix."(".(($part-1)*400*71+1)."-".($part*400*71).").txt";
        open( SEQ, $directory."/".$file_name ) or die("Could not find the file");
        my $seq="";
        while(<SEQ>)
        {
        	chomp;
                $seq .= $_;
        }
        my $start_true = $pos_start_2  -$flanking_len - 1 ;
        my $length_true = 2 * $flanking_len + $length ;
        my $up_length_true = $flanking_len ;
        if($start_true<0)
        {
        	$length_true = $pos_start_2 + $length + $flanking_len - 1 ;
                $start_true = 0 ;
                $up_length_true = $pos_start_2 - 1 ;
        }
        my $seq_all = substr($seq, $start_true, $length_true);
        my $seq_up = substr($seq, $start_true, $up_length_true);
        my $seq_str = substr($seq, $pos_start_2-1, $length);
        my $seq_down = substr($seq, $pos_start_2+$length-1, $flanking_len);
        print "</td></tr><tr align=center>";
        print "<td width=10%>Upstream ",$flanking_len,"bp</td><td>";
        print "<textarea cols=70 rows=7 style='font-size:-1' readonly>";
        &print_seq("Upstream_".$flanking_len."bp", $seq_up, 'fasta');
        print "</textarea></td></tr><tr align=center><td>STR</td><td>";
        print "<textarea cols=70 rows=3 style='color:#0000FF; font-size:-1' readonly>";
        &print_seq("Short_Tandem_Repeat", $seq_str, 'fasta');
        print "</textarea></td></tr><tr align=center><td>Downstream ",$flanking_len,"bp</td><td>";
        print "<textarea cols=70 rows=7 style='font-size:-1' readonly>";
        &print_seq("Downstream_".$flanking_len."bp", $seq_down, 'fasta');
        print "</textarea></td></tr><tr align=center><td>FASTA</td><td>";
        print "<textarea cols=70 rows=13 style='color:#666666; font-size:-1' readonly>";
        &print_seq('All_sequence',$seq_all,'fasta');
        print "</textarea></td></tr><tr align=center><td>GENBANK</td><td>";
        print "<textarea cols=80 rows=18 style='color:#666666; font-size:-1' readonly>";
        &print_seq('All_Sequence',$seq_all,'genbank');
        print "</textarea>";
        print "</td></tr></table>";
}
sub print_seq
{
	my($id,$seq_all,$format) = @_;
        my $seq_obj = Bio::Seq->new(-display_id=>$id, -seq=>$seq_all);
        my $seq = Bio::SeqIO->new(-format=>$format);
        $seq->write_seq($seq_obj);
}