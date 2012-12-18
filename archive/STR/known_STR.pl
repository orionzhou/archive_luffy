#!/usr/bin/perl -w
use strict;
#motif-repeat parameters:
#specify motif length, minimum number of repeats.
my @specs = ( [2,9], [3,6], [4,5] ,[5,5] );
my $mydirectory = "F:/Genius/STR/perl+php+fasta/reported_STRs";
opendir ( MYDIR , $mydirectory ) or die("could not open that DIR");
my @dir_array = readdir(MYDIR);
my @dirs = ();
foreach my $dir (@dir_array)
{
	if($dir =~ /^([0-9]{3})_[a-zA-Z0-9]+\.txt$/ )
	{
		#print ($dir,"\t");
		push @dirs,$dir;
	}
}
closedir(MYDIR);
@dirs = sort(@dirs);
#print ("There are totally ",scalar(@dirs)," files to handle, so, just be patient...\n\n");
my $count = 0; my $j=1;
open( RESULT, ">F:/Genius/STR/perl+php+fasta/reported_STRs/result.txt")
         or die ("could not create result file");
foreach my $current_file (@dirs)
{
	open( FILE, $mydirectory."/".$current_file)
	 	or die("could not open corresponding file");
        my $content="";
	while(<FILE>)
	{
        	$content .= $_ ;
	}
        my @paragraph_arr = split(/\n\n/,$content);
        my @seq_arr; my @name_arr; my $count=0;
        foreach my $paragraph (@paragraph_arr)
        {
        	if($paragraph)
                {
        		my @line_arr = split(/\n/,$paragraph,2);
                	$name_arr[$count] = $line_arr[0];
                	$seq_arr[$count] = $line_arr[1];
                	$count++;
                }
        }
        print $current_file,": ",$count," sequence\n";
        my @id_arr = split(/_/,$current_file);
        print RESULT join("\t",@id_arr,$count),"\t";
        my $i = 0;
        foreach my $name (@name_arr)
        {
        	my ($start_seq) = $name =~ /([0-9]+)\-[0-9]+$/ ;
                print "Seq ",($i+1),": from $start_seq\n";
                &find_STR($seq_arr[$i],$start_seq);
                $i++;
        if($id_arr[1] eq "DYS390.txt" )
        {
        	print $seq_arr[$i-1];
        }
        }
        print "\n";
        print RESULT "\n";
        $j++;
}

sub find_STR
{
	my ($sequence,$start_seq) = @_;
	$sequence =~ s/\W//g; #concatenate multi-line sequence
	study($sequence);     #This is absolutely necessary?
	#my $seqlength = length($sequence);
	my $ssr_number = 1;   #track multiple ssrs within a single sequence
	my %locations;        #track location of SSRs as detected
	my $i;
	#print STDOUT join("\t","SSR_id","Length","Motif","Repeats","Start","End","Start1","End1"),"\n";
	for($i=0; $i<scalar(@specs); $i++)
	{
		#test each spec against sequence
		my $motiflength = $specs[$i]->[0];
		my $minreps = $specs[$i]->[1] - 1;
		my $regexp = "(([gatc]{$motiflength})\\2{$minreps,})";
		while ($sequence =~ /$regexp/ig)
		{
			my $motif = uc($2); my $ssr = $1;
			if (&homo_or_pseudo_di($motif,$motiflength)==0)  #comment out this line to report monomers
			{
				my $ssrlength = length($ssr);          #overall STR length
				my $repeats = $ssrlength/$motiflength; #number of rep units
		       		my $end = pos($sequence);              #where STR ends
				pos($sequence) = $end - $motiflength;  #see docs
				my $start = $end - $ssrlength + 1;     #where STR starts
                                my $start_global = $start + $start_seq - 1;
                                my $end_global = $end + $start_seq - 1;
				if (&novel($start, \%locations))
		       		{
					#count STR only once
					print STDOUT join("\t", $ssr_number++, $motiflength,
                                         $motif, $repeats ,$start_global, $end_global, $start, $end), "\n";
                                        print RESULT join("\t",$motif,$repeats,$start_global,$end_global),"\t";
				}
			}
		}
	}
}
sub homo_or_pseudo_di
{
	#return true if motif is repeat of single nucleotide
	my ($motif,$motiflength) = @_;
	my ($reps) = $motiflength - 1;
	return 1 if ($motif =~ /([gatc])\1{$reps}/i || $motif =~ /^([gatc][gatc])\1{1,}$/i );
	return 0;
}
sub novel
{
	my($position, $locationsref) = @_;
	if(defined $locationsref->{$position})
	{
		return undef;
	}
	else
	{
		$locationsref->{$position} = 1;
		return 1;
	}
}