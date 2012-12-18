#!/usr/bin/perl

$/ = ">";
$| = 1;

#motif-repeat parameters:
#specify motif length, minimum number of repeats.
#modify according to researcher's preferences
my @specs = ([3,6],  #trinucl. with >= 6 repeats
	     [4,5]); #tetranucl. with >= 5 repeats

$mydirectory = "/var/www/html/webroot/personal/genius/Y_STR/handled_data";
opendir ( MYDIR , $mydirectory ) or die("could not open that DIR");
@dir_array = readdir(MYDIR);
@dirs = ();
foreach $dir (@dir_array)
{
	if($dir =~ /^[0-9]{3}_NT_[0-9]+\.[0-9]{1,2}_part[0-9]{4}\([0-9\-]+\)\.txt$/  )
	{
		#print ($dir,"\t");
		push @dirs,$dir;
	}
}
closedir(MYDIR);

#print ("There are totally ",scalar(@dirs)," files to handle, so, just be patient...\n\n");
my $seqcount = 0;$j=1;
foreach $current_file (@dirs)
{
print STDOUT ("here begins file ",$j++," : $current_file \n");
open( RESULT, ">>result.txt") or die ("could not create result file");
open( FILE, $mydirectory."/".$current_file ) or die("could not open corresponding file");
while(<FILE>)
  {
	#FASTA formatted sequences as input
	chomp;
	#my ($titleline, $sequence) = split(/\n/,$_,2);
	#next unless ($sequence && $titleline);
	my ($id) = $current_file =~ /(NT_[0-9\.]+_part[0-9]{4})/;
	my $sequence = $_;
	$seqcount++;
	#my ($id) = $titleline =~ /^(\S+)/;
	$sequence =~ s/\s//g; #concatenate multi-line sequence
	study($sequence);     #This is absolutely necessary?
	#my $seqlength = length($sequence);
	my $ssr_number = 1;   #track multiple ssrs within a single sequence
	my %locations;        #track location of SSRs as detected
	my $i;
	print STDOUT join("\t","id\t\t","SSR_id","Length","Motif","Repeats","Start","End"),"\n";
	for($i=0; $i<scalar(@specs); $i++)
	{
		#test each spec against sequence
		my $motiflength = $specs[$i]->[0];
		my $minreps = $specs[$i]->[1] - 1;
		my $regexp = "(([gatc]{$motiflength})\\2{$minreps,})";
		while ($sequence =~ /$regexp/ig)
		{
			my $motif = uc($2); my $ssr = $1;
			#reject "aaaaaaaaa", "ggggggggggg", etc.
			if (&homo_or_pseudo_di($motif,$motiflength)==0)  #comment out this line to report monomers
			{
			my $ssrlength = length($ssr);          #overall STR length
			my $repeats = $ssrlength/$motiflength; #number of rep units
			my $end = pos($sequence);              #where STR ends
			pos($sequence) = $end - $motiflength;  #see docs
			my $start = $end - $ssrlength + 1;     #where STR starts
			if (&novel($start, \%locations))
			{
				#count STR only once
				print RESULT join("\t", $id, $ssr_number, $motiflength, $motif, $repeats, $start, $end), "\n";
				print STDOUT join("\t", $id, $ssr_number++, $motiflength, $motif, $repeats, $start, $end), "\n";
			}
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
	return 1 if ($motif =~ /([gatc])\1{$reps}/i || $motif =~ /([gatc][gatc])\1{1,}/i );
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
