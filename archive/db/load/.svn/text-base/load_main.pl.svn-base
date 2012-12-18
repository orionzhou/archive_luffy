#!/usr/bin/perl -w
#根据in文件内的序列文件，对mt_main表进行插入，可以自动提取NCBI输出的混合id中的Acc No.和gi号，插入时以Accession Number为主要依据决定是插入还是更新序列数据
unshift (@INC, "/home/orion/Scripts/DB/INC");
require "db_connect.pl";

use strict;
use Bio::Seq;
use Bio::SeqIO;

our $dbh;
my $working_dir = "./";
my $seq_file = "in_seq";

my $seq_obj = Bio::SeqIO->new(-file=>$working_dir.$seq_file, -format=>'fasta')
	or die("Cannot open sequence file");
while(my $seq = $seq_obj->next_seq()) {
	$seq->id =~ /\|([A-Z]{2}\d{4,8})\.(\d+)\|/;
	if(!$1 || !$2 ) {
		die("Error in resolving acc_number from ".$seq->id."\n");
	} else {
		my $acc = $1;
		my $version = $2;
		$seq->id =~ /gi\|(\d{7,10})\|/;
		my $gi = $1;
		&insert_acc($acc,$version,$gi,$seq);
	}

}
$dbh->disconnect();

sub insert_acc() {
	my ($acc,$version,$gi,$seq) = @_;
	my $query = "select acc,version,id from mt_main where acc=\"$acc\"";
	my $sth = $dbh->prepare($query);
	$sth->execute();
	my @row = $sth->fetchrow_array();
	$sth->finish();
	if(!@row) {
		#do insert
		my $query = "insert into mt_main (acc,version,gi,length,seq) VALUES (\"$acc\", \"$version\", \"$gi\", ".$seq->length().", \"".$seq->seq()."\")";
		my $rows = $dbh->do($query);
		if(!defined($rows)) {
			die("Error in inserting $acc\n");
		} else {
			print "$acc inserted\n";
		}
	} elsif($row[1] < $version) {
		#do update
		$query = "update mt_main set version=\"$version\", gi=\"$gi\", length=".$seq->length().", seq=\"".$seq->seq()."\" where acc=\"$acc\"";
		my $rows = $dbh->do($query);
		if(!defined($rows)) {
			die("Error in updating $acc\n");
		} else {
			print "$acc updated\n";
		}
	} else {
		print "$acc - no need to update\n";
	}
}