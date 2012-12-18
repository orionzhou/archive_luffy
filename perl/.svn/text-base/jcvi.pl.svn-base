#!/usr/bin/perl -w
use strict; use Init; use Common; use Bio::Tools::GFF; use Localdb; use Path::Class;
use Parser; use Gff; use Crp; use Mapping; use GeneModel; use Writefile; use Readfile;
use List::Util qw/min max sum/; use POSIX qw/ceil floor/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
my ($feDb, $refDb) = ("mt_35_crp_jcvi", "mt_35");
my $dirW = dir($DIR_Misc2, "jcvi");
my $f01 = file($dirW, "01_id_mapping.txt");
my $f02 = file($dirW, "02_raw.gff");
my $f03 = file($dirW, "03.gff");
#gff_conv_jcvi(-mapping=>$f01, -in=>$f02, -out=>$f03);
my $f04 = file($dirW, "04.gff");
#gff_conv_loc(-in=>$f03, -out=>$f04, -db=>$refDb);
my $f11 = file($dirW, "11.gtb");
#gff2Gtb(-in=>$f04, -out=>$f11);
my $f12 = file($dirW, "12.gff");
#gtb2Gff(-in=>$f11, -out=>$f12);
my $ld = Localdb->new(-db=>$feDb);
#$ld->loadGff(-files=>$f12, -empty=>1);



