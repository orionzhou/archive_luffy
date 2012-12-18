#!/usr/bin/perl -w
use strict; use Init; use Common; use Run; use Localdb; use Readfile; use Path::Class;
use Parser; use Gff; use Crp; use Mapping;
use List::Util qw/min max sum/; use POSIX qw/ceil floor/;
my ($feDb, $refDb) = ("mt_35", "mt_35");
my $f00 = file($DIR_Misc2, "crp", "06_groups.txt");

my $DIR_Work = dir($DIR_Misc2, "crp", "30_blast");
$DIR_Work = dir($DIR_Misc2, "crp", "31_neighbor");
my $type = "up1k";
my $f01 = file($DIR_Genome, $refDb, "01_seq_$type.fa");
#writeSeqByOpt(-type=>$type, -out=>$f01, -opt=>1, -fedb=>$feDb, -refdb=>$refDb);
my $f02 = file($DIR_Work, "02_ncr_$type.fa");
#writeSeqByOpt(-type=>$type, -out=>$f02, -opt=>2, -fedb=>$feDb, -refdb=>$refDb, -in=>$f00);
my $f03 = file($DIR_Work, "03_blast_out.txt");
#run_blast(-in=>$f02, -out=>$f03, -program=>'blastn', -db=>'mt_35_up1k', -aln=>0);
my $f04 = file($DIR_Work, "04_blast_filtered.txt");
#blastFilter(-in=>$f03, -out=>$f04, -strand=>1);
my $f11 = file($DIR_Work, "11_candidates.txt");
#getCandidates(-in=>$f00, -opt=>1, -fedb=>$feDb, -fblast=>$f04, -out=>$f11);
#getCandidates(-in=>$f00, -opt=>2, -fedb=>$feDb, -out=>$f11);

my $f12 = file($DIR_Work, "12_homology.txt");
my $d13 = dir($DIR_Work, "13_aln");
#findHomology(-in=>$f11, -out=>$f12, -alndir=>$d13, -fedb=>$feDb, -refdb=>$refDb);
my $f14 = file($DIR_Work, "14_homology.txt");
#addExp(-in=>$f12, -out=>$f14, -ds=>'rnaseq', -opt=>[['qId', 2], ['hId', 5]]);

$DIR_Work = dir($DIR_Misc2, "crp", "40_dotplot");
my $f41 = file($DIR_Work, "01_anchors.txt");
my $f42 = file($DIR_Work, "02_loc.txt");
#anchorToLoc(-in=>$f41, -out=>$f42, -fedb=>$feDb, -refdb=>$refDb);
doCombo(-in=>$f42, -dir=>$DIR_Work);

