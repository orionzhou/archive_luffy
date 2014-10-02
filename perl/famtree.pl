#!/usr/bin/perl
use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin";
use Bio::SeqIO;
use Bio::DB::Fasta;
use Common;
use Data::Dumper;
use File::Path qw/make_path remove_tree/;
use File::Basename;
use List::Util qw/min max sum/;

my @orgs = qw/HM101 HM004 HM010 HM018 HM034 HM050 HM056 HM058 
  HM060 HM095 HM129 HM185 HM324 HM340/;
my @fams = qw/crp0010 crp0110 crp0355 crp0675/;

my $dir = "$ENV{'misc2'}/genefam";
-d $dir || make_path($dir);
chdir $dir || die "cannot chdir to $dir\n";

#collect_fam_crp(\@orgs, '11.crp.tbl');

build_fam_tree('11.crp.tbl');
sub build_fam_tree {
  my ($fi) = @_;
  -d "21_seq" || make_path("21_seq");
  -d "23_aln" || make_path("23_aln");
  -d "24_phb" || make_path("24_phb");
  -d "26_fig" || make_path("26_fig");
  my $t = readTable(-in => $fi, -header => 1);
  $t->sort("family", 1, 0, "org", 1, 0, "chr", 1, 0, "beg", 0, 0);

  my $ref = group($t->colRef("family"));
  for my $fam (sort(keys(%$ref))) {
    my ($idx, $cnt) = @{$ref->{$fam}};
    $fam = lc($fam);
    open(my $fho, ">21_seq/$fam.fas") or die "cannot write $fam.fas\n";
    for my $i ($idx..$idx+$cnt-1) {
      my ($id, $fam2, $seq) = map {$t->elm($i, $_)} qw/id family sequence/;
      print $fho ">$id\n$seq\n";
    }
    close $fho;
    runCmd("clustalo -i 21_seq/$fam.fas -o 23_aln/$fam.aln \\
      --outfmt=clu --force --full --full-iter");
    runCmd("clustalw2 -infile=23_aln/$fam.aln -bootstrap=1000 \\
      -outputtree=phylip -bootlabels=node -clustering=NJ -kimura");
    runCmd("mv 23_aln/$fam.phb 24_phb/$fam.phb");
  }
}
sub collect_fam_crp {
  my ($orgs, $fo) = @_;
  open(my $fho, ">$fo") or die "cannot write $fo\n";

  my $flag = 1;
  for my $org (@$orgs) {
    my $di = "$ENV{'misc4'}/spada.crp.".$org;
    my $fi = "$di/61_final.tbl";
    my $ti = readTable(-in => $fi, -header => 1);
    if($flag) {
      print $fho join("\t", "org", 'loc',  $ti->header)."\n";
      $flag = 0;
    }
     
    for my $i (0..$ti->lastRow) {
      my ($id, $fam, $chr, $beg, $end, $srd, $e, 
        $score_sp, $score_hmm, $score_aln, $seq) = $ti->row($i);
      my @ps = split("_", $id);
      my ($fam2, $chr2, $pos, $cnt) = @ps;
      my $nid = join("_", lc($org), lc($chr), $pos, $cnt);
      my $loc = "$chr:$beg-$end";
      print $fho join("\t", $org, $loc, $nid, $fam, 
        $chr, $beg, $end, $srd, $e,
        $score_sp, $score_hmm, $score_aln, $seq)."\n";
    }
  }
  close $fho;
}

