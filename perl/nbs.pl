#!/usr/bin/perl
use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin";
use Bio::DB::Fasta;
use Common;
use Data::Dumper;
use File::Path qw/make_path remove_tree/;
use File::Basename;
use List::Util qw/min max sum/;

my $fam = "PF00931";
my @orgs = qw/HM034 HM056 HM340.APECCA/;
my $dir = "/home/youngn/zhoup/Data/misc2/genefam/$fam";
-d $dir || make_path($dir);
chdir $dir || die "cannot chdir to $dir\n";

extract_seq_by_pfam("$dir/01.fas", \@orgs, $fam);
runCmd("clustalo -i 01.fas -o 05.aln --outfmt=clu --force --hmm-in=../$fam.hmm --full");
runCmd("aln2phy.pl -i 05.aln -o 05.phy -l 30");
runCmd("/home/youngn/zhoup/Source/PhyML-3.1/phyml -i 05.phy -d aa");
runCmd("mv 05.phy_phyml_tree.txt 06.nwk");
runCmd("rm 05.phy_phyml*");
 
my $f17 = "$dir/17.phb";
#convertTree(-in=>$f021, -out=>$f022, -informat=>'newick', -outformat=>'newick');

sub extract_seq_by_pfam {
  my ($fo, $orgs, $fam) = @_;
  open(my $fho, ">$fo") or die "cannot write $fo\n";
 
  my $di = '/home/youngn/zhoup/Data/genome';
  for my $org (@$orgs) {
    my $fd = "$di/$org/21.fas"; 
    my $fp = "$di/$org/23.pfam.tsv"; 
    my $tp = readTable(-in => $fp, -header => 0);
    my $db = Bio::DB::Fasta->new($fd);
  
    for my $i (0..$tp->nofRow-1) {
      my ($id, $md5, $len, $prog, $pid, $desc, $beg, $end, $e) = $tp->row($i);
      if($pid eq $fam) {
        my $nid = "$org-$id-$beg-$end";
        my $seq = $db->seq($id, $beg, $end);
        print $fho ">$nid\n$seq\n";
      }
    }
  }
  close $fho;
}

=cut
#run_hmmscan(-in=>$f03, -out=>$f04);
my $f05 = file($dir, "05_domain.txt");
#hmmParse(-in=>$f04, -out=>$f05);
#hmmscanTiling(-in=>$f04, -e=>1e-5, -inc=>"LRR", -out=>$f05);

my $f06 = dir($dir, "06_meme");
my $param_meme = {flag=>["protein"], mod=>"anr", nmotifs=>3, maxsize=>1000000};
#run_meme(-in=>$f03, -out=>$f06, -param=>$param_meme);

my $f11 = file($dir, "../../11_mt1_anno.txt");
my $f12  = file($dir, '12_anno.tlf');
#annForNBS(-in=>$f11, -out=>$f12, -db=>'mt_nbs', -ids=>\@ids, -opt=>2);

my $f081  = file($dir, '08_anno.txt');
my $f082  = file($dir, '08_anno.tlf');
#annProbe(-out=>$f081, -ids=>$param->{ids}, -opt=>1);
#annProbe(-out=>$f082, -ids=>$param->{ids}, -opt=>2);
my $f091  = file($dir, '09_anno.txt');
my $f09   = file($dir, '09_anno.tlf');

#annHclust(-out=>$f09, -in1=>$f081, -in2=>$f091);
my $f101 = file($dir, '10_exp.nwk');
=cut

sub extract_gtb_by_ids {
    my ($f00, $f_gtb, $f01) = @_;
    my $t0 = readTable(-in=>$f00, -header=>1);
    my $h = readGtb(-in=>$f_gtb, -opt=>2);

    open(FH, ">$f01");
    print FH join("\t", qw/id parent chr beg end strand locE locI locC loc5 loc3 phase source conf cat1 cat2 cat3 note/)."\n";
    for my $i (0..$t0->nofRow-1) {
        my ($id_g, $cat) = $t0->row($i);
        my $id = "$id_g.1";
        if(exists $h->{$id}) {
            my @ps = @{$h->{$id}};
            $ps[16] = $cat;
            print FH join("\t", @ps)."\n";
        } else {
            print "$id_g not found in [$f_gtb]\n";
        }
    }
    close FH;
}



