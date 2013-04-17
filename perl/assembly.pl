#!/usr/bin/perl -w
use strict;
use lib ($ENV{"SCRIPT_HOME_PERL"});
use Common;
use Seq;
use Align;
use Blast;
use Eutils;
use Data::Dumper;
use File::Path qw/make_path remove_tree/;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;

=begin
my $dir = "/home/youngn/zhoup/Data/misc3/hm056/mummer";
my $f_seq = "$dir/../hm056.fa";
my $fl = "$dir/../hm056_seqlen.tbl";
#$src/MUMmer3.23/nucmer -p chr1-4 ../mt4_chr1-4.fa ../HM056.ALLPATHSLG.Jan2013.fasta
#$src/MUMmer3.23/nucmer -p chr5-8 ../mt4_chr5-8.fa ../HM056.ALLPATHSLG.Jan2013.fasta
#$src/MUMmer3.23/show-coords chr1-4.delta > chr1-4.coords
#$src/MUMmer3.23/show-coords chr5-8.delta > chr5-8.coords
#cat chr1-4.coords chr5-8.coords > 01_mummer.coords
my $f01 = "$dir/01_mummer.coords";
my $f02 = "$dir/02_coords.tbl";
#mummer_coords2tbl($f01, $f02);
my $f05 = "$dir/05_tiled.tbl";
#mummer_tiling($f02, $f05);
=end
=cut

my $acc = "hm056";
#my $acc = "hm340";
my $dir = "/home/youngn/zhoup/Data/misc3/$acc";
my $f_seq = "$dir/01_assembly.fa";
my $f_len = "$dir/11_seqlen.tbl";
my $f_gap = "$dir/12_gaploc.tbl";
#print "seqlen.pl -out $f_len $f_seq\n";
#print "seqgap.pl -out $f_gap -min 100 $f_seq\n";

my $d21 = "$dir/21_blastn";
make_path($d21) unless -d $d21;
my $f21_02 = "$d21/02_raw.tbl";
my $f21_05 = "$d21/05_tiled.tbl";
runCmd("blastn -db \$data/db/blast/mt4.0 -outfmt 6 -query $f_seq -out $f21_02\n", 1);
runCmd("blastTiling -i $f21_02 -o $f21_05");

# run R script assembly.R
my $f21_21 = "$d21/21_scaffold_status.tbl";
my $f21_22 = "$d21/22_unmapped.tbl";

my $f21_31 = "$d21/31_unmapped.fa";
#print "seqextract.pl -out $f21_31 -id $f21_22 $f_seq\n";
my $f21_32 = "$d21/32_raw.tbl";
my $f21_35 = "$d21/35_tiled.tbl";
#print "blastn -db /project/db/blast/current/nt -outfmt 6 -query $f21_31 -out $f21_32\n";
#print "blastTiling -i $f21_32 -o $f21_35\n";
my $f21_36 = "$d21/36_annotated.tbl";
#annotate_nr_hits($f21_35, $f21_36);

# run R script assembly.R
my $f21_41 = "$d21/41_scaffold_status.tbl";
my $f21_42 = "$d21/42_unknown.tbl";

my $f21_51 = "$d21/51_unknown.fa";
#print "seqextract.pl -out $f21_51 -id $f21_42 $f_seq\n";

sub annotate_nr_hits {
    my ($fi, $fo) = @_;
    my $t = readTable(-in=>$fi, -header=>1);
    my @gis;
    for my $i (0..$t->nofRow-1) {
        my $id = $t->elm($i, "hId");
        if($id =~ /gi\|(\d+)\|/) {
            push @gis, $1;
        } else {
            die "unknown hId: $id\n";
        }
    }
    my ($cats, $species) = get_gi_taxonomy(@gis);
    $t->addCol($cats, "category");
    $t->addCol($species, "species");
    open(FH, ">$fo") or die "cannot write to $fo\n";
    print FH $t->tsv(1);
    close FH;
}


__END__

