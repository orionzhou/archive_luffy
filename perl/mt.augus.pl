#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  mt.augus.pl - use augustus to annotate an Medicago genome

=head1 SYNOPSIS
  
  mt.augus.pl [-help] [-org organism]

  Options:
    -h (--help)   brief help message
    -g (--org)    genmome ID (organism) to process

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin";
use Common;
use Location;
use Data::Dumper;
use File::Path qw/make_path remove_tree/;
use File::Basename;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;

my ($org, $ncpu) = ('', 16);
my $help_flag;
#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "org|g=s"  => \$org,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$org;

my $dir = "$ENV{'genome'}/$org/augustus";
-d $dir || make_path($dir);
chdir $dir || die "cannot chdir $dir\n";

my $fg = "../11_genome.fas";
-s $fg || die "$fg is not there\n";

get_hints_ortholog();
get_hints_rnaseq();
run_aug();
postprocess_aug();
pipe_pfam();

sub get_hints_ortholog {
  my $f_gtb = "$ENV{'genome'}/HM101/51.gtb";
  my $f_ref = "$ENV{'genome'}/HM101/11_genome.fas";
  my $f_gax = "$ENV{'misc3'}/$org\_HM101/23_blat/31.9/gax.gz";
  my $f_snp = "$ENV{'misc3'}/$org\_HM101/23_blat/31.9/snp.gz";
  runCmd("liftover.gtb2hint.pl -i $f_gtb -r $f_ref -x $f_gax -s $f_snp -o 05.hm101.gff");
}
sub get_hints_rnaseq {
  my $dirin = "$ENV{'misc2'}/rnaseq/mt";
  
  my $t = readTable(-in => "$dirin/21.tbl", -header => 1);
  my $h;
  for my $i (0..$t->lastRow) {
    my ($org, $orgr) = $t->row($i);
    !exists $h->{$org} || die "$org has >=2 mappings\n";
    $h->{$org} = $orgr;
  }
  exists $h->{$org} || die "no RNA-seq for $org\n";
  my $orgr = $h->{$org};
  
  my $bam_in = "$dirin/22_tophat/$orgr\_$org/accepted_hits.bam";
  -s $bam_in || die "$bam_in is not there\n";
  runCmd("bamtools filter -isPrimaryAlignment 1 -isMapped 1 -isMateMapped 1 -in $bam_in -out 11.f.bam");
  runCmd("samtools view -H 11.f.bam > 12.header.txt");

  runCmd("samtools sort 11.f.bam 14.sf");
  runCmd("bam2hints --intronsonly --in=14.sf.bam --out=15.rnaseq.gff");
#  runCmd("bam2wig $bam_in | \$soft/augustus/scripts/wig2hints.pl > 15.ep.gff");

  runCmd("rm 11.f.bam 12.header.txt 14.sf.bam");
}
sub run_aug {
  -s "05.hm101.gff" || die "no 05.hm101.gff";
  -s "15.rnaseq.gff" || die "no 15.rnaseq.gff";
  runCmd("cat 05.hm101.gff 15.rnaseq.gff > 20.hints.gff");
  runCmd("sort -k1,1 -k4,4n 20.hints.gff -o 20.hints.gff");

  my $dig = getDigits($ncpu);
  runCmd("ln -sf $fg 21.fas");
  runCmd("pyfasta split -n $ncpu 21.fas");
  my $cfg_aug = "$ENV{'AUGUSTUS_CONFIG_PATH'}/extrinsic/extrinsic.mt.cfg";
  my $end = $ncpu - 1;
  runCmd("seq 0 $end | xargs -i printf \"%0".$dig."d\\n\" {} | \\
    parallel -j $ncpu augustus --species=medicago \\
    --extrinsicCfgFile=$cfg_aug \\
    --alternatives-from-evidence=true \\
    --allow_hinted_splicesites=atac \\
    --introns=on --genemodel=partial \\
    --strand=both --gff3=on \\
    --hintsfile=20.hints.gff \\
    --outfile=21.{}.gff 21.{}.fas");
  runCmd("cat 21.*.gff | join_aug_pred.pl > 21.gff");
  runCmd("rm 21.fas.* 21.*.fas 21.*.gff");
}
sub postprocess_aug {
  runCmd("gff.augus.pl -i 21.gff -o - | gff2gtb.pl -i - -o 22.gtb");
  runCmd("gtb.rmutr.pl -i 22.gtb -o 23.gtb");
  runCmd("gtb.dedup.pl -i 23.gtb -o 25.dedup.gtb");
  runCmd("ln -sf 25.dedup.gtb 31.gtb");
  runCmd("gtb2gff.pl -i 31.gtb -o 31.gff");
  runCmd("gtb2fas.pl -i 31.gtb -d $fg -o 31.fas");
}
sub pipe_pfam {  
##runCmd("interproscan.sh -appl PfamA -dp -i 31.fas -f tsv -o 33.pfam.tsv");
##runCmd("pfam2tbl.pl -i 33.pfam.tsv -o 34.pfam.tbl -e 1 -l 10");
#  -s "33.txt" && runCmd("cp 33.txt 34.tbl.1.txt");
  runCmd("pfam.scan.pl -i 31.fas -o 34.tbl");
  runCmd("gtb.addpfam.pl -i 31.gtb -p 34.tbl -o 41.gtb"); 
}



__END__

