#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  genomestat.pl - generating genome statistics

=head1 SYNOPSIS
  
  genomestat.pl [-help] [-o organism]

  Options:
    -h (-help)   brief help message
    -o (--org)   organism

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use File::Path qw/make_path remove_tree/;
use Common;

my ($org, $k) = ('', 60);
my ($beg, $end, $chr) = (1, 1, 'chr1');
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "org|o=s"  => \$org,
  "kmer|k=i"  => \$k,
  "beg|b=i"  => \$beg,
  "end|e=i"  => \$end,
  "chr|c=s"  => \$chr,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$org;

my $dir = "/home/youngn/zhoup/Data/genome/$org/18_stat_k$k";
-d $dir || make_path($dir);
chdir $dir || die "cannot chdir to $dir\n";

##winslide.pl -i ../15.sizes -step 1 -size $K | seqdbret.pl -d ../11_genome.fa -o 01.fa
#seqtile.pl -i ../11_genome.fa -o 01.fa -step 1 -size $K

#echo -e "$CHR" | seqdbret.pl -d ../11_genome.fa | \
#  seqtile.pl -step 1 -size $K | \
#  bowtie2 -x $data/db/bowtie2/$org -f -U - \
#    --end-to-end --fast -k 100 -p 16 --reorder | \
#  sammapp.pl -mis 3 | tbl2bed.pl -o 15_mapp_$CHR.bed
#cat 15_mapp_*.bed > 15_mapp.bed
#bedGraphToBigWig 15_mapp.bed ../15.sizes 15_mapp.bw


#seqgc.pl -i 01.fa -o 11_gc.tbl
#tbl2bed.pl -i 11_gc.tbl -o 11_gc.bed
#bedGraphToBigWig 11_gc.bed ../15.sizes 11_gc.bw

#seq $BEG $END | xargs -i printf "%02d\\n" {} | \
#  parallel hybrid-ss-min -n DNA 00_split/part{} -o 17_deltag/part{}
#seq $BEG $END | xargs -i printf "17_deltag/part%02d\\n" {} | \
#  xargs -i rm {}.ct {}.run {}.37.ext {}.37.plot
#seq $BEG $END | xargs -i printf "part%02d\\n" {} | \
#  parallel -q sh -c "unafoldo.pl -i 17_deltag/{}.dG -s 00_split/{} | tbl2bed.pl -o 17_deltag/{}.bed"
#cat 17_deltag/part*.bed > 17_deltag.bed
#bedGraphToBigWig 17_deltag.bed ../15.sizes 17_deltag.bw

__END__
