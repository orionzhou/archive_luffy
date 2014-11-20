#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  comp.ucsc.pl - generate UCSC trackDb.txt

=head1 SYNOPSIS
  
  comp.ucsc.pl [-help] 

  Options:
    -h (--help)   brief help message

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin";
use Common;
use Data::Dumper;
use File::Path qw/make_path remove_tree/;
use File::Basename;
use List::Util qw/min max sum/;

my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
) or pod2usage(2);
pod2usage(1) if $help_flag;

my @qrys = qw/
  HM058 HM125 HM056 HM129 HM060
  HM095 HM185 HM034 HM004 HM050 
  HM023 HM010 HM022 HM324 HM340
/;
my $tgt = "HM101";

my $dir = "$ENV{'misc3'}/comp.ucsc";
-d $dir || make_path($dir);
chdir $dir || die "cannot chdir to $dir\n";

my $ft = "trackDb.txt";

my @tkeys = qw/track group type shortLabel longLabel bigDataUrl
  visibility priority color/;
my @tracks = (
  ['gap', 'map', 'bigBed 3', 'Gaps', 'Assembly Gaps', '16.gap.bb', 'squish', 1, '153,0,76'],
  ['gene', 'genes', 'bigBed 12', 'Gene', 'All Genes', '51.gtb.bb', 'pack', 1, '', 'itemRgb=on'],
  ['mappability_60mer', 'genomestat', 'bigWig 0 1', 'mapp_60mer', '60mer Mappability', '18_stat_k60/15_mapp.bw', 'full', 1, '64,64,64', 'maxHeightPixels=30:20:10']
);

my @cols = ('30,144,255', '60,179,113', '255,102,178', '64,64,64');
for my $i (0..$#qrys) {
  my $qry = $qrys[$i];
  my $col = $cols[$i%4];
  push @tracks, ["$qry\_blat", 'blat', 'bigBed 12', 
    "$qry\_blat", "$qry Blat",
    "$qry/31.9/gal.bb", "squish", $i+1, $col];
  push @tracks, ["$qry\_blatx", 'blatx', 'bigBed 6', 
    "$qry\_blatx", "$qry Blatx",
    "$qry/31.9/gax.bb", "hide", $i+1, $col];
  push @tracks, ["$qry\_blat_snp", 'blat_snp', 'bigBed 4', 
    "$qry\_blat_snp", "$qry Blat SNP",
    "$qry/31.9/snp.bb", "dense", $i+1, $col];
  
  my $diri = "$ENV{'misc3'}/$qry\_$tgt/23_blat/31.9";
  runCmd("rm -rf $qry");
  -d "$qry/31.9" || make_path("$qry/31.9");
  runCmd("cp $diri/gal.bb $qry/31.9/");
  runCmd("cp $diri/gax.bb $qry/31.9/");
  runCmd("cp $diri/snp.bb $qry/31.9/");
}

open(my $fht, ">$ft") or die "cannot write $ft\n";
for (@tracks) {
  my @tvals = @$_;
  @tvals >= @tkeys || die Dumper(@tvals)." < ".@tkeys." values\n";
  for my $i (0..$#tkeys) {
    $tvals[$i] eq "" && next;
    print $fht join(" ", $tkeys[$i], $tvals[$i])."\n";
  }
  if(@tvals == @tkeys) {
    print $fht "\n";
    next;
  }
  for my $j (scalar(@tkeys)..$#tvals) {
    my ($tkey, $tval) = split("=", $tvals[$j]);
    print $fht join(" ", $tkey, $tval)."\n";
  }
  print $fht "\n";
}
close $fht;

-s "all.tgz" && runCmd("rm all.tgz");
runCmd("tar czf all.tgz *");
__END__

