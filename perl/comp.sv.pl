#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  comp.sv.pl - 

=head1 SYNOPSIS
  
  comp.sv.pl [-help] [-qry qry-genome] [-tgt tgt-genome]

  Options:
    -h (--help)   brief help message

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use File::Path qw/make_path remove_tree/;
use File::Basename;
use Tabix;
use Bio::DB::Fasta;
use Common;
use Location;
use Gtb;
use Gal;
use List::Util qw/min max sum/;

my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
) or pod2usage(2);
pod2usage(1) if $help_flag;

my @qrys = qw/
  HM058 HM125 HM056 HM129 HM060
  HM095 HM185 HM034 HM004 HM050 
  HM023 HM010 HM022 HM324 HM340
/;
@qrys = qw/HM004/;
my $tgt = "HM101";

my $dir = "$ENV{'misc3'}/gene_sv";
-d $dir || make_path($dir);
chdir $dir || die "cannot chdir to $dir\n";

my $tdir = "$ENV{'genome'}/$tgt";
runCmd("awk 'BEGIN{OFS=\"\\t\"} {if(\$6==\"cds\") print \$1, \$2-1, \$3}' \\
  $tdir/51.tbl | sortBed -i stdin | uniq > $tgt.bed");

for my $qry (@qrys) {
my $qdir = "$ENV{'genome'}/$qry";
runCmd("awk 'BEGIN{OFS=\"\\t\"} {if(\$6==\"cds\") print \$1, \$2-1, \$3}' \\
  $qdir/51.tbl | sortBed -i stdin | uniq > $qry.bed");

my $cdir = "$ENV{'misc3'}/$qry\_$tgt/23_blat/31.9";
runCmd("intersectBed -wao -a $qry.bed -b $cdir/sv.ins.bed | \\
  sortBed -i stdin > $qry.ins.bed");
runCmd("intersectBed -wao -a $qry.bed -b $cdir/sv.gan.bed | \\
  sortBed -i stdin > $qry.gan.bed");
runCmd("intersectBed -wao -a $tgt.bed -b $cdir/sv.del.bed | \\
  sortBed -i stdin > $qry.del.bed");
runCmd("intersectBed -wao -a $tgt.bed -b $cdir/sv.los.bed | \\
  sortBed -i stdin > $qry.los.bed");
runCmd("intersectBed -wao -a $qry.bed -b $cdir/sv.tlc.gan.bed | \\
  sortBed -i stdin > $qry.tlc.gan.bed");
runCmd("intersectBed -wao -a $tgt.bed -b $cdir/sv.tlc.los.bed | \\
  sortBed -i stdin > $qry.tlc.los.bed");

}

__END__

