#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  comp.batch.pl - batch genome comparison

=head1 SYNOPSIS
  
  comp.batch.pl [-help] 

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

for my $qry (@qrys) {
  my $tag = $qry;
  $tag =~ s/HM//i;
  my $dir = "$ENV{'misc3'}/$qry\_$tgt/23_blat";
  my $qsize = "$ENV{'genome'}/$qry/15.sizes";
  my $tsize = "$ENV{'genome'}/$tgt/15.sizes";
  chdir $dir || die "cannot chdir to $dir\n";
#  runCmd("comp.vcf.pl -q $qry -t $tgt");
  runCmd("idm.refine.pl -q $qry -t $tgt -i 31.9/idm -o 31.9/sv");
#  runCmd("comp.ortho.pl -q $qry -t $tgt");
#  print "qsub spada -N spada.$tag -v ORG=$qry -l qos=weightlessqos\n";
}

__END__

