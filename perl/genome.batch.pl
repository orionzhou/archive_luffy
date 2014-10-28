#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  genome.batch.pl - batch processing genome seqs

=head1 SYNOPSIS
  
  genome.batch.pl [-help] 

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

my @orgs = qw/
  HM004 HM010 HM018 HM022 HM034 
  HM050 HM056 HM058 HM060 HM095 
  HM125 HM129 HM185 HM324 HM340
  HM056.AP HM340.AP/;

for my $org (@orgs) {
#  print "qsub rm -N rm.$org -v ORG=$org -l qos=weightlessqos\n";
#  runCmd("mt.nbs.pl -g $org");
#  runCmd("gff.rename.pl -i \$genome/$org/51.gff -m \$genome/$org/raw.fix.fas.map -o \$misc2/coge/$org.gff");
#  print "qsub spada -N spada.$org -v ORG=$org -l qos=weightlessqos\n";
#  runCmd("mt.augus.pl -g $org");
#  runCmd("mt.anno.pl -g $org");
#  runCmd("genome.fas.pl -g $org");
  runCmd("genome.db.pl -g $org");
}

__END__

