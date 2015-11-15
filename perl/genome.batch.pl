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
use Medicago;

my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
) or pod2usage(2);
pod2usage(1) if $help_flag;

print "#qsub rnaseq\n";

for my $org ($tname, @qnames) {
  my $dir = "$ENV{'genome'}/$org";
  chdir $dir || die "cannot chdir to $dir\n";
#  runCmd("genome.fas.pl -g $org");
#  runCmd("genome.db.pl -g $org");
#  print "qsub rm -N rm.$org -v ORG=$org -l qos=weightlessqos\n";
#  runCmd("mt.augus.pl -g $org");
  runCmd("mt.anno.pl -g $org");
#  print "qsub mtanno -N mtanno.$org -v ORG=$org -l qos=weightlessqos");
#  runCmd("gff.rename.pl -i \$genome/$org/51.gff -m \$genome/$org/raw.fix.fas.map -o \$misc2/coge/$org.gff");
}

__END__

