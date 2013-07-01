#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gtb_compare.pl - make comparison of two Gtb files

=head1 SYNOPSIS
  
  gtb_compare.pl [-help] [-query query-Gtb] [-target target-Gtb] [-comp comparison-file] [-out merged-Gtb]

  Options:
      -help    brief help message
      -query   query Gtb file
      -target  target Gtb file
      -comp    comparison result file
      -out     output merged Gtb file

=head1 DESCRIPTION

  This program compares gene models in two Gtb files, gives a summary and write a Gtb with non-redundant models

=cut
  
#### END of POD documentation.
#-----------------------------------------------------------------------------
use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;
use Gtb;
use CompareModel;

my ($fq, $ft, $fc, $fo) = ('') x 4;
my $help_flag;

#----------------------------------- MAIN -----------------------------------#
GetOptions(
    "help|h"  => \$help_flag,
    "query|qry|q=s"  => \$fq,
    "target|tgt|t=s" => \$ft,
    "comp|c=s" => \$fc,
    "out|o=s"  => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fq || !$ft || !$ft;

my $f_tmp = "$fc.tmp";
compare_models($fq, $ft, $f_tmp);
model_eval($f_tmp, $ft, $fc);
if($fo) {
    my ($sn_nt, $sp_nt, $sn_ex, $sp_ex) = get_sn_sp($fc);
    printf "sn_nt=%.03f sp_nt=%.03f sn_ex=%.03f sp_ex=%.03f\n", $sn_nt, $sp_nt, $sn_ex, $sp_ex;
    gtb_rmdup($fq, $ft, $fc, $fo);
}
system("rm $f_tmp");

sub gtb_rmdup {
    my ($f_qry, $f_tgt, $f_cmp, $fo) = @_;
    my $h;
    my $tc = readTable(-in=>$f_cmp, -header=>1);
    for my $i (0..$tc->nofRow-1) {
        my ($idQ, $tag) = $tc->row($i);
        $h->{$idQ} = 1 if $idQ ne "" && $tag != 1;
    }

    my $tt = readTable(-in=>$f_tgt, -header=>1);
    my $tq = readTable(-in=>$f_qry, -header=>1);
    $tt->addCol([("") x $tt->nofRow], "seq") if $tt->nofCol == 18;
    $tq->addCol([("") x $tq->nofRow], "seq") if $tq->nofCol == 18;
    for my $i (0..$tq->nofRow-1) {
        my ($idQ) = $tq->row($i);
        if(exists $h->{$idQ}) {
            $tt->addRow( $tq->rowRef($i) );
        }
    }

    open(FHO, ">$fo") or die "cannot write to $fo\n";
    print FHO $tt->tsv(1);
    close FHO;
}

__END__
