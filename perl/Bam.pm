package Bam; 
use strict;
use Init;
use Common;
use List::Util qw/min max sum/; use POSIX qw/ceil floor/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw/$picard $gatk $svtoolkit
    bam_sort check_bam/;
@EXPORT_OK = qw//;

our $picard = $ENV{"src"}."/picard/dist";
our $gatk = $ENV{"src"}."/gatk/dist";
our $svtoolkit = $ENV{"src"}."/svtoolkit";
sub bam_sort {
    my ($fi, $fo_pre) = @_;
    my $cmd = "java -Xmx4g -jar $picard/SortSam.jar \\
        VALIDATION_STRINGENCY=LENIENT TMP_DIR=$DIR_Tmp \\
        INPUT=$fi OUTPUT=$fo_pre";
    runCmd($cmd);
}
sub check_bam {
    my ($fi) = @_;
    my $tag = 1;
    $tag = 0 if ! -s $fi;
    my $cmd = "samtools view -H $fi";
    open(JJ, $cmd." 2>&1 |") || die "Failed: $! in \n$cmd\n";
    my ($tag_sq, $tag_rg, $tag_pg) = (0, 0, 0);
    while ( <JJ> ){
        chomp;
        $tag = 0 if /EOF marker is absent/;
        $tag = 0 if /fail to open/;
        $tag_sq = 1 if /^\@SQ/;
        $tag_rg = 1 if /^\@RG/;
        $tag_pg = 1 if /^\@PG/;
    }
    $tag = 0 if $tag_sq+$tag_rg+$tag_pg < 3;
    return $tag;
}
  



1;
__END__
  
