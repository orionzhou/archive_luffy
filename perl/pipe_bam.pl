#!/usr/bin/perl
use strict;
use lib ($ENV{"SCRIPT_HOME_PERL"});
use Getopt::Long;
use Data::Dumper;
use File::Basename;
use File::Path qw/make_path remove_tree/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;

use InitPath;
use Common;
use Bam;

my ($program, $lb) = ('', '');
GetOptions('lib|l=s'=>\$lb, 'program|p=s'=>\$program);

my $f_ref = "$DIR_genome/Mtruncatula_4.0/01_refseq.fa";
my $f_lst = "$DIR_misc3/ncgr_fastq/11_library.tbl";
my $dir = "$DIR_misc3/hapmap_mt40/15_pipe_bam";
my $dir_fq = "$DIR_misc3/hapmap_mt40/11_pipe_bwa/06_pos_sorted";

pipe_bam_run($dir, $f_lst, $dir_fq, $lb) if $program eq "run";
sub read_lib_list {
    my ($f_lst) = @_;
    my $t = readTable(-in=>$f_lst, -header=>1);
    my $h;
    for my $i (0..$t->nofRow-1) {
        my ($sm, $lb, $idxs, $rl, $pi, $pl, $pd) = $t->row($i);
        my @idxs = split(",", $idxs);
        $h->{$lb} = [$sm, \@idxs, $rl, $pi, $pl, $pd];
    }
    return $h;
}
sub pipe_bam_run {
    my ($dir, $f_lst, $dir_fq, $lb) = @_;
    my $h = read_lib_list($f_lst);
   
    die "no library named '$lb'\n" unless exists $h->{$lb};
    my ($sm, $idxs, $rl, $pi, $pl) = @{$h->{$lb}};
    print join(", ", @$idxs)."\n";
    my @fis;
    for my $idx (@$idxs) {
        my $fi = "$dir_fq/$idx.bam";
        die "$idx [$fi] is not there\n" unless -s $fi;
        push @fis, $fi;
    }

    my $cmd;
    my $d05 = "$dir/05_dup_marked";
    make_path($d05) unless -d $d05;
    my $input_str = join(" ", map {"INPUT=$_"} @fis);
    my ($f05, $f05b) = ("$d05/$lb.bam", "$d05/$lb.txt");
    $cmd = "java -Xmx12g -jar $picard/MarkDuplicates.jar \\
        VALIDATION_STRINGENCY=LENIENT TMP_DIR=$DIR_tmp \\
        $input_str OUTPUT=$f05 METRICS_FILE=$f05b";
    runCmd($cmd);
    runCmd("samtools index $f05");
    
    my $d07 = "$dir/07_realigned";
    make_path($d07) unless -d $d07;
    my ($f07, $f07b) = ("$d07/$lb.bam", "$d07/$lb.intervals");
    $cmd = "java -Xmx12g -Djava.io.tmpdir=$DIR_tmp -jar $gatk/GenomeAnalysisTK.jar \\
        -T RealignerTargetCreator -I $f05 -R $f_ref -o $f07b";
    runCmd($cmd);
    $cmd = "java -Xmx12g -Djava.io.tmpdir=$DIR_tmp -jar $gatk/GenomeAnalysisTK.jar \\
        -T IndelRealigner -I $f05 -R $f_ref -targetIntervals $f07b -o $f07 \\
        -LOD 0.4 --maxReadsForRealignment 20000 --maxReadsInMemory 200000";

    my $d12 = "$dir/12_stat";
    make_path($d12) unless -d $d12;
    my $f12 = "$d12/$sm";
    runCmd("bamISD -i $f07 -o $f12");
}

sub sort_readname {
    my ($dir, $lb) = @_;
    my $d13 = "$dir/13_dup_removed";
    my $d21 = "$dir/21_rn_sorted";
    my $cmd = "java -Xmx12g -jar $picard/SortSam.jar \\
        VALIDATION_STRINGENCY=LENIENT TMP_DIR=$DIR_tmp \\
        SORT_ORDER=queryname INPUT=$d13/$lb.bam OUTPUT=$d21/$lb.bam";
    runCmd($cmd);
}
  
stat_dup($dir) if $program eq "stat_dup";
stat_isd($dir) if $program eq "stat_isd";
sub stat_isd {
    my ($dir) = @_;
    my $f01 = "$dir/01_sample.tbl";
    my $d12 = "$dir/12_stat";
    my $t = readTable(-in=>$f01, -header=>1);
    my $to1 = Data::Table->new([], [qw/rg total unmapped unpaired paired unpaired_dup unpaired_uniq paired_dup paired_uniq paired_proper/]);
    my $to2 = Data::Table->new([], [qw/rg is cnt/]);
    for my $i (0..$t->nofRow-1) {
        my ($idx, $sm, $lbs, $rgs, $pi, $pl) = $t->row($i);
        my $f12a = "$d12/$sm.tbl";
        my $ta = readTable(-in=>$f12a, -header=>1);
        for my $j (0..$ta->nofRow-1) {
            $to1->addRow($ta->rowRef($j));
        }

        my $f12b = "$d12/$sm\_isd.tbl";
        my $tb = readTable(-in=>$f12b, -header=>1);
        for my $j (0..$tb->nofRow-1) {
            $to2->addRow($tb->rowRef($j));
        }
    }
    open(FH1, ">$dir/13_stat.tbl");
    open(FH2, ">$dir/13_stat_isd.tbl");
    print FH1 $to1->tsv(1);
    print FH2 $to2->tsv(1);
    close FH1;
    close FH2;
} 
sub stat_dup {
    my ($dir) = @_;
    my $f01 = "$dir/01_sample.tbl";
    my $d03 = "$dir/03_pos_sorted";
    my $d05 = "$dir/05_dup_marked";
    my $f11 = "$dir/11_stat_dup.tbl";

    my $t = readTable(-in=>$f01, -header=>1);
    open(FHO, ">$f11");
    print FHO join("\t", qw/sm lb unpaired_reads_examined read_pairs_examined unmapped_reads unpaired_read_duplicates read_pair_duplicates read_pair_optical_duplicates percent_duplicates estimated_library_size/)."\n";
    for my $i (0..$t->nofRow-1) {
        my ($idx, $sm, $lbs, $rgs, $pi, $pl) = $t->row($i);
        my ($f05, $f05b) = ("$d05/$sm.bam", "$d05/$sm.txt");
        unless( -s $f05b ) {
            print FHO join("\t", $sm, ("") x 9)."\n";
            next;
        }
        open(FHI, "<$f05b");
        while(<FHI>) {
            next unless /^HM/;
            chomp;
            my @ps = split("\t", $_, -1);
            print FHO join("\t", $sm, @ps)."\n";
        }
        close FHI;
    }
    close FHO;
}




