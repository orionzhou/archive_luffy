#!/usr/bin/perl
use strict;
use lib ($ENV{"SCRIPT_HOME_PERL"});
use Getopt::Long;
use Data::Dumper;
use File::Basename;
use File::Path qw/make_path remove_tree/;

use InitPath;
use Common;
use Bam;

my ($opt, $rn, $lb, $sm, $beg, $end) = (('') x 4, 0, -1);
GetOptions(
    'opt|p=s' => \$opt,
    'rn|r=s'  => \$rn,
    'lb|l=s'  => \$lb,
    'sm|s=s'  => \$sm,
    'beg|b=i' => \$beg, 'end|e=i' => \$end, 
);

my $f_ref = "$DIR_genome/Mtruncatula_4.0/01_refseq.fa";
my $f_bwa = "$DIR_db/bwa/mt_40";
my $f_rn = "$DIR_misc3/ncgr_fastq/04_fastq_stats.tbl";
my $f_lb = "$DIR_misc3/ncgr_fastq/11_library.tbl";
my $f_sm = "$DIR_misc3/ncgr_fastq/21_sample.tbl";
my $dir = "$DIR_misc3/hapmap_mt40/11_pipe_mapping";
my $d03 = "$dir/03_bwa";
my $d06 = "$dir/06_pos_sorted";
my $f09 = "$dir/09_status_run.tbl";
my $d11 = "$dir/11_dedup";
my $f19 = "$dir/19_status_lib.tbl";
my $d21 = "$dir/21_realigned";
my $f29 = "$dir/29_status_sample.tbl";
my $d31 = "$dir/31_stat";

if($opt eq "run") {
    pipe_run($dir, $f_rn, $f_bwa, $rn, $beg, $end);
} elsif($opt eq "lib") {
    pipe_lib($dir, $f_lb, $d06, $lb);
} elsif($opt eq "sample") {
    pipe_sample($dir, $f_sm, $d11, $f_ref, $sm);
} elsif($opt eq "update") {
    status_update($dir, $f_rn, $f_lb, $f_sm, $f09, $f19, $f29);
} else {
    die "unknown option: $opt\n";
}

sub pipe_run {
    my ($dir, $f_rn, $f_bwa, $rni, $beg, $end) = @_;
    my $d03 = "$dir/03_bwa";
    make_path($d03) unless -d $d03;
    my $d06 = "$dir/06_pos_sorted";
    make_path($d06) unless -d $d06;
    
    my $dir_abs = dirname($f_rn);
    my $t = readTable(-in=>$f_rn, -header=>1);
    for my $i (0..$t->nofRow-1) {
        my ($idx, $sm, $lb, $rn, $pl, $rl, $pi, $dir_rel, $n_seq, $encoding) = $t->row($i);
        next if $rn ne $rni;
#        next if $idx < $beg || $idx > $end;
        my $f1a = "$dir_abs/$dir_rel/$rn.1.fq.gz";
        my $f1b = "$dir_abs/$dir_rel/$rn.2.fq.gz";
        die "$f1a is not there\n" unless -s $f1a;
        die "$f1b is not there\n" unless -s $f1b;

        my $f3a = "$d03/$rn.1.sai";
        my $f3b = "$d03/$rn.2.sai";
        my $f3  = "$d03/$rn.bam";
        my $tag = ($encoding < 1.8 & $encoding >= 1.3) ? "-I" : "";
        runCmd("bwa aln -t 4 -n 0.01 $tag $f_bwa $f1a > $f3a", 1);
        runCmd("bwa aln -t 4 -n 0.01 $tag $f_bwa $f1b > $f3b", 1);

        my $tag_is = ($pi == 6500) ? "-a 10000" : ($pi == 3000) ? "-a 6000" : "";
        my $f6 = "$d06/$rn.bam";
        runCmd("bwa sampe $f_bwa $tag_is \\
            -r '\@RG\\tID:$rn\\tSM:$sm\\tLB:$lb\\tPL:$pl\\tPU:lane' \\
            $f3a $f3b $f1a $f1b | samtools view -Sb - > $f3", 1);
        runCmd("java -Xmx8g -jar $picard/FixMateInformation.jar \\
            TMP_DIR=$DIR_tmp VALIDATION_STRINGENCY=LENIENT \\
            INPUT=$f3 OUTPUT=$f6 SORT_ORDER=coordinate", 1);
    }
}
sub pipe_lib {
    my ($dir, $f_lb, $dirI, $lb) = @_;

    my $t = readTable(-in=>$f_lb, -header=>1);
    my $h;
    for my $i (0..$t->nofRow-1) {
        my ($sm, $lb, $rns, $idxs) = $t->row($i);
        my @runs = split(",", $rns);
        $h->{$lb} = \@runs;
    }
    die "no library named '$lb'\n" unless exists $h->{$lb};
   
    my @rns = @{$h->{$lb}};
    print join("\n", @rns)."\n";
    my @fis;
    for my $rn (@rns) {
        my $fi = "$dirI/$rn.bam";
        die "$rn [$fi] is not there\n" unless -s $fi;
        push @fis, $fi;
    }

    my $d11 = "$dir/11_dedup";
    make_path($d11) unless -d $d11;
    my $input_str = join(" ", map {"INPUT=$_"} @fis);
    my ($f11, $f11b) = ("$d11/$lb.bam", "$d11/$lb.dup.txt");
    runCmd("java -Xmx10g -jar $picard/MarkDuplicates.jar \\
        VALIDATION_STRINGENCY=LENIENT TMP_DIR=$DIR_tmp \\
        REMOVE_DUPLICATES=true \\
        $input_str OUTPUT=$f11 METRICS_FILE=$f11b", 1);
    runCmd("samtools index $f11", 1);
    runCmd("bamStat -i $d11/$lb.bam -o $d11/$lb", 1);
}
sub pipe_sample {
    my ($dir, $f_sm, $dirI, $f_ref, $sm) = @_;

    my $t = readTable(-in=>$f_sm, -header=>1);
    my $h;
    for my $i (0..$t->nofRow-1) {
        my ($sm, $lbs) = $t->row($i);
        my @lbs = split(",", $lbs);
        $h->{$sm} = \@lbs;
    }
    die "no sample named '$sm'\n" unless exists $h->{$sm};
   
    my @lbs = @{$h->{$sm}};
    print join("\n", @lbs)."\n";
    my @fis;
    for my $lb (@lbs) {
        my $fi = "$dirI/$lb.bam";
        die "$lb [$fi] is not there\n" unless -s $fi;
        push @fis, $fi;
    }

    my $d21 = "$dir/21_realigned";
    make_path($d21) unless -d $d21;
    my $input_str = join(" ", map {"-I $_"} @fis);
    my ($f21, $f21b) = ("$d21/$sm.bam", "$d21/$sm.intervals");
    runCmd("java -Xmx12g -Djava.io.tmpdir=$DIR_tmp -jar $gatk/GenomeAnalysisTK.jar \\
        -T RealignerTargetCreator -R $f_ref -o $f21b", 1);
    runCmd("java -Xmx12g -Djava.io.tmpdir=$DIR_tmp -jar $gatk/GenomeAnalysisTK.jar \\
        -T IndelRealigner $input_str \\
        -R $f_ref -targetIntervals $f21b -o $f21 \\
        -LOD 0.4 --maxReadsForRealignment 20000 --maxReadsInMemory 200000", 1);

    my $d31 = "$dir/31_stat";
    make_path($d31) unless -d $d31;
    my $f31 = "$d31/$sm";
    runCmd("bamStat -i $f21 -o $f31", 1);
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

sub status_update {
    my ($dir, $f_run, $f_lib, $f_sam, $fo1, $fo2, $fo3) = @_;
    
    my $d03 = "$dir/03_bwa";
    my $d06 = "$dir/06_pos_sorted";
    open(FH1, ">$fo1") or die "Cannot open $fo1 for writing\n";
    print FH1 join("\t", qw/idx sm lb rn 03_bwa 06_pos_sorted/)."\n";
    my $tr = readTable(-in=>$f_run, -header=>1);
    for my $i (0..$tr->nofRow-1) {
        my ($idx, $sm, $lb, $rn, $pl, $rl, $pi, $dir_rel, $n_seq, $encoding) = $tr->row($i);
        my $tag03 = check_bam("$d03/$rn.bam");
        my $tag06 = check_bam("$d06/$rn.bam");
        $tag06 = 0 if $tag03 == 0;
        print FH1 join("\t", $idx, $sm, $lb, $rn, $tag03, $tag06)."\n";
        print "  checking run $rn\r";
    }
    print "\n";
    close FH1;
   
    my $d11 = "$dir/11_dup_marked";
    open(FH2, ">$fo2") or die "Cannot open $fo2 for writing\n";
    print FH2 join("\t", qw/sm lb 11_dup_marked rns/)."\n";
    my $tl = readTable(-in=>$f_lib, -header=>1);
    for my $i (0..$tl->nofRow-1) {
        my ($sm, $lb, $rns, $idxs) = $tl->row($i);
        my $tag11 = check_bam("$d11/$lb.bam");
        print FH2 join("\t", $sm, $lb, $tag11, $rns, $idxs)."\n";
        print "  checking lib $lb\r";
    }
    print "\n";
    close FH2;

    my $d21 = "$dir/21_realigned";
    open(FH3, ">$fo3") or die "Cannot open $fo3 for writing\n";
    print FH3 join("\t", qw/sm 21_realigned lbs/)."\n";
    my $ts = readTable(-in=>$f_sam, -header=>1);
    for my $i (0..$ts->nofRow-1) {
        my ($sm, $lbs) = $ts->row($i);
        my $tag21 = check_bam("$d21/$sm.bam");
        print FH3 join("\t", $sm, $tag21, $lbs)."\n";
        print "  checking sample $sm\r";
    }
    print "\n";
    close FH3;
}


#my $cmd1 = "novoalign -d $DIR_db/novoalign/mt_35 -f $f1a -f2 $f1b -r E 10 -t 99 -o SAM | samtools view -Sb - > $f1";

