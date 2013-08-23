#!/usr/bin/perl
use strict;
use FindBin;
use lib $FindBin::Bin;
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

my $f_ref = "$DIR_genome/Mtruncatula_4.0/11_genome.fa";
my $f_bwa = "$DIR_db/bwa/mt_40";
my $f_rn = "$DIR_misc3/ncgr_fastq/04_fastq_stats.tbl";
my $f_lb = "$DIR_misc3/ncgr_fastq/11_library.tbl";
my $f_sm = "$DIR_misc3/ncgr_fastq/21_sample.tbl";
my $dir = "$DIR_misc3/hapmap_mt40/11_pipe_mapping";

my $d03 = "$dir/03_aln";
my $d04 = "$dir/04_fixmate_sort";

my $d12 = "$dir/12_dedup";
my $d19 = "$dir/19_remap_dedup";

my $d21 = "$dir/21_realigned";

my $d31 = "$dir/31_bcf_raw";

my $f51 = "$dir/51_status_run.tbl";
my $f52 = "$dir/52_status_lib.tbl";
my $f53 = "$dir/53_status_sample.tbl";

my $f61 = "$dir/61_stat_raw.tbl";
my $f62 = "$dir/62_isd.tbl";

if($opt eq "run") {
    pipe_run($dir, $f_rn, $f_bwa, $rn, $beg, $end);
} elsif($opt eq "rerun") {
    pipe_rerun($dir, $f_rn, $f_bwa, $rn);
} elsif($opt eq "lib") {
    pipe_lib($dir, $f_lb, $lb);
} elsif($opt eq "sample") {
    pipe_sample($dir, $f_sm, $f_ref, $sm);
} elsif($opt eq "vnt") {
    pipe_vnt($dir, $f_sm, $f_ref, $sm);
} elsif($opt eq "update") {
    status_update($dir, $f_rn, $f_lb, $f_sm, $f51, $f52, $f53);
} elsif($opt eq "stat") {
    stats($dir, $f_sm, $f61, $f62);
} else {
    die "unknown option: $opt\n";
}

sub pipe_run {
    my ($dir, $f_rn, $f_bwa, $rni, $beg, $end) = @_;
    my $d03 = "$dir/03_aln";
    make_path($d03) unless -d $d03;
    my $d04 = "$dir/04_fixmate_sort";
    make_path($d04) unless -d $d04;
    
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

        my $tag = ($encoding < 1.8 & $encoding >= 1.3) ? "-I" : "";
        runCmd("bwa aln -t 4 -n 0.01 $tag $f_bwa $f1a > $d03/$rn.1.sai", 1);
        runCmd("bwa aln -t 4 -n 0.01 $tag $f_bwa $f1b > $d03/$rn.2.sai", 1);

        my $tag_is = ($pi == 6500) ? "-a 10000" : ($pi == 3000) ? "-a 6000" : "";
        runCmd("bwa sampe $f_bwa \\
            -P -r '\@RG\\tID:$rn\\tSM:$sm\\tLB:$lb\\tPL:$pl\\tPU:lane' \\
            $d03/$rn.1.sai $d03/$rn.2.sai $f1a $f1b \\
            | samtools view -Sb - > $d03/$rn.bam", 1);
        runCmd("java -Xmx8g -jar $picard/FixMateInformation.jar \\
            TMP_DIR=$DIR_tmp VALIDATION_STRINGENCY=LENIENT \\
            INPUT=$d03/$rn.bam OUTPUT=$d04/$rn.bam SORT_ORDER=coordinate", 1);
        runCmd("samtools index $d04/$rn.bam", 1);
        runCmd("bamStat -i $d04/$rn.bam -o $d04/$rn", 1);
    }
}
sub pipe_lib {
    my ($dir, $f_lb, $lb) = @_;
    my $dirI = "$dir/04_fixmate_sort";

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

    my $d11 = "$dir/11_markdup";
    make_path($d11) unless -d $d11;
    my $input_str = join(" ", map {"INPUT=$_"} @fis);
    runCmd("java -Xmx10g -jar $picard/MarkDuplicates.jar \\
        VALIDATION_STRINGENCY=LENIENT TMP_DIR=$DIR_tmp \\
        $input_str OUTPUT=$d11/$lb.bam METRICS_FILE=$d11/$lb.dup.txt", 1);
    
    my $d12 = "$dir/12_dedup";
    make_path($d12) unless -d $d12;
    runCmd("bamtools filter -in $d11/$lb.bam -script \\
        $DIR_code/conf/dedup.json -out $d12/$lb.bam", 1);
    runCmd("bamStat -i $d12/$lb.bam -o $d12/$lb", 1);
    runCmd("samtools index $d12/$lb.bam", 1);
}
sub pipe_rerun {
    my ($dir, $f_rn, $f_bwa, $rni) = @_;

    my $t = readTable(-in=>$f_rn, -header=>1);
    my $h;
    for my $i (0..$t->nofRow-1) {
        my ($sm, $lb, $rns, $idxs) = $t->row($i);
        $h->{$lb} = $sm;
    }
    die "no library named '$lb'\n" unless exists $h->{$lb};
    my $sm = $h->{$lb};
 
    my $d12 = "$dir/12_dedup";
    my $d15 = "$dir/15_remap_reads";
    make_path($d15) unless -d $d15;
    runCmd("bamtools filter -in $d12/$lb.bam -script \\
        $DIR_code/conf/filter_lipe.json -out $d15/$lb.bam", 1);
    runCmd("samtools sort -n $d15/$lb.bam $d15/$lb.namesorted", 1);
    runCmd("$DIR_src/bamUtil/bin/bam bam2FastQ --in $d15/$lb.namesorted.bam \\
        --readname --firstOut $d15/$lb.1.fq --secondOut $d15/$lb.2.fq \\
        --unpairedOut $d15/$lb.fq", 1);
    
    my $d16 = "$dir/16_remap_aln";
    make_path($d16) unless -d $d16;
    runCmd("bwa aln -t 4 -n 0.01 $f_bwa $d15/$lb.1.fq > $d16/$lb.1.sai", 1);
    runCmd("bwa aln -t 4 -n 0.01 $f_bwa $d15/$lb.2.fq > $d16/$lb.2.sai", 1);
    runCmd("bwa sampe $f_bwa -r \\
        '\@RG\\tID:$lb\\tSM:$sm\\tLB:$lb\\tPL:ILLUMINA\\tPU:lane' \\
        $d16/$lb.1.sai $d16/$lb.2.sai $d15/$lb.1.fq $d15/$lb.2.fq \\
        | samtools view -Sb - > $d16/$lb.bam", 1);
    
    my $d17 = "$dir/17_remap_sorted";
    make_path($d17) unless -d $d17;
    runCmd("java -Xmx8g -jar $picard/FixMateInformation.jar \\
        TMP_DIR=$DIR_tmp VALIDATION_STRINGENCY=LENIENT \\
        INPUT=$d16/$lb.bam  OUTPUT=$d17/$lb.bam SORT_ORDER=coordinate", 1);
    
    my $d18 = "$dir/18_remap_markdup";
    make_path($d18) unless -d $d18;
    runCmd("java -Xmx10g -jar $picard/MarkDuplicates.jar \\
        TMP_DIR=$DIR_tmp VALIDATION_STRINGENCY=LENIENT \\
        REMOVE_DUPLICATES=true \\
        INPUT=$d17/$lb.bam OUTPUT=$d18/$lb.bam METRICS_FILE=$d18/$lb.dup.txt", 1);
    
    my $d19 = "$dir/19_remap_dedup";
    make_path($d19) unless -d $d19;
    runCmd("bamtools filter -in $d18/$lb.bam -script \\
        $DIR_code/conf/dedup.json -out $d19/$lb.bam", 1);
    runCmd("samtools index $d19/$lb.bam", 1);
    runCmd("bamStat -i $d19/$lb.bam -o $d19/$lb", 1);
}
sub pipe_sample {
    my ($dir, $f_sm, $f_ref, $sm) = @_;
    my $dirI1 = "$dir/12_dedup";
    my $dirI2 = "$dir/19_remap_dedup";

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
        my $fi = ($lb =~ /lipe/i) ? "$dirI2/$lb.bam" : "$dirI1/$lb.bam";
        die "$lb [$fi] is not there\n" unless -s $fi;
        push @fis, $fi;
    }

    my $d21 = "$dir/21_realigned";
    make_path($d21) unless -d $d21;
    my $input_str = join(" ", map {"-I $_"} @fis);
    runCmd("java -Xmx10g -Djava.io.tmpdir=$DIR_tmp -jar $gatk/GenomeAnalysisTK.jar \\
        -T RealignerTargetCreator -R $f_ref $input_str \\
        -o $d21/$sm.intervals", 1);
    runCmd("java -Xmx10g -Djava.io.tmpdir=$DIR_tmp -jar $gatk/GenomeAnalysisTK.jar \\
        -T IndelRealigner -R $f_ref $input_str \\
        -targetIntervals $d21/$sm.intervals -o $d21/$sm.bam \\
        -LOD 0.4 --maxReadsForRealignment 20000 --maxReadsInMemory 200000", 1);
    runCmd("bamStat -i $d21/$sm.bam -o $d21/$sm", 1);
}
sub pipe_vnt {
    my ($dir, $f_sm, $f_ref, $sm) = @_;
    my $d21 = "$dir/21_realigned";

    my $d31 = "$dir/31_bcf_raw";
    make_path($d31) unless -d $d31;
    runCmd("samtools mpileup -gD -f $f_ref $d21/$sm.bam -q 10 > $d31/$sm.bcf", 1);
}
 
sub stats {
    my ($dir, $f_sm, $fo1, $fo2) = @_;
    my $dirI = "$dir/21_realigned";
    my $t = readTable(-in=>$f_sm, -header=>1);
    my $to1 = Data::Table->new([], [qw/rg total unmapped unpaired unpaired_dedup unpaired_uniq  paired paired_dedup paired_uniq paired_proper/]);
    my $to2 = Data::Table->new([], [qw/rg is cnt/]);
    for my $i (0..$t->nofRow-1) {
        my ($sm, $lbs) = $t->row($i);
        my $fa = "$dirI/$sm.tbl";
        next unless -s $fa;
        my $ta = readTable(-in=>$fa, -header=>1);
        for my $j (0..$ta->nofRow-1) {
            $to1->addRow($ta->rowRef($j));
        }

        my $fb = "$dirI/$sm.isd.tbl";
        my $tb = readTable(-in=>$fb, -header=>1);
        for my $j (0..$tb->nofRow-1) {
            $to2->addRow($tb->rowRef($j));
        }
    }
    
    open(FH1, ">$fo1") or die "cannot open $fo1 for writing\n";
    open(FH2, ">$fo2") or die "cannot open $fo2 for writing\n";
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
    
    my $d03 = "$dir/03_aln";
    my $d06 = "$dir/06_fixmate_sort";
    open(FH1, ">$fo1") or die "Cannot open $fo1 for writing\n";
    print FH1 join("\t", qw/idx sm lb rn aln fixmate_sort/)."\n";
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
   
    my $d11 = "$dir/11_dedup";
    my $d18 = "$dir/18_remap_dedup";
    open(FH2, ">$fo2") or die "Cannot open $fo2 for writing\n";
    print FH2 join("\t", qw/sm lb dedup rns/)."\n";
    my $tl = readTable(-in=>$f_lib, -header=>1);
    for my $i (0..$tl->nofRow-1) {
        my ($sm, $lb, $rns, $idxs) = $tl->row($i);
        my $fi = ($lb =~ /lipe/i) ? "$d18/$lb.bam" : "$d11/$lb.bam";
        my $tag = check_bam($fi);
        print FH2 join("\t", $sm, $lb, $tag, $rns, $idxs)."\n";
        print "  checking lib $lb\r";
    }
    print "\n";
    close FH2;

    my $d21 = "$dir/21_realigned";
    open(FH3, ">$fo3") or die "Cannot open $fo3 for writing\n";
    print FH3 join("\t", qw/sm realigned lbs/)."\n";
    my $ts = readTable(-in=>$f_sam, -header=>1);
    for my $i (0..$ts->nofRow-1) {
        my ($sm, $lbs) = $ts->row($i);
        my $tag = check_bam("$d21/$sm.bam");
        print FH3 join("\t", $sm, $tag, $lbs)."\n";
        print "  checking sample $sm\r";
    }
    print "\n";
    close FH3;
}


#my $cmd1 = "novoalign -d $DIR_db/novoalign/mt_35 -f $f1a -f2 $f1b -r E 10 -t 99 -o SAM | samtools view -Sb - > $f1";

