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

my ($program, $beg, $end) = ('', 0, -1);
GetOptions('beg|b=i'=>\$beg, 'end|e=i'=>\$end, 'program|p=s'=>\$program);

my $f_ref = "$DIR_genome/Mtruncatula_4.0/01_refseq.fa";
my $f_bwa = "$DIR_db/bwa/mt_40";
my $f_lst = "$DIR_misc3/ncgr_fastq/04_fastq_stats.tbl";
my $dir = "$DIR_misc3/hapmap_mt40";
$dir = "$dir/11_pipe_bwa";

pipe_bwa_run($dir, $beg, $end, $f_lst, $f_bwa) if $program eq "run";
pipe_bwa_update($dir, $f_lst) if $program eq "update";
sub pipe_bwa_run {
    my ($dir, $beg, $end, $f_lst, $f_bwa) = @_;
    my $d03 = "$dir/03_bwa";
    make_path($d03) unless -d $d03;
    my $d06 = "$dir/06_pos_sorted";
    make_path($d06) unless -d $d06;
    
    my $dir_abs = dirname($f_lst);
    my $t = readTable(-in=>$f_lst, -header=>1);
    for my $i (0..$t->nofRow-1) {
        my ($idx, $sm, $lb, $run, $pl, $rl, $pi, $dir_rel, $n_seq, $encoding) = $t->row($i);
        next if $idx < $beg || $idx > $end;
        my $id = sprintf "$sm\_$lb\_%02d", $run;
        my $f1a = "$dir_abs/$dir_rel/$id.1.fq.gz";
        my $f1b = "$dir_abs/$dir_rel/$id.2.fq.gz";
        die "$f1a is not there\n" unless -s $f1a;
        die "$f1b is not there\n" unless -s $f1b;

        my $f3a = "$d03/$id.1.sai";
        my $f3b = "$d03/$id.2.sai";
        my $f3  = "$d03/$id.bam";
        my $tag = ($encoding < 1.8 & $encoding >= 1.3) ? "-I" : "";
        runCmd("bwa aln -t 4 -n 0.06 $tag $f_bwa $f1a > $f3a", 1);
        runCmd("bwa aln -t 4 -n 0.06 $tag $f_bwa $f1b > $f3b", 1);

        my $f6 = "$d06/$id.bam";
        runCmd("bwa sampe $f_bwa -r \\
            '\@RG\\tID:$idx\\tSM:$sm\\tLB:$lb\\tPL:$pl' \\
            $f3a $f3b $f1a $f1b | samtools view -Sb - > $f3", 1);
        runCmd("java -Xmx8g -jar $picard/FixMateInformation.jar \\
            TMP_DIR=$DIR_tmp VALIDATION_STRINGENCY=LENIENT \\
            INPUT=$f3 OUTPUT=$f6 SORT_ORDER=coordinate", 1);
        runCmd("samtools index $f6", 1);
    }
}
sub pipe_bwa_update {
    my ($dir, $f_lst) = @_;
    my $d03 = "$dir/03_bwa";
    my $d06 = "$dir/06_pos_sorted";
    my $f07 = "$dir/07_status.tbl";
    
    open(FH, ">$f07") or die "Cannot open $f07 for writing\n";
    print FH join("\t", qw/idx sm lb run 03_bwa 06_pos_sorted/)."\n";
    my $t = readTable(-in=>$f_lst, -header=>1);
    for my $i (0..$t->nofRow-1) {
        my ($idx, $sm, $lb, $run, $pl, $rl, $pi, $dir_rel, $n_seq, $encoding) = $t->row($i);
        my $id = sprintf "$sm\_$lb\_%02d", $run;
        my $tag03 = check_bam("$d03/$id.bam");
        my $tag06 = check_bam("$d06/$id.bam");
        $tag06 = 0 if $tag03 == 0;
        print FH join("\t", $idx, $sm, $lb, $run, $tag03, $tag06)."\n";
        printf "  checking %4d\r", $idx;
    }
    print "\n";
    close FH;
}


#my $cmd1 = "novoalign -d $DIR_db/novoalign/mt_35 -f $f1a -f2 $f1b -r E 10 -t 99 -o SAM | samtools view -Sb - > $f1";

