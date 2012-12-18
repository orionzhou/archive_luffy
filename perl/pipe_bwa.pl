#!/usr/bin/perl
use strict;
use Init;
use Common;
use Medicago;
use Bam;
use Getopt::Long;
use Data::Dumper;
use List::Util qw/min max sum/; use POSIX qw/ceil floor/;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;

my ($program, $beg, $end) = ('', 0, -1);
GetOptions('b=i'=>\$beg, 'e=i'=>\$end, 'p=s'=>\$program);

my $f_genome = "$DIR_Genome/mt_35/41_genome.fa";
my ($ids) = get_mt_ids("acc84");

my $dir = "$DIR_Misc3/hapmap/11_pipe_bwa";
run_pipe_bwa($dir, $beg, $end) if $program eq "pipe_bwa";
update_status($dir) if $program eq "update";
sub run_pipe_bwa {
    my ($dir, $beg, $end) = @_;
    my $f01 = "$dir/01_fastq.tbl";
    my $d03 = "$dir/03_bwa";
    my $d06 = "$dir/06_pos_sorted";
    
    my $t = readTable(-in=>$f01, -header=>1);
    for my $i ($beg..$end) {
        my ($idx, $sm, $rg, $lb, $id_lb, $id_run, $pl, $rl, $pi, $f1a, $f1b) = $t->row($i);
        die "$f1a is not there\n" unless -s $f1a;
        die "$f1b is not there\n" unless -s $f1b;
        my ($f3a, $f3b) = map {"$d03/$rg.$_.sai"} (1,2);
        my $f3 = "$d03/$rg.bam";
        my $f6 = "$d06/$rg.bam";
        runCmd("bwa aln -t 4 -n 8 \$data/db/bwa/mt_35 $f1a > $f3a");
        runCmd("bwa aln -t 4 -n 8 \$data/db/bwa/mt_35 $f1b > $f3b");
=cut
=cut
        runCmd("bwa sampe \$data/db/bwa/mt_35 -r '\@RG\\tID:$rg\\tSM:$sm\\tLB:$lb\\tPL:$pl' $f3a $f3b $f1a $f1b | samtools view -Sb - > $f3"."\n");
        runCmd("java -Xmx7g -jar $picard/FixMateInformation.jar "
            . "TMP_DIR=$DIR_Tmp VALIDATION_STRINGENCY=LENIENT "
            . "INPUT=$f3 OUTPUT=$f6 SORT_ORDER=coordinate", 1);
        runCmd("samtools index $f6", 1);
    }
}
sub update_status {
    my ($dir) = @_;
    my $f01 = "$dir/01_fastq.tbl";
    my $d03 = "$dir/03_bwa";
    my $d06 = "$dir/06_pos_sorted";
    my $f07 = "$dir/07_status.tbl";
    
    open(FH, ">$f07");
    print FH join("\t", qw/idx sm rg lb id_run pi 03_bwa 06_pos_sorted/)."\n";
    my $t = readTable(-in=>$f01, -header=>1);
    for my $i (0..$t->nofRow-1) {
        my ($idx, $sm, $rg, $lb, $id_lb, $id_run, $pl, $pi, $f1a, $f1b) = $t->row($i);
        my $tag03 = check_bam("$d03/$rg.bam");
        my $tag06 = check_bam("$d06/$rg.bam");
        $tag06 = 0 if $tag03 == 0;
        print FH join("\t", $idx, $sm, $rg, $lb, $id_run, $pi, $tag03, $tag06)."\n";
        printf "  checking %s: %3d/%3d done...\r", $rg, $i+1, $t->nofRow;
    }
    print "\n";
    close FH;
}


$dir = "$DIR_Misc3/hapmap/13_pipe_novoalign";
run_pipe_novoalign($dir, $beg, $end) if $program eq "pipe_novoalign";
sub run_pipe_novoalign {
    my ($f_fq, $dir_fq, $dir, $beg, $end) = @_;
    my $t = readTable(-in=>$f_fq, -header=>1);
    my $d01 = "$dir/01";
    my $d02 = "$dir/02_pos_sorted";
    for my $i ($beg..$end) {
        my ($id, $cnt_run, $runs) = $t->row($i);
        for my $run (split(" ", $runs)) {
            my $pre = "$id\_$run";
            my ($f1a, $f1b) = map {"$dir_fq/$pre.$_.fq.gz"} (1,2);
            die "$f1a is not there\n" unless -s $f1a;
            die "$f1b is not there\n" unless -s $f1b;
            my $f1 = "$d01/$pre.bam";
            my $cmd1 = "novoalign -d $DIR_Db/novoalign/mt_35 -f $f1a -f2 $f1b -r E 10 -t 99 -o SAM | samtools view -Sb - > $f1";
            print "$cmd1\n";
#      runCmd($cmd1, 1);
            
            my $f2 = "$d02/$pre.bam";
            my $cmd2 = "java -Xmx4g -jar $picard/FixMateInformation.jar "
                . "TMP_DIR=$DIR_Tmp VALIDATION_STRINGENCY=LENIENT "
                . "INPUT=$f1 OUTPUT=$f2 SORT_ORDER=coordinate";
            print "$cmd2\n";
#      runCmd($cmd2, 1);
#      runCmd("samtools index $f2", 1);
        }
    }
}


