#!/usr/bin/perl
use strict;
use Init;
use Common;
use Bam;
use Getopt::Long;
use Data::Dumper;
use List::Util qw/min max sum/; use POSIX qw/ceil floor/;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;

my ($program, $beg, $end) = ('', 0, -1);
GetOptions('beg|b=i'=>\$beg, 'end|e=i'=>\$end, 'program|p=s'=>\$program);

my $f_genome = "$DIR_Genome/mt_35/41_genome.fa";
my $dir = "$DIR_Misc3/hapmap/15_pipe_bam";

my $f01 = "$dir/01_sample.tbl";
#sample_info("$dir/../09_fastq.tbl", $f01);
sub sample_info {
    my ($fi, $fo) = @_;
    my $ti = readTable(-in=>$fi, -header=>1);
    open(FH, ">$fo");
    print FH join("\t", qw/idx sm lbs rgs pi pl/)."\n";
    my $ref = group($ti->colRef("sm"));
    my $i = 1;
    for my $sm (sort(keys(%$ref))) {
        my ($idx, $cnt) = @{$ref->{$sm}};
        my $ts = $ti->subTable([$idx..$idx+$cnt-1]);
        die ">1 sample for $sm\n" unless uniq($ts->col("sm")) == 1;
        my $sm = $ts->elm(0, "sm");
        my @lbs = $ts->col("lb");
        my $lb = join(",", uniq(@lbs));
        my @pis = $ts->col("pi");
        my $pi = join(",", uniq(@pis));
        my @rgs = $ts->col("rg");
        my $rg = join(",", @rgs);
        my $pl = $ts->elm(0, "pl");
        print FH join("\t", $i++, $sm, $lb, $rg, $pi, $pl)."\n";
    }
    close FH;
}

pipe_bam($dir, $beg, $end) if $program eq "pipe_bam";
sub pipe_bam {
    my ($dir, $beg, $end) = @_;
    my $f01 = "$dir/01_sample.tbl";
    my $d03 = "$dir/03_pos_sorted";
    my $d05 = "$dir/05_dup_marked";
    my $d07 = "$dir/07_realigned";
    my $d12 = "$dir/12_stat";
    system("mkdir -p $d05") unless -d $d05;
    system("mkdir -p $d07") unless -d $d07;
    system("mkdir -p $d12") unless -d $d12;
    my $t = readTable(-in=>$f01, -header=>1);
    for my $i ($beg..$end) {
        my ($idx, $sm, $lbs, $rgs, $pi, $pl) = $t->row($i-1);
        my $cmd;

        my @rgs = split ",", $rgs;
        my $input_str = join(" ", map {"INPUT=$d03/$_.bam"} @rgs);
        my ($f05, $f05b) = ("$d05/$sm.bam", "$d05/$sm.txt");
        $cmd = "java -Xmx12g -jar $picard/MarkDuplicates.jar \\
            VALIDATION_STRINGENCY=LENIENT TMP_DIR=$DIR_Tmp \\
            $input_str OUTPUT=$f05 METRICS_FILE=$f05b";
=cut
        runCmd($cmd);
        runCmd("samtools index $f05");
=cut
        
        my ($f07, $f07b) = ("$d07/$sm.bam", "$d07/$sm.intervals");
=cut
        $cmd = "java -Xmx12g -Djava.io.tmpdir=$DIR_Tmp -jar $gatk/GenomeAnalysisTK.jar \\
            -T RealignerTargetCreator -I $f05 -R $f_genome -o $f07b";
        runCmd($cmd);
        $cmd = "java -Xmx12g -Djava.io.tmpdir=$DIR_Tmp -jar $gatk/GenomeAnalysisTK.jar \\
            -T IndelRealigner -I $f05 -R $f_genome -targetIntervals $f07b -o $f07 \\
            -LOD 0.4 --maxReadsForRealignment 20000 --maxReadsInMemory 200000";
        runCmd($cmd);

=cut
        my $f12 = "$d12/$sm";
        runCmd("bamISD -i $f07 -o $f12");
    }
}

sort_readname($dir, $beg, $end) if $program eq "sort_readname";
sub sort_readname {
    my ($dir, $beg, $end) = @_;
    my $d13 = "$dir/13_dup_removed";
    my $d21 = "$dir/21_rn_sorted";
    for my $i ($beg..$end) {
        my $id = sprintf "HM%03d", $i;
        my $cmd = "java -Xmx4g -jar $picard/SortSam.jar \\
            VALIDATION_STRINGENCY=LENIENT TMP_DIR=$DIR_Tmp \\
            SORT_ORDER=queryname INPUT=$d13/$id.bam OUTPUT=$d21/$id.bam";
        runCmd($cmd);
    }
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




