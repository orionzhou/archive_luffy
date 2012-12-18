#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Init;
use Common; 
use Seq;
use Medicago;
use Bam;
use Data::Dumper;
use List::Util qw/min max sum/; use POSIX qw/ceil floor/; use Math::Round qw/round/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
my ($program, $beg, $end) = ('', 0, -1);
GetOptions('beg|b=i'=>\$beg, 'end|e=i'=>\$end, 'program|p=s'=>\$program);

my $f_genome = "$DIR_Genome/mt_35/41_genome.fa";
my ($ids) = get_mt_ids("acc84");

my $dir = "$DIR_Misc3/hapmap/03_pipe_ncgr";
my $d01 = "$dir/01_raw";
my $f02 = "$dir/02_readgroups.tbl";
my $d03 = "$dir/03_readname_sorted";
my $d04 = "$dir/04_processed";
my $d05 = "$dir/05_pos_sorted";
my $d09 = "$dir/09_fastq";

#collect_rg($ids, $d01, $f02, $d09);
pipe_ncgr($dir, $beg, $end) if $program eq "pipe_ncgr";

sub collect_rg {
    my ($ids, $d01, $f02, $d09) = @_;
    my $t = Data::Table->new([], [qw/sm rg lb id_lb id_run pl rl pi fq1 fq2 ds fn/]);
    for my $id (@$ids) {
        my $fi = "$d01/$id.bam";
        die "$fi is not there\n" unless -s $fi;
        open(PP, "samtools view -H $fi | grep '\@RG' |") or die "error piping\n";
        my $hl = {};
        while( <PP> ) {
            chomp;
            my @ps = split "\t";
            next unless $ps[0] eq '@RG';
            my ($fn, $pi, $ds) = ("") x 3;
            for my $ele (@ps[1..$#ps]) {
                my ($tag, $value) = split(":", $ele); 
                $fn = $value if $tag eq "ID";
                $pi = $value if $tag eq "PI";
                $ds = $value if $tag eq "DS";
            }
            my ($id_lb, $id_run) = (1, 1);
            if(exists($hl->{$pi})) {
                ($id_lb, $id_run) = @{$hl->{$pi}};
                $hl->{$pi}->[1] = ++$id_run;
            } else {
                $id_lb = scalar(keys(%$hl)) + 1;
                $hl->{$pi} = [$id_lb, $id_run];
            }
            my $lb = sprintf "%s_%02d", $id, $id_lb;
            my $rg = sprintf "%s_%02d_%02d", $id, $id_lb, $id_run;
            my ($fq1, $fq2) = map {"$d09/$rg.$_.fq.gz"} (1,2);
            $t->addRow([$id, $rg, $lb, $id_lb, $id_run, 'ILLUMINA', 90, $pi, $fq1, $fq2, $ds, $fn]);
        }
        printf "  collecting \@RG tags for %s...\r", $id;
    }
    print "\n";
    $t->sort('rg', 1, 0);
    open(FH, ">$f02");
    print FH $t->tsv(1);
    close FH;
}
sub pipe_ncgr {
    my ($dir, $beg, $end) = @_;
    my $d01 = "$dir/01_raw";
    my $f02 = "$dir/02_readgroups.tbl";
    my $d03 = "$dir/03_readname_sorted";
    my $d09 = "$dir/09_fastq";
    for my $i ($beg..$end) {
        my $id = sprintf "HM%03d", $i;
        my $f01 = "$d01/$id.bam";
#    runCmd("samtools sort -m 5000000000 -n $f01 $d03/$id");
        my $f03 = "$d03/$id.bam";
        runCmd("bam2Fastq -i $f03 -o $d09 -s $id -m $f02", 1);
=cut
        runCmd("samtools sort -m 5000000000 -n $d01/$id.bam $d03/$id", 1);
        runCmd("bamPreprocess -i $f03 -o $d04/$id", 1);
        runCmd("samtools sort -m 5000000000 $d04/$id $d05/$id", 1);
        runCmd("samtools index $d05/$id.bam", 1);
=cut
    }
}

my $f_lst_add = "/project/youngn/ncgr_fastq/file_list.tbl";
my $f_fq = "$DIR_Misc3/hapmap/09_fastq.tbl";
merge_fq_list($f02, $f_lst_add, $f_fq) if $program eq "merge_fq_list";
sub merge_fq_list {
    my ($fi, $f_add, $fo) = @_;
    my $ti = readTable(-in=>$fi, -header=>1);
    my $ta = readTable(-in=>$f_add, -header=>1);
    my $h;
    for my $i (0..$ti->nofRow-1) {
        my ($sm, $rg, $lb, $id_lb, $id_run, $pl, $rl, $pi, $fq1, $fq2, $ds, $fn) = $ti->row($i);
        die "[$rg/fastq1] $fq1 is not there\n" unless -s $fq1;
        die "[$rg/fastq2] $fq2 is not there\n" unless -s $fq2;
        $h->{$sm}->{$pi} ||= [$id_lb, []];
        push @{$h->{$sm}->{$pi}->[1]}, $id_run;
    }
    for my $i (0..$ta->nofRow-1) {
        my ($sm, $rg, $lb, $id_lb, $id_run, $pl, $rl, $pi, $fq1, $fq2, $ds, $fn) = $ta->row($i);
        die "[$rg/fastq1] $fq1 is not there\n" unless -s $fq1;
        die "[$rg/fastq2] $fq2 is not there\n" unless -s $fq2;
        if(exists($h->{$sm}->{$pi})) {
            $id_lb = $h->{$sm}->{$pi}->[0];
            my $ids_run = $h->{$sm}->{$pi}->[1];
            $id_run = max(@$ids_run) + 1;
            push @$ids_run, $id_run;
        } else {
            my @ids_lb = map {$_->[0]} values(%{$h->{$sm}});
            $id_lb = max(@ids_lb) + 1;
            $id_run = 1;
            $h->{$sm}->{$pi} = [$id_lb, [$id_run]];
        }
        $pl ||= "ILLUMINA";
        $lb = sprintf "%s_%02d", $sm, $id_lb;
        $rg = sprintf "%s_%02d_%02d", $sm, $id_lb, $id_run;
        $ti->addRow([$sm, $rg, $lb, $id_lb, $id_run, $pl, $rl, $pi, $fq1, $fq2, $ds, $fn]);
    }
    $ti->sort("rg", 1, 0);
    $ti->addCol([1..$ti->nofRow], "idx", 0);
    open(FH, ">$fo");
    print FH $ti->tsv(1);
    close FH;
}
sub process_fastq {
    my ($d10, $beg, $end) = @_;
    for ($beg..$end) {
        my $acc = sprintf "HM%03d", $_;
        my $fi1 = "$d10/$acc/s_3_1_sequence.txt.gz";
        my $fi2 = "$d10/$acc/s_3_2_sequence.txt.gz";
        open(my $fhi1, "zcat $fi1 |") or die "cannot open $fi1\n";
        open(my $fhi2, "zcat $fi2 |") or die "cannot open $fi2\n";
        my $i = 0;
        my $h = {};
        while( !eof($fhi1) && !eof($fhi2) ) {
            my (@ls1, @ls2);
            for my $i (0..3) {
                push @ls1, substr(readline($fhi1), 0, -1);
                push @ls2, substr(readline($fhi2), 0, -1);
            }
            $ls1[0] =~ /^\@([\w\:\#]+)\/([12]) run=([\w\_]+)/;
            my ($id, $pair, $run) = ($1, $2, $3);
            $ls2[0] =~ /^\@([\w\:\#]+)\/([12]) run=([\w\_]+)/;
            die "error matching reads \n $id $1 $pair $2 $run $3\n" unless $id eq $1 && $run eq $3 && $pair == 1 && $2 == 2;
            $ls1[2] =~ /^\+([\w\:\#]+)\/([12])/;
            die "error matching reads\n" unless $id eq $1 && $2 == 1; 
            $ls2[2] =~ /^\+([\w\:\#]+)\/([12])/;
            die "error matching reads\n" unless $id eq $1 && $2 == 2; 

            if(exists($h->{$run})) {
                my ($fho1, $fho2) = @{$h->{$run}};
                print $fho1 join("\n", "\@$id/1", $ls1[1], "+", $ls1[3])."\n";
                print $fho2 join("\n", "\@$id/1", $ls2[1], "+", $ls2[3])."\n";
            } else {
                $i ++;
                my $fo1 = sprintf "$d10/$acc\_%02d.1.fq", 20+$i;
                my $fo2 = sprintf "$d10/$acc\_%02d.2.fq", 20+$i;
                my $fho1 = new IO::File $fo1, "w";
                my $fho2 = new IO::File $fo2, "w";
                $h->{$run} = [$fho1, $fho2];
            }
        }
        for my $j (1..$i) {
            my $fo1 = sprintf "$d10/$acc\_%02d.1.fq", 20+$j;
            my $fo2 = sprintf "$d10/$acc\_%02d.2.fq", 20+$j;
            runCmd("gzip $fo1");
            runCmd("gzip $fo2");
        }
    }
}

sub run_mpileup {
#run_mpileup($d05, "$dir/21_bcf", "$dir/22_bcf_uniq", $ids);
    my ($dir_bam, $d01, $d02, $ids) = @_;
    my @chrs = map {"chr$_"} (1..8);
#  @chrs = map {"chr$_"} (5);
    for my $id (@$ids) {
        print "working on $id\n";
        my $f00 = file($dir_bam, "$id.bam");
        my $f01 = file($d01, "$id.bcf");
        my $f02 = file($d02, "$id.bcf");

        my @f1_chrs;
        my @f2_chrs;
        for my $chr (@chrs) {
            my $f1_chr = file($d01, "$id\_$chr.bcf");
            my $f2_chr = file($d02, "$id\_$chr.bcf");
            runCmd2("samtools mpileup -ugf $f_genome -r $chr $f00 | bcftools view -bcg - > $f1_chr");
            runCmd2("samtools mpileup -q1 -ugf $f_genome -r $chr $f00 | bcftools view -bcg - > $f2_chr");
            push @f1_chrs, $f1_chr;
            push @f2_chrs, $f2_chr;
        }
        runCmd2("bcftools cat ".join(" ", @f1_chrs)." > $f01");
        runCmd2("bcftools index $f01");
        runCmd2("bcftools cat ".join(" ", @f2_chrs)." > $f02");
        runCmd2("bcftools index $f02");
        runCmd2("rm ".join(" ", @f1_chrs));
        runCmd2("rm ".join(" ", @f2_chrs));
    }
}



