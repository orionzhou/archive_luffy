#!/usr/bin/perl
use strict;
use lib ($ENV{"SCRIPT_HOME_PERL"});
use Data::Dumper;
use File::Basename;
use File::Path qw/make_path remove_tree/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;

use InitPath;
use Common;

my $dir = "$DIR_misc3/ncgr_fastq";
my $f01 = "$dir/01_fastq.tbl";
#check_fastq($f01);
my $d03 = "$dir/03_fastq_stats";
#run_fastqc($f01, $d03);
my $f04 = "$dir/04_fastq_stats.tbl";
get_fastq_stats($f01, $d03, $f04);
my $f11 = "$dir/11_library.tbl";
my $f21 = "$dir/21_sample.tbl";
#cat_lib_info($f01, $f11, $f21);

sub check_fastq {
    my ($fi, $fo) = @_;
    my $dir = dirname($fi);
    my $ti = readTable(-in=>$fi, -header=>1);
    for my $i (0..$ti->nofRow-1) {
        my ($idx, $sm, $lb, $run, $pl, $pi, $dir_rel) = $ti->row($i);
        my $id = sprintf "$sm\_$lb\_%02d", $run;
        my $fq1 = "$dir/$dir_rel/$id.1.fq.gz";
        my $fq2 = "$dir/$dir_rel/$id.2.fq.gz";
        die "$fq1 is not there\n" unless -s $fq1;
        die "$fq2 is not there\n" unless -s $fq2;
    }
}

sub run_fastqc {
    my ($fi, $do) = @_;
    my $dir = dirname($fi);
    my $ti = readTable(-in=>$fi, -header=>1);
    my $f_bin = "/home/youngn/zhoup/Source/FastQC/fastqc";
    for my $i (0..$ti->nofRow-1) {
        my ($idx, $sm, $lb, $run, $pl, $pi, $dir_rel) = $ti->row($i);
        my $id = sprintf "$sm\_$lb\_%02d", $run;
        my $fq1 = "$dir/$dir_rel/$id.1.fq.gz";
        my $fq2 = "$dir/$dir_rel/$id.2.fq.gz";
        die "$fq1 is not there\n" unless -s $fq1;
        die "$fq2 is not there\n" unless -s $fq2;
        runCmd("$f_bin -o $do -t 2 $fq1 $fq2", 1);
    }
}
sub parse_fastqc {
    my ($fi) = @_;
    open(FHI, "<$fi") or die "cannot open $fi for reading\n";
    my $h;
    while(<FHI>) {
        chomp;
        my @ps = split "\t";
        if(/^Total Sequences/) {
            $h->{'n_seq'} = $ps[1];
        } elsif(/^Encoding/) {
            if($ps[1] =~ /illumina ([\d\.]+)/i) {
                $h->{'encoding'} = $1;
            } else {
                die "unknown illumina encoding: $ps[1]\n";
            }
        } elsif(/^Sequence length/) {
            $h->{'seq_len'} = $ps[1];
        } elsif(/^\%GC/) {
            $h->{'pct_GC'} = $ps[1];
        } elsif(/Total Duplicate Percentage/) {
            $h->{'pct_dup'} = sprintf "%.02f", $ps[1];
        }
    }
    close FHI;
    return $h;
}
sub get_fastq_stats {
    my ($fi, $dir, $fo) = @_;
    my $ti = readTable(-in=>$fi, -header=>1);
    $ti->addCol([('') x $ti->nofRow], "rl", 5);
    $ti->addCol([('') x $ti->nofRow], "n_seq");
    $ti->addCol([('') x $ti->nofRow], "encoding");
    $ti->addCol([('') x $ti->nofRow], "pct_GC");
    $ti->addCol([('') x $ti->nofRow], "pct_dup1");
    $ti->addCol([('') x $ti->nofRow], "pct_dup2");
    for my $i (0..$ti->nofRow-1) {
        my ($idx, $sm, $lb, $run, $pl, $pi, $dir_rel) = $ti->row($i);
        my $id = sprintf "$sm\_$lb\_%02d", $run;
        my $fs1 = "$dir/$id.1.fq_fastqc/fastqc_data.txt";
        my $fs2 = "$dir/$id.2.fq_fastqc/fastqc_data.txt";
        next unless -s $fs1;
        my $s1 = parse_fastqc($fs1);
        my $s2 = parse_fastqc($fs2);
        die "n_seq not equal for $idx:$id\n" unless $s1->{'n_seq'} == $s2->{'n_seq'}; 
        die "encoding not same for $idx:$id\n" unless $s1->{'encoding'} == $s2->{'encoding'}; 
        die "seq_len not equal for $idx:$id\n" unless $s1->{'seq_len'} == $s2->{'seq_len'}; 
        $ti->setElm($i, "rl", $s1->{'seq_len'});
        $ti->setElm($i, "n_seq", $s1->{'n_seq'});
        $ti->setElm($i, "encoding", $s1->{'encoding'});
        $ti->setElm($i, "pct_GC", $s1->{'pct_GC'});
        $ti->setElm($i, "pct_dup1", $s1->{"pct_dup"});
        $ti->setElm($i, "pct_dup2", $s2->{"pct_dup"});
    }
    open(FHO, ">$fo") or die "cannot write to $fo\n";
    print FHO $ti->tsv(1);
    close FHO;
}
 
sub cat_lib_info {
    my ($fi, $fo1, $fo2) = @_;
    my $ti = readTable(-in=>$fi, -header=>1);
    my $h;
    for my $i (0..$ti->nofRow-1) {
        my ($idx, $sm, $lb, $run) = $ti->row($i);
        $h->{$sm}->{$lb} ||= [];
        push @{$h->{$sm}->{$lb}}, [$run, $idx];
    }

    open(FHO1, ">$fo1") or die "cannot write to $fo1\n";
    print FHO1 join("\t", qw/sm lb runs idxs/)."\n";
    for my $sm (sort(keys(%$h))) {
        for my $lb (sort(keys(%{$h->{$sm}}))) {
            my $runs = join(",", map {$_->[0]} @{$h->{$sm}->{$lb}});
            my $idxs = join(",", map {$_->[1]} @{$h->{$sm}->{$lb}});
            print FHO1 join("\t", $sm, $lb, $runs, $idxs)."\n";
        }
    }
    close FHO1;
    
    open(FHO2, ">$fo2") or die "cannot write to $fo2\n";
    print FHO2 join("\t", qw/sm lbs/)."\n";
    for my $sm (sort(keys(%$h))) {
        my $lbs = join(",", keys(%{$h->{$sm}}));
        print FHO2 join("\t", $sm, $lbs)."\n";
    }
    close FHO2;
}


