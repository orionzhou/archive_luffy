#!/usr/bin/perl
use strict; use Init; use Common; use Readfile; use Localdb;
use Bio::Seq; use Path::Class; use Data::Dumper; use Medicago;
use List::Util qw/min max sum/; use POSIX qw/ceil floor/;
use Time::HiRes qw/gettimeofday tv_interval/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
my $refDb = "mt_35";
my $opt = "acc84";
my $ids = get_acc_ids($opt);

my $dirW = dir($DIR_Repo, $refDb, "30_vnt_$opt", "11_snps");
my $chrs = [map { "chr".$_ } (1..8)];
my $d01 = dir($dirW, "01_ssp");
my $d02 = dir($dirW, "02_plain");
#sspExpand($d01, $d02, $ids, $chrs);
my $d05 = dir($dirW, "05_ped");
#ssp2Ped($d01, $d05, $chrs);
my $d11 = dir($dirW, "11_merged");
#ssp2Tbl($d01, $d11, $ids, $chrs, $refDb);
my $d31 = dir($dirW, "31_vcf");
my $d32 = dir($dirW, "32_snpEff");
#snpEff($d01, $d31, $d32, $chrs);

sub sspExpand {
    my ($dirI, $dirO, $ids, $chrs) = @_;
    for my $chr (@$chrs) {
#    next unless $chr eq "chr1";
        my $fi = file($dirI, "$chr.txt");
        for my $id (@$ids) {
            system("mkdir -p $dirO/$id") unless -d "$dirO/$id";
        }
        runCmd("sspConvert -i $fi -o $dirO -f plain");
    }
}
sub ssp2Tbl {
    my ($dirI, $dirO, $ids, $chrs, $refDb) = @_;
    system("mkdir -p $dirO") unless -d $dirO;
    my @fns;
    for my $chr (@$chrs) {
        my $fi = file($dirI, "$chr.txt");
        my $fo = file($dirO, "$chr.tbl");
        runCmd("sspConvert -i $fi -o $fo -f table");
        runCmd('sed -i \'s/^\([0-9\.][0-9\.]*\)/'.$chr.'\t\1\t\1/\' '.$fo);
        push @fns, $fo;
    }
    my $fo = file($dirO, "01.tbl.gz");
    runCmd("cat ".join(" ", @fns)." | bgzip > $fo");
    runCmd("tabix -f -s1 -b2 -e3 $fo");
    runCmd("rm ".join(" ", @fns));
}
sub ssp2Ped {
    my ($dirI, $dirO, $chrs, $refDb) = @_;
    system("mkdir -p $dirO") unless -d $dirO;
    for my $chr (@$chrs) {
        my $fi = file($dirI, "$chr.txt");
        my $fo = file($dirO, "$chr.ped");
        runCmd("sspFilter -i $fi -o $fo -b 1 -f ped");
    }
}
sub snpEff {
    my ($d01, $d31, $d32, $chrs) = @_;
    system("mkdir -p $d31") unless -d $d31;
    system("mkdir -p $d32") unless -d $d32;
#  $chrs = [map { "chr".$_ } (4)];
    for my $chr (@$chrs) {
        my $f01 = "$d01/$chr.txt";
        my $f31 = "$d31/$chr.vcf";
        my $f32 = "$d32/$chr.tbl";
#    runCmd("ssp2Vcf -i $f01 -o $f31 -r HM000 -c $chr");
        my $cmd = "java -Xmx5g -jar \$src/snpEff_2_0_3/snpEff.jar eff mt_35 -c \$m/conf/snpEff.config -i vcf -ud 0 $f31 -s $d32/$chr.html > $f32";
        system($cmd);
    }
}


