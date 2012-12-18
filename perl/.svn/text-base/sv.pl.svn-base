#!/usr/bin/perl -w
use strict; use Init; use Common; use Getopt::Long; use Run; 
use Bio::Seq; use Bio::SeqIO; use Graph; use Bam;
use Readfile; use Writefile; use Annotate; use Align; use Parser; use Medicago;
use Gff; use CircosConf; use Draw; use Parser; use Seq; use Convert; use GeneModel;
use Time::HiRes qw/gettimeofday tv_interval/; use Data::Dumper;use Path::Class; 
use List::Util qw/min max sum/; use POSIX qw/ceil floor/; use Math::Round qw/round/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
my $refDb = "mt_35";
my $f_genome = file($DIR_Genome, $refDb, "41_genome.fa");
my ($program, $beg, $end) = ('', 0, -1);
GetOptions('b=i'=>\$beg, 'e=i'=>\$end, 'p=s'=>\$program);
my $ids = get_acc_ids("acc26");

my $dir = dir($DIR_Repo, $refDb, "40_sv");
my $d01 = file($dir, "01_pos_sorted"); # ln -sf ../11_pipe_bwa/13_dup_removed 01_pos_sorted
my $d02 = file($dir, "02_rn_sorted"); # ln -sf ../11_pipe_bwa/21_rn_sorted 02_rn_sorted
my $f03 = file($dir, "03_stat.tbl"); # ln -sf ../15_pipe_bam/13_stat.tbl 03_stat.tbl
run_hydra($ids, $dir) if $program eq "hydra";
run_pindel($ids, $dir, $beg, $end) if $program eq "pindel";
run_crest($ids, $dir, $beg, $end) if $program eq "crest";
sum_crest($ids, $dir) if $program eq "sum_crest";
run_genome_strip($ids, $dir, $beg, $end) if $program eq "genome_strip";
run_tigrasv($ids, $dir) if $program eq "tigrasv";

sub run_hydra {
    my ($dir) = @_;
    my $d01 = file($dir, "01_pos_sorted"); 
    my $d02 = file($dir, "02_rn_sorted"); 
    my $f03 = file($dir, "03_stat.tbl"); 
    my $t = readTable(-in=>$f03, -header=>1);
    
    my $d11 = file($dir, "11_stretched");
    my $d13 = file($dir, "13_hydra");
    my $d14 = file($dir, "14_hydra_filtered");
    
    my $hydra = "/project/youngn/zhoup/Source/Hydra-Version-0.5.3/bin/hydra";
    for my $id (@$ids) {
        my $t2 = $t->match_pattern("\$_->[0] eq '$id'");
        die "cannot find $id in table\n" unless $t2->nofRow == 1;
        my ($id, $runs, $is_mld, $is_mno) = map {$t2->elm(0, $_)} qw/id runs is_mld is_mno/;
        my $rgs_str = join(" ", map {sprintf "-g %s_%02d", $id, $_} split(" ", $runs));
        runCmd("bamPickStretched -i $d01/$id.bam -o $d11/$id.bed -m $is_mno $rgs_str -r chr5");
        runCmd("$hydra -in $d11/$id.bed -out $d13/$id -mld $is_mld -mno $is_mno -ms 5");
    }
}
sub run_pindel {
    my ($ids, $dir, $beg, $end) = @_;
    my $d01 = file($dir, "01_pos_sorted"); 
    my $d02 = file($dir, "02_rn_sorted"); 
    my $f03 = file($dir, "03_stat.tbl"); 
    my $t = readTable(-in=>$f03, -header=>1);
  
    my $d21 = file($dir, "21_orphan");
    my $d31 = file($dir, "31_pindel");
    
    my $chr = "chr5";
    for my $i ($beg..$end) {
        my $id = sprintf "HM%03d", $i;
        my $t2 = $t->match_pattern("\$_->[0] eq '$id'");
        die "cannot find $id in table\n" unless $t2->nofRow == 1;
        my ($runs, $is_mld, $is_mno) = map {$t2->elm(0, $_)} qw/runs is_mld is_mno/;
#    runCmd("bamPickOrphan -i $d01/$id.bam -o $d21/01_read_id/$id.txt -r $chr", 1);
#    runCmd("bamPickReads -i $d21/01_read_id/$id.txt -b $d02/$id.bam -o $d21/02_info/$id.txt -t $id", 1);
    }
    
    my $f_str = join(" ", map {sprintf "21_orphan/11_pindel/%s.txt", $_} @$ids);
#  print "cat $f_str > 21_orphan/11_pindel.txt\n";
#  print "pindel -f $f_genome -p 21_orphan/11_pindel.txt -o $d24/01 -c $chr -T 4\n";
#  print "pindel2vcf -p 31_pindel/01_D -r ../01_reference/41_genome.fa -R Mt3.5 -d 20110501\n";
#  parse_pindel("$d31/01_D", "$d31/11.tbl");

#  print "cov_window -i $d11/04_bp.tbl -o $d11/11_cov.tbl -t acc26 -c 1\n";
}
sub run_crest {
    my ($ids, $dir, $beg, $end) = @_;
    my $d01 = file($dir, "01_pos_sorted"); 
    my $d02 = file($dir, "02_rn_sorted"); 
    my $d33 = file($dir, "33_crest");
    my $f_genome_2bit = file($DIR_Db, "blat/mt_35.2bit");

#gfClient `cat $m/pbs/host` 1986 $DIR_Db/blat test.fa test.psl
    my $chr = "chr5";
    for my $i ($beg..$end) {
        my $id = sprintf "HM%03d", $i;
        my $cmd = "extractSClip.pl -i $d01/$id.bam --ref_genome $f_genome -o $d33/01_in -r $chr";
#    runCmd($cmd);
        $cmd = "CREST.pl -f $d33/01_in/$id.bam.$chr.cover -d $d01/$id.bam \\
            --ref_genome $f_genome --2bitdir $DIR_Db/blat -t $f_genome_2bit -o $d33/03_predSV \\
            --blatserver elmob111 --blatport 1986 -r $chr";
#    runCmd($cmd, 0);
    }
} 
sub sum_crest {
    my ($ids, $dirP) = @_;
    my $dir = dir($dirP, "33_crest");
    my $d03 = dir($dir, "03_predSV");
    my $f11 = file($dir, "11_sum.tbl");
    my $fh = new IO::File $f11, "w";
    print $fh join("\t", qw/acc chr_l pos_l strand_l n_sc_l chr_r pos_r strand_r n_sc_r type cov_l cov_r len_l len_r pct_idty_l pct_reads_nu_l pct_idty_r pct_reads_nu_r beg_con chr_b beg_con_map end_con chr_e end_con_map seq beg_con2 chr_b2 beg_con_map2 end_con2 chr_e2 end_con_map2 seq2/)."\n";
    opendir(my $dh, $d03) || die "cannot open $d03\n";
    for my $fn (sort readdir $dh) {
        if($fn =~ /^(HM\d+)\.bam\.predSV/) {
            my $acc = $1;
            my $fhi = new IO::File "$d03/$fn", "r";
            while(<$fhi>) {
                chomp;
                my @ps = split("\t", $_);
                if(@ps == 24) {
                    push @ps, ("") x 7;
                } elsif(@ps == 31) {
                } else {
                    die "cannot handle ".@ps." columns\n";
                }
                die join("\t", @ps) if @ps ne 31;
                print $fh join("\t", $acc, @ps)."\n";
            }
        }
    }
    closedir $dh;
}
sub run_tigrasv {
    my ($ids, $dir, $beg, $end) = @_;
    my $d01 = file($dir, "01_pos_sorted"); 
    my $d02 = file($dir, "02_rn_sorted"); 
    my $d37 = file($dir, "37_tigra");
    print "./tigra-sv -R $f_genome $d37/02_sv.tbl -f $d37/01_bamlist.tbl -o $d37/03.fa -r $d37/04_ref.fa\n";
}
sub run_genome_strip {
    my ($ids, $dir, $beg, $end) = @_;
    my $d01 = file($dir, "01_pos_sorted"); 
    my $d02 = file($dir, "02_rn_sorted"); 
    my $d35 = file($dir, "35_genome_strip");
    my $f_genome = file($DIR_Db, "bwa/mt_35.fasta");
    my $f_genome_mask = file($DIR_Db, "bwa/mt_35.mask.fasta");
    
    my $chr = "chr5";
    my $lib = "$svtoolkit/lib";
    my $qscript = "$svtoolkit/qscript";
    for my $i ($beg..$end) {
        my $id = sprintf "HM%03d", $i;
        my $cmd = "java -Xmx3g -cp $lib/gatk/Queue.jar:$lib/SVToolkit.jar:$lib/gatk/GenomeAnalysisTK.jar \\
            org.broadinstitute.sting.queue.QCommandLine \\
            -S $qscript/SVPreprocess.q -S $qscript/SVQScript.q \\
            -gatk $lib/gatk/GenomeAnalysisTK.jar \\
            -cp $lib/SVToolkit.jar:$lib/gatk/GenomeAnalysisTK.jar \\
            -configFile \$conf/svtoolkit/mt.txt -tempDir $DIR_Tmp \\
            -md $d35/01_preprocess -R $f_genome -genomeMaskFile $f_genome_mask \\
            -I $d01/$id.bam -run";
#    runCmd($cmd); 
        $cmd = "java -Xmx3g -cp $lib/gatk/Queue.jar:$lib/SVToolkit.jar:$lib/gatk/GenomeAnalysisTK.jar \\
            org.broadinstitute.sting.queue.QCommandLine \\
            -S $qscript/SVDiscovery.q -S $qscript/SVQScript.q \\
            -gatk $lib/gatk/GenomeAnalysisTK.jar \\
            -cp $lib/SVToolkit.jar:$lib/gatk/GenomeAnalysisTK.jar \\
            -configFile \$conf/svtoolkit/mt.txt -tempDir $DIR_Tmp \\
            -runDirectory $d35 \\
            -md $d35/01_preprocess -R $f_genome -genomeMaskFile $f_genome_mask \\
            -I $d01/$id.bam -O $d35/02_alt/$id.vcf \\
            -minimumSize 100 -maximumSize 1000000 -windowSize 10000000 -windowPadding 10000 -run";
#    runCmd($cmd); 
        print $cmd."\n";
        $cmd = "java -Xmx3g -Djava.library.path=$svtoolkit/bwa \\
            -cp $lib/gatk/Queue.jar:$lib/SVToolkit.jar:$lib/gatk/GenomeAnalysisTK.jar \\
            org.broadinstitute.sting.queue.QCommandLine \\
            -S $qscript/SVAltAlign.q -S $qscript/SVQScript.q \\
            -gatk $lib/gatk/GenomeAnalysisTK.jar \\
            -cp $lib/SVToolkit.jar:$lib/gatk/GenomeAnalysisTK.jar \\
            -configFile \$conf/svtoolkit/mt.txt -tempDir $DIR_Tmp \\
            -md $d35/01_preprocess -runDirectory $d35 -R $f_genome \\
            -vcf $d35/02_alt/$id.vcf \\
            -I $d01/$id.bam -O $d35/03_alt_align/$id.bam -run";
    }
}




