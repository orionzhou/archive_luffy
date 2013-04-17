#!/usr/bin/perl -w
use strict;
use lib ($ENV{"SCRIPT_HOME_PERL"});
use Path::Class;
use Getopt::Long;
use File::Basename;
use File::Path qw/make_path remove_tree/;
use Data::Dumper;
use InitPath;
use Common;
use Bam;
use Medicago;

my ($opt, $sm) = (('') x 2);
GetOptions('opt|p=s'=>\$opt, 'sm|s=s'=>\$sm);

my $sms = get_mt_ids('acc26');
my $f_ref = "$DIR_genome/Mtruncatula_4.0/01_refseq.fa";
my $f_bwa = "$DIR_db/bwa/mt_40";
my $f_sm = "$DIR_misc3/ncgr_fastq/21_sample.tbl";

my $dir = "$DIR_misc3/hapmap_mt40/40_sv";
my $d01 = "$dir/01_pos_sorted"; # ln -sf ../11_pipe_bwa/13_dup_removed 01_pos_sorted
my $d02 = "$dir/02_rn_sorted"; # ln -sf ../11_pipe_bwa/21_rn_sorted 02_rn_sorted
my $f03 = "$dir/03_stat.tbl"; # ln -sf ../15_pipe_bam/13_stat.tbl 03_stat.tbl
if($opt eq "hydra") {
    run_hydra($dir, $sms);
} elsif($opt eq "pindel") {
    run_pindel($dir, $sms);
} elsif($opt eq "crest") {
    run_crest($dir, $sms);
} elsif($opt eq "genome_strip") {
    run_genome_strip($dir, $sms);
} elsif($opt eq "tigrasv") {
    run_tigrasv($dir, $sms);
} else {
    die "unknown option: $opt\n";
}

sub run_hydra {
    my ($dir, $sms) = @_;
    my $d01 = "$dir/01_pos_sorted"; 
    my $d02 = "$dir/02_rn_sorted"; 
    my $f03 = "$dir/03_stat.tbl"; 
    my $t = readTable(-in=>$f03, -header=>1);
    
    my $d11 = "$dir/11_stretched";
    my $d13 = "$dir/13_hydra";
    my $d14 = "$dir/14_hydra_filtered";
    
    my $hydra = "\$src/Hydra-Version-0.5.3/bin/hydra";
    for my $sm (@$sms) {
        my $t2 = $t->match_pattern("\$_->[0] eq '$sm'");
        die "cannot find $sm in table\n" unless $t2->nofRow == 1;
        my ($sm, $rns, $is_mld, $is_mno) = map {$t2->elm(0, $_)} qw/sm rns is_mld is_mno/;
        my $rgs_str = join(" ", map {"-g $_"} split(" ", $rns));
        runCmd("bamPickStretched -i $d01/$sm.bam -o $d11/$sm.bed -m $is_mno $rgs_str -r chr5");
        runCmd("$hydra -in $d11/$sm.bed -out $d13/$sm -mld $is_mld -mno $is_mno -ms 5");
    }
}
sub run_pindel {
    my ($dir, $sms) = @_;
    my $d01 = "$dir/01_pos_sorted"; 
    my $d02 = "$dir/02_rn_sorted"; 
    my $f03 = "$dir/03_stat.tbl"; 
    my $t = readTable(-in=>$f03, -header=>1);
  
    $ENV{PATH} = "$ENV{PATH}:$pindel";
    my $d21 = "$dir/21_orphan";
    my $d31 = "$dir/31_pindel";
    
    my $chr = "chr5";
    for my $sm (@$sms) {
        my $t2 = $t->match_pattern("\$_->[0] eq '$sm'");
        die "cannot find $sm in table\n" unless $t2->nofRow == 1;
        my ($rns, $is_mld, $is_mno) = map {$t2->elm(0, $_)} qw/rns is_mld is_mno/;
        runCmd("bamPickOrphan -i $d01/$sm.bam -o $d21/01_read_id/$sm.txt -r $chr", 1);
        runCmd("bamPickReads -i $d21/01_read_id/$sm.txt -b $d02/$sm.bam -o $d21/02_info/$id.txt -t $id", 1);
    }
    
    my $f_str = join(" ", map {sprintf "21_orphan/11_pindel/%s.txt", $_} @$sms);
    print "cat $f_str > 21_orphan/11_pindel.txt\n";
    print "pindel -f $f_ref -p 21_orphan/11_pindel.txt -o $d24/01 -c $chr -T 4\n";
    print "pindel2vcf -p 31_pindel/01_D -r ../01_reference/41_genome.fa -R Mt3.5 -d 20110501\n";
#    parse_pindel("$d31/01_D", "$d31/11.tbl");

#    print "cov_window -i $d11/04_bp.tbl -o $d11/11_cov.tbl -t acc26 -c 1\n";
}
sub run_crest {
    my ($dir, $sms) = @_;
    my $d01 = "$dir/01_pos_sorted"; 
    my $d02 = "$dir/02_rn_sorted"; 
    my $f03 = "$dir/03_stat.tbl"; 
    my $t = readTable(-in=>$f03, -header=>1);
    
    $ENV{PATH} = "$ENV{PATH}:$crest";
    my $f_ref_2bit = "$DIR_db/blat/Mtruncatula_4.0.2bit";
    my $d33 = "$dir/33_crest";

#gfClient `cat $m/pbs/host` 1986 $DIR_Db/blat test.fa test.psl
    my $chr = "chr5";
    for my $sm (@$sms) {
        runCmd("extractSClip.pl -i $d01/$sm.bam --ref_genome $f_ref \\
            -o $d33/01_in -r $chr", 1);
        runCmd("CREST.pl -f $d33/01_in/$sm.bam.$chr.cover -d $d01/$sm.bam \\
            --ref_genome $f_ref --2bitdir $DIR_db/blat -t $f_ref_2bit -o $d33/03_predSV \\
            --blatserver elmob111 --blatport 1986 -r $chr", 1);
    }
    
    my $d33_03 = "$d33/03_predSV";
    my $f33_11 = "$d33/11_sum.tbl";
    open(FHO, ">$f33_11") or die "cannot write to $f33_11\n";
    print FHO join("\t", qw/acc chr_l pos_l strand_l n_sc_l chr_r pos_r strand_r n_sc_r type cov_l cov_r len_l len_r pct_idty_l pct_reads_nu_l pct_idty_r pct_reads_nu_r beg_con chr_b beg_con_map end_con chr_e end_con_map seq beg_con2 chr_b2 beg_con_map2 end_con2 chr_e2 end_con_map2 seq2/)."\n";
    
    opendir(my $dh, $d33_03) || die "cannot open $d33_03\n";
    for my $fn (sort readdir $dh) {
        if($fn =~ /^(HM\d+)\.bam\.predSV/) {
            my $acc = $1;
            my $fi = "$d33_03/$fn";
            open(FHI, "<$fi") or die "cannot read $fi\n";
            while(<FHI>) {
                chomp;
                my @ps = split("\t", $_);
                if(@ps == 24) {
                    push @ps, ("") x 7;
                } elsif(@ps == 31) {
                } else {
                    die "cannot handle ".@ps." columns\n";
                }
                die join("\t", @ps) if @ps ne 31;
                print FHO join("\t", $acc, @ps)."\n";
            }
            close FHI;
        }
    }
    closedir $dh;
    close FHO;
}
sub run_tigrasv {
    my ($dir, $sms) = @_;
    my $d01 = "$dir/01_pos_sorted"; 
    my $d02 = "$dir/02_rn_sorted"; 
    my $f03 = "$dir/03_stat.tbl"; 
    my $t = readTable(-in=>$f03, -header=>1);
    
    my $d37 = "$dir/37_tigrasv"; 
    print "tigra-sv -R $f_ref $d37/02_sv.tbl -f $d37/01_bamlist.tbl -o $d37/03.fa -r $d37/04_ref.fa\n";
}
sub run_genome_strip {
    my ($dir, $sms) = @_;
    my $d01 = "$dir/01_pos_sorted"; 
    my $d02 = "$dir/02_rn_sorted"; 
    my $f03 = "$dir/03_stat.tbl"; 
    my $t = readTable(-in=>$f03, -header=>1);
    
    my $f_ref_mask = "$DIR_db/bwa/mt_40.mask.fasta";
    my $d35 = "$dir/35_genome_strip"; 
    
    my $chr = "chr5";
    my $lib = "$svtoolkit/lib";
    my $qscript = "$svtoolkit/qscript";
    for my $sm (@$sms) {
        runCmd("java -Xmx3g -cp $lib/gatk/Queue.jar:$lib/SVToolkit.jar:$lib/gatk/GenomeAnalysisTK.jar \\
            org.broadinstitute.sting.queue.QCommandLine \\
            -S $qscript/SVPreprocess.q -S $qscript/SVQScript.q \\
            -gatk $lib/gatk/GenomeAnalysisTK.jar \\
            -cp $lib/SVToolkit.jar:$lib/gatk/GenomeAnalysisTK.jar \\
            -configFile \$conf/svtoolkit/mt.txt -tempDir $DIR_tmp \\
            -md $d35/01_preprocess -R $f_ref -genomeMaskFile $f_ref_mask \\
            -I $d01/$sm.bam -run", 1);
        runCmd("java -Xmx3g -cp $lib/gatk/Queue.jar:$lib/SVToolkit.jar:$lib/gatk/GenomeAnalysisTK.jar \\
            org.broadinstitute.sting.queue.QCommandLine \\
            -S $qscript/SVDiscovery.q -S $qscript/SVQScript.q \\
            -gatk $lib/gatk/GenomeAnalysisTK.jar \\
            -cp $lib/SVToolkit.jar:$lib/gatk/GenomeAnalysisTK.jar \\
            -configFile \$conf/svtoolkit/mt.txt -tempDir $DIR_tmp \\
            -runDirectory $d35 \\
            -md $d35/01_preprocess -R $f_ref -genomeMaskFile $f_ref_mask \\
            -I $d01/$sm.bam -O $d35/02_alt/$sm.vcf \\
            -minimumSize 100 -maximumSize 1000000 \\
            -windowSize 10000000 -windowPadding 10000 -run", 1);
        runCmd("java -Xmx3g -Djava.library.path=$svtoolkit/bwa \\
            -cp $lib/gatk/Queue.jar:$lib/SVToolkit.jar:$lib/gatk/GenomeAnalysisTK.jar \\
            org.broadinstitute.sting.queue.QCommandLine \\
            -S $qscript/SVAltAlign.q -S $qscript/SVQScript.q \\
            -gatk $lib/gatk/GenomeAnalysisTK.jar \\
            -cp $lib/SVToolkit.jar:$lib/gatk/GenomeAnalysisTK.jar \\
            -configFile \$conf/svtoolkit/mt.txt -tempDir $DIR_tmp \\
            -md $d35/01_preprocess -runDirectory $d35 -R $f_ref \\
            -vcf $d35/02_alt/$sm.vcf \\
            -I $d01/$sm.bam -O $d35/03_alt_align/$sm.bam -run", 1);
    }
}




