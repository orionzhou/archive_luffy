#!/usr/bin/perl -w
use strict;
use lib ($ENV{"SCRIPT_HOME_PERL"});
use InitPath;
use Common;
use Gff;
use Gtb;
use Seq;
use Align;
use Data::Dumper;
use Time::HiRes qw/gettimeofday tv_interval/;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
sub test_time {
    my $t0 = [gettimeofday];
    my $t1 = [gettimeofday];
    printf "***%.01f min\n", tv_interval($t0, $t1) / 60;
}
sub gmSnp {
    my $DIR_Work = dir($DIR_Misc2, 'gm_snp2');
    my $paramR = [
        {fname=>'01_LG_G_1_10.txt', head=>1, left=>1},
        {fname=>'01_LG_G_1_10_noSSRs.txt', head=>3, left=>1},
        {fname=>'summer07_F.txt', head=>3, left=>1},
        {fname=>'summer07_G.txt', head=>3, left=>1},
        {fname=>'summer07_A2.txt', head=>3, left=>1},
        {fname=>'summer07_J.txt', head=>3, left=>1},
        {fname=>'summer07_O.txt', head=>3, left=>1},
        {fname=>'01_LG_G_oct08.txt', head=>2, left=>1} ];
    for my $param (@$paramR[2]) {
        my $fIn = file($DIR_Work, $param->{fname});
        #writeGmSnp($fIn, $param);
    }
}
sub testDeepMerge {
    my $a = [ [[6,15]], [[25,35]], [[1,10], [21,30]], [[50,60], [70,80]], [[60,64]], [[90,100]] ];
    my $b = posMergeDeep($a);
    for (@$b) {
        my ($locA, $idxs) = @$_;
        my $locStr = join(" ", map {join("_", @$_)} @$locA);
        print join("\t", $locStr, join(",", @$idxs))."\n";
    }
}
sub test_pg_connect {
    use DBI;
    my $pg_db = 'chado';
    my ($pg_host, $pg_user, $pg_pw) = ($ENV{'PGSQLH'}, $ENV{'PGSQLU'}, $ENV{'PGSQLP'});
    my $dbh = DBI->connect("dbi:Pg:dbname=$pg_db;host=$pg_host", $pg_user, $pg_pw)
        or die "cannot connect to pgsql\n";
}
sub convId2Url {
    my $dir = dir($DIR_In, "seminar_March_7");
    my @fns = qw/frog_a bovine_a human_a mouse_a frog_b bovine_b human_b mouse_b/;
    my $seqNH = Bio::SeqIO->new(-file=>">$dir/hemoglobin_na.fa", -format=>"fasta");
    my $seqAH = Bio::SeqIO->new(-file=>">$dir/hemoglobin_aa.fa", -format=>"fasta");
    for my $fn (@fns) {
        my $f_gb = file($dir, "01_seqs/$fn.gb");
        my $seqH = Bio::SeqIO->new(-file=>$f_gb, -format=>"genbank");
        my $seq = $seqH->next_seq();
        print join("\t", $seq->id, $seq->length)."\n";
        my @fes_cds = grep {$_->primary_tag eq "CDS"} $seq->get_SeqFeatures();
        die " not 1 CDS in $fn\n" unless @fes_cds == 1;
        my $fe = $fes_cds[0];
        my $seqStr = $seq->subseq($fe->start, $fe->end);
        my $seqCds = Bio::Seq->new(-id=>$fn, -seq=>$seqStr);
        $seqCds = $seqCds->revcom() if $fe->strand == -1;
        $seqNH->write_seq($seqCds);
        print "\t".join("\t", $seqCds->id, $seqCds->length)."\n";
        my $seqPro = $seqCds->translate();
        $seqAH->write_seq($seqPro);
        print "\t".join("\t", $seqPro->id, $seqPro->length)."\n";
    }
}
sub cp_fq_scratch {
#cp_fq_scratch();
    my $fi = "$DIR_Misc3/hapmap/09_fastq.tbl";
    my $do = "/project/scratch/zhoup/data_for_taehyun/fastq";
    my $t = readTable(-in=>$fi, -header=>1);
    my %h_sm = map { $_=>1 } qw/HM004 HM010 HM050 HM056 HM101/;
    for my $i (0..$t->nofRow-1) {
        my ($idx, $sm, $rg, $lb, $id_lb, $id_run, $pl, $rl, $pi, $f1a, $f1b) = $t->row($i);
        next unless exists $h_sm{$sm};
        my $fo1 = "$do/$rg.1.fq.gz";
        my $fo2 = "$do/$rg.2.fq.gz";
        runCmd("cp $f1a $fo1");
        runCmd("cp $f1b $fo2");
    }
}
sub fasta_mask {
#my $fi = "/project/youngn/zhoup/Data/misc3/spada/Zmays/01_genome/01_refseq.fa";
#my $fo = "/project/youngn/zhoup/Data/misc3/spada/Zmays_masked/01_genome/01_refseq.fa";
    my ($fi, $fo) = @_;
    open(FHI, "<$fi");
    open(FHO, ">$fo");
    while(<FHI>) {
        chomp;
        if(/^\>/) {
            print FHO $_."\n";
            print $_."\n";
        } else {
            $_ =~ s/[atcgn]/N/g;
            print FHO $_."\n";
        }
    }
    close FHI;
    close FHO;
}
sub double_genome {
#my $fi = "/project/youngn/zhoup/Data/misc3/spada/Athaliana/01_genome/01_refseq.fa";
#my $fo = "/project/youngn/zhoup/Data/misc3/spada/Athaliana_large/01_genome/01_refseq.fa";
#double_genome($fi, $fo);
    my ($fi, $fo) = @_;
    my $seqHI = Bio::SeqIO->new(-file=>"<$fi", -format=>'fasta');
    my $seqHO = Bio::SeqIO->new(-file=>">$fo", -format=>"fasta");
    while(my $seq = $seqHI->next_seq()) {
        print join("\t", $seq->id, $seq->length)."\n";
        $seqHO->write_seq($seq);
        my $seqN = Bio::Seq->new(-id=>$seq->id."_N", -seq=>"N"x$seq->length);
        $seqHO->write_seq($seqN);
    }
    $seqHI->close();
    $seqHO->close();
}



