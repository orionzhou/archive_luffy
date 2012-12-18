#!/usr/bin/perl -w
use strict;
use lib ($ENV{"m"}."/mt2/modules");
use InitPath;
use Common;
use Bio::SeqIO;
use Gff;
use Seq;
use Data::Dumper;
use List::Util qw/min max sum/; use POSIX qw/ceil floor/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
sub sra2Fq {
    my ($fi, $fo) = @_;
    my $cmd = "awk '{  if(NR % 4 == 1) {
                                        print \"@\" \$2 > \"$fo\"
                                    } else if(NR % 4 == 3) {
                                        print \"+\" \$2 > \"$fo\"
                                    } else if(NF == 1) {
                                        print > \"$fo\"
                                    } else {
                                        print \"error format:\" \$0 > \"/dev/stderr\"
                                        exit 1
                                    }
                                } END {print NR/4 \" reads processed\"}' $fi";
    system($cmd);
}
sub download_SRA {
    my ($d02, $url1, $url2, $url3) = @_; 
    for my $j (1..2) {
        my $url = "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/static/$url1/$url2/$url3\_$j.fastq.bz2";
        my $pre = "$url2\_$url3\_$j";
        my $f1 = file($d02, "$pre\_orig.fq.bz2");
        system("wget $url -O $f1");
        my $f2 = file($d02, "$pre\_orig.fq");
        system("bunzip2 $f1");
        my $f3 = file($d02, "$pre.fq");
        sra2Fq($f2, $f3);
        my $f4 = file($d02, "$pre.fq.gz");
        system("gzip $f3");
        system("rm $f2");
    }
}
sub getRInput {
#getRInput(-outdir=>$d01, -fwin=>$f02, -param=>$p);
    my ($outDir, $f, $param) = rearrange(['outdir', 'fwin', 'param'], @_);
    my $inDir = dir('/export/lab/analysis/antoine/WindowsRhoData/LDHatFiles100kbwin/');
    my $refDb = $param->{refdb};
    my $fh = new IO::File $f, "w";
    print $fh join("\t", qw/round chr wStart wEnd cntSnp sStart sEnd/)."\n";
    for my $i (1..8) {
        my $chr = "chr$i";
        my $wh = getWindows(-chr=>$i, -db=>$refDb, -winsize=>100_000, -winstep=>100_000, -opt=>2);
        my $cnt = 0;
        for my $j (0..@$wh-1) {
            my ($wS, $wE) = @{$wh->[$j]};
            my $round = sprintf("%04d", $j+1);
            my $round2 = sprintf("%01d_%03d", $i, $j+1);
            my $f1 = file($inDir, "$chr\_win_$round.locs");
            my ($cntSnp, $sS, $sE) = (0, "", "");
            if(-s $f1) {
                my $f2 = file($inDir, "$chr\_win_$round.sites");
                die "$f2 is not there\n" unless -s $f2;
                my $fh1 = new IO::File $f1, "r";
                my $line = readline($fh1);
                $line = readline($fh1);
                my @ps = split("\t", $line);
                $cntSnp = @ps;
                $sS = $wS + $ps[0] - 1;
                $sE = $wS + $ps[-1] - 1;
                die "snpEnd[$sE] > winEnd[$wE] at $f1\n" if $sE > $wE;
                my $fo1 = file($outDir, "$round2\_ldhat_loc.txt");
                my $fo2 = file($outDir, "$round2\_ldhat.txt");
                system("cp $f1 $fo1");
                system("cp $f2 $fo2");
            }
            $cnt += $cntSnp;
            print $fh join("\t", $round2, $chr, $wS, $wE, $cntSnp, $sS, $sE)."\n";
        }
        print join("\t", $chr, $cnt)."\n";
    }
}
sub fas2GecoTxt {
#my $refDb = "mt_30";
#my $fi = file($DIR_Genome, $refDb, "41_genome.fa");
#my $dirO = dir($DIR_Data, "db/geco", $refDb);
#fasta2Txt($fi, $dirO);
    my ($fi, $dirO) = @_;
    system("mkdir -p $dirO") unless -d $dirO;
    my $seqH = Bio::SeqIO->new(-file=>$fi, -fasta=>"fasta");
    while(my $seq = $seqH->next_seq() ) {
        my $fo = file($dirO, $seq->id);
        my $fh = new IO::File $fo, "w";
        print $fh $seq->seq."\n";
        close $fh;
    }
}
sub goldengate_MergeSnps {
=cut
my $dir = "/export/lab/data/mt35_illumina/variants300/snps/GoldenGate/illumina";
my @names = qw/Gentz Joelle Marie Monteros random/;
my $f10 = "$dir/10_snps.tbl";
goldengate_MergeSnps($dir, \@names, $f10);
my $f_ref = "$DIR_Genome/mt_35/41_genome.fa";
my $f20 = "$dir/20_illumna.csv";
goldengate_Convert2Illumina($f10, $f_ref, $f20);
=cut
    my ($dir, $names, $fo) = @_;

    my $h;
    for my $name (@$names) {
        my $fi = "$dir/snps_$name.txt";
        open(FHI, "<$fi");
        my $line = <FHI>;
        
        while(<FHI>) {
            chomp;
            my @ps = split "\t";
            my ($locus, $chr, $pos, $ref, $maf, $cntA, $cntC, $cntG, $cntT, $n) = @ps[0..2,4,5,7..11];
            $chr = "chr$chr" if $chr =~ /^\d+$/;
            
            my $r = {'A'=>$cntA, 'C'=>$cntC, 'G'=>$cntG, 'T'=>$cntT};
            my @alleles = grep {$r->{$_} > 0} keys(%$r);
            die "not biallelic: $chr:$pos\n" unless @alleles == 2;
            my @vars = grep {$_ ne $ref} @alleles;
            die "not biallelic: $name $chr:$pos\n" unless @vars == 1;
            my $var = $vars[0];

            my $stat = [$ref, $var, $maf, $n, $locus, $name];
            $h->{$chr}->{$pos} ||= [];
            push @{$h->{$chr}->{$pos}}, $stat;
        }
        close FHI;
    }

    open(FHO, ">$fo");
    print FHO join("\t", qw/chr pos ref var locus source/)."\n";
    for my $chr (sort(keys(%$h))) {
        my $h2 = $h->{$chr};
        for my $pos (sort {$a<=>$b} keys(%$h2)) {
            my @stats = @{$h2->{$pos}};
            my @refs = map {$_->[0]} @stats;
            die Dumper(@stats) unless uniq(@refs) == 1;
            my @vars = map {$_->[1]} @stats;
            die Dumper(@stats) unless uniq(@vars) == 1;
            my @mafs = map {$_->[2]} @stats;
            die Dumper(@stats) unless uniq(@mafs) == 1;
            my @ns = map {$_->[3]} @stats;
            die Dumper(@stats) unless uniq(@ns) == 1;
            my $locus = join(" ", map {$_->[-2]} @stats);
            my $source = join(" ", map {$_->[-1]} @stats);
            print FHO join("\t", $chr, $pos, $refs[0], $vars[0], $mafs[0], $ns[0], $locus, $source)."\n";
        }
    }
    close FHO;
}
sub goldengate_Convert2Illumina {
=cut
my $dir = "$DIR_Misc3/Mt_goldengate";
my $fi = "$dir/20120718_in.tbl";
my $fo = "$dir/20120718_illumina.csv";
my $f_ref = "$DIR_Data/genome/mt_35/41_genome.fa";
goldengate_Convert2Illumina($fi, $f_ref, $fo);
=cut
    my ($fi, $f_ref, $fo) = @_;
    my $sep = ",";
    open(FHO, ">$fo") || die "cannot open $fi for reading\n";
    print FHO join($sep, qw/Locus_Name Sequence Target_Type Genome_Build_Version Chromosome Coordinate Source Source_Version Sequence_Orientation Plus_Minus/)."\n";

    open(FHI, "<$fi");
    my $line = <FHI>;
    while(<FHI>) {
        chomp;
        my @ps = split "\t";
        my ($chr, $pos, $ref, $var, $locus, $src) = @ps;
        $chr = "chr$chr" if $chr =~ /^\d+$/;
        my $loc = [[$pos-60, $pos+60]];
        my $seq = seqRet($loc, $chr, "+", $f_ref);

        my $seq_up = substr($seq, 0, 60);
        my $ref2 = substr($seq, 60, 1);
        my $seq_dw = substr($seq, 61, 60);
        die "not ref: $ref != $ref2\n" unless $ref eq $ref2;
        
        my $seq_str = "$seq_up\[$ref/$var\]$seq_dw";
        print FHO join($sep, $locus, $seq_str, "SNP", "Mt3.5", $chr, $pos, $src, 0, "Forward", "Plus")."\n";
    }
    close FHI;
    close FHO;
}
sub change_zm_fasta_id {
#my $dir = "/project/youngn/zhoup/Data/misc3/spada/Zmays/01_genome";
#change_zm_fasta_id("$dir/ZmB73_RefGen_v2.masked.fasta", "$dir/01_refseq.fa");
    my ($fi, $fo) = @_;
    open(FHI, "<$fi");
    open(FHO, ">$fo");
    while(<FHI>) {
        chomp;
        if(/^\>(.+)$/) {
            if($1 =~ /chromosome (\w+)$/) {
                my $id = $1;
                if($id =~ /^\d+$/) {
                    $id = "chr$id";
                } elsif($id eq "UNKNOWN") {
                    $id = "chrU";
                } elsif($id eq "mitochondrion") {
                    $id = "Mt";
                } elsif($id eq "chloroplast") {
                    $id = "Pt"
                } else {
                    die "unknonw chr ID: $id\n";
                }
                print FHO ">$id\n";
                print "  $id\n";
            } else {
                die "unknown chr ID: $1\n";
            }
        } else {
            print FHO $_."\n";
        }
    }
    close FHI;
    close FHO;
}
sub get_fas_gaps {
    my ($fi, $fo) = @_;
    my $seqI = Bio::SeqIO->new(-file=>"<$fi", -format=>"fasta");

    open(FHO, ">$fo") or die "cannot open $fo for writing\n";
    print FHO join("\t", qw/id beg end len type/)."\n";
    while(my $seqObj = $seqI->next_seq()) {
        my ($id, $seq, $seqlen) = ($seqObj->id, $seqObj->seq, $seqObj->length);
        print FHO join("\t", $id, 1, $seqlen, $seqlen, 'chr')."\n";
        while($seq =~ /N+/ig) {
            my ($beg, $end) = ($-[0]+1, $+[0]);
            my $len = $end - $beg + 1;
            print FHO join("\t", $id, $beg, $end, $len, 'gap')."\n";
        }
    }
}
my $dir = "$DIR_Genome/Mtruncatula_4.0";
my $f01 = "$dir/01_refseq.fa";
my $f51 = "$dir/51_gap_loc.tbl";
get_fas_gaps($f01, $f51);


