#!/usr/bin/perl -w
use strict;
use FindBin;
use lib $FindBin::Bin;
use InitPath;
use Common;
use Bio::SeqIO;
use Gff;
use Seq;
use Data::Dumper;
use List::Util qw/min max sum/; use POSIX qw/ceil floor/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
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
sub tair9_updates {
#my $dir = "$DIR_Misc2/crp.ssp";
#tair9_updates("$dir/tair9_updates_raw.tbl", "$dir/tair9_updates.tbl");
    my ($fi, $fo) = @_;
    open(FHO, ">$fo") or die "cannot open $fo for writing\n";
    print FHO join("\t", qw/chr pos1 pos2 nt_old nt_new/)."\n";
    
    my $hb;
    my $t = readTable(-in=>$fi, -header=>1);
    for my $i (1..$t->nofRow) {
        my ($chr, $version, $pos1, $pos2, $type, $nt, $nt_new) = $t->row($i-1);
        if(!defined($hb->{$chr})) {
            print FHO join("\t", $chr, 1, 1, '', '')."\n";
            $hb->{$chr} = 1;
        }

        my $len;
        if($nt =~ /^Nx(\d+)$/) {
            $len = $1;
        } else {
            $len = length($nt);
        }

        if($type eq "substitution") {
            print FHO join("\t", $chr, $pos1, $pos2, $nt, $nt_new)."\n";
        } elsif($type eq "insertion") {
            print FHO join("\t", $chr, $pos1, $pos2, '', '')."\n";
            print FHO join("\t", $chr, $pos1+1, $pos2+1+$len, '', '')."\n";
        } elsif($type eq "deletion") {
            print FHO join("\t", $chr, $pos1, $pos2, '', '')."\n";
            print FHO join("\t", $chr, $pos1+1+$len, $pos2+1, '', '')."\n";
        } else {
            die "unknown type: $type\n";
        }
    }
    close FHO;
}



