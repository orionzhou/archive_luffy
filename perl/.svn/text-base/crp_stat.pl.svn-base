#!/usr/bin/perl -w
use strict;
use Common;
use Seq;
use List::Util qw/min max sum/;

$ENV{'PATH'} = join(":", "/project/youngn/zhoup/Source/signalp-3.0",
    "/home/msi/zhoup/bin/x86_64/bin",
    $ENV{'PATH'});

my $dir = "/project/youngn/zhoup/Data/misc2/crp_misc";
my $d10 = "$dir/10_sigp";
my $f10_01 = "$d10/01.gtb";
my $f10_06 = "$d10/06.fa";
my $f10_11 = "$d10/11_sigp.tbl";
#crp_stat_sigp($f10_01, $f10_06, $f10_11);

sub crp_stat_sigp {
    my ($f_gtb, $f_seq, $fo) = @_;
    my $tg = readTable(-in=>$f_gtb, -header=>1);
    my $h_seq = readSeq($f_seq, 2);
  
    open(FH, ">$fo");
    print FH join("\t", qw/id parent cat2 cat3 note sigp sp_pos/)."\n";
    for my $i (0..$tg->nofRow-1) {
        my ($id, $pa, $chr, $beg, $end, $strand, $locE, $locI, $locC, $loc5, $loc3, $phase, $source, $conf, $cat1, $cat2, $cat3, $note) = $tg->row($i);
        my $seq_pro = $h_seq->{$pa}->seq;
        die "no sequence for $id\n" unless $seq_pro;
        $seq_pro =~ s/\*$//;
        my $seq = Bio::Seq->new(-id=>$id, -seq=>$seq_pro);
        my ($sp, $sp_pos) = run_signalp($seq, 3);
        print FH join("\t", $id, $pa, $cat2, $cat3, $note, $sp, $sp_pos)."\n";
    }
    close FH;
}
sub run_signalp {
    my ($seq, $version) = @_;
    my $fTmp = "sigp.fa";
    writeSeq($seq, $fTmp);
    my ($prob1, $prob2, $site) = ("") x 3;
    open(OUT, "signalp -t euk $fTmp |") or die "failed: $!\n";
    my ($pos, $score, $sp) = ("", "", 0);
    if($version == 3) {
        my ($score1, $score2);
        while( <OUT> ) {
            chomp;
            if(/Signal peptide probability: ([\d\.]+)/i) {
                $score1 = $1;
            } elsif(/Signal anchor probability: ([\d\.]+)/i) {
                $score2 = $1;
            } elsif(/Max cleavage site probability: [\d\.]+ between pos\. (\d+)/i) {
                $pos = $1+1;
            }
        }
        $score = max($score1, $score2);
        $sp = ($score >= 0.450) ? 1 : 0;
        $pos = "" if $sp == 0;
    } elsif($version == 4) {
        while( <OUT> ) {
            chomp;
            next unless /^tmp/;
            my ($id, $Cmax, $posC, $Ymax, $posY, $Smax, $posS, $Smean, $D, $tag, $Dmaxcut, $network) = split " ";
            ($pos, $score, $sp) = ($posY, $D, 1) if $tag eq "Y";
        }
    } else {
        die "unsupported version: $version\n";
    }
    system("rm $fTmp");
    return ($sp, $pos);
}



