package Bed;
use strict; 
use Data::Dumper;
use File::Basename;
use Common; 
use Gff;
use List::Util qw/min max sum/; use POSIX qw/ceil floor/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw/bed_refine bedpe2Gff/;
@EXPORT_OK = qw//;

sub bedpe2Gff {
    my ($fIn, $fOut, $acc) = rearrange(['in', 'out', 'acc'], @_);
    my $fInH = new IO::File $fIn, "r";
    my $fOutH = new IO::File $fOut, "w";
    printGffHead($fOutH);
    while( <$fInH> ) {
        chomp;
        next unless $_;
        my ($c1, $s1, $e1, $c2, $s2, $e2, $id, $numDistinctPairs, $strand1, $strand2,
            $meanEditDist1, $meanEditDist2, $meanMappings1, $meanMappings2, $size,
            $numMappings, $allWeightedSupport, $finalSupport, $finalWeightedSupport,
            $numUniqPairs, $numAnchoredPairs, $numMultiplyMappedPairs) = split("\t");
        die "not on same chr: $id [ $c1, $c2 ]\n" unless $c1 eq $c2;
        $c1 = "chr$c1";
        $c2 = "chr$c2";
        my $note = "size:".pretty($size)." editDist:$meanEditDist1-$meanEditDist2"
            . " numDistinctPairs:$numDistinctPairs";
        my $tagH = {ID=>$id, Name=>"$acc\_$id", Note=>$note};
        my $ref = {
            seqid=>$c1, source=>'hydra', type=>"read_pair", start=>min($s1, $s2), 
            end=>max($e1, $e2), score=>$finalWeightedSupport, tags=>$tagH
        };
#    print $fOutH fe2GffLine($ref);
        my $ref1 = {seqid=>$c1, type=>"read", start=>$s1, end=>$e1, strand=>$strand1, tags=>{Parent=>$id}};
        my $ref2 = {seqid=>$c2, type=>"read", start=>$s2, end=>$e2, strand=>$strand2, tags=>{Parent=>$id}};
        $ref1->{tags}->{Note} = $note;
        $ref2->{tags}->{Note} = $note;
        print $fOutH fe2GffLine($ref1);
        print $fOutH fe2GffLine($ref2);
    }
}
sub bed_refine {
    my ($fi, $fo) = @_;
    my $ft = dirname($fo)."/tmp.bed";
    runCmd("mergeBed -i $fi -nms > $ft", 1);
    
    open(FH1, "<$fi");
    my $h;
    while(<FH1>) {
        chomp;
        next if /^\#/;
        my ($chr, $beg, $end, $id, $score, $strand, $note) = split("\t");
        $h->{$id} = [$chr, $beg, $end, $score, $strand, $note];
    }
    close FH1;

    open(FH2, "<$ft");
    open(FH3, ">$fo");
    print FH3 "#track name=gene_models useScore=0\n";
    while(<FH2>) {
        chomp;
        next if /^\#/;
        my ($chr, $beg, $end, $id_str) = split("\t");
        my @ids = split(";", $id_str);
        if(@ids > 1) {
            @ids = sort {$h->{$a}->[2]-$h->{$a}->[1] <=> $h->{$b}->[2]-$h->{$b}->[1]} @ids;
        }
        my $id = $ids[-1];
        die "cannot find info of $id\n" unless exists $h->{$id};
        my @stats = @{$h->{$id}};
        print FH3 join("\t", @stats[0..2], $id, @stats[3..5])."\n";
    }
    close FH2;
    close FH3;
    system("rm $ft");
}


1;
__END__
