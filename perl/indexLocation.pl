#!/usr/bin/perl
use strict; use Init; use Common; use Localdb;
use Bio::Seq; use Path::Class; use Qry; use Data::Dumper;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
#store the file position of gene starts in db table, so that retrieve variants
#in a particular gene only need to look up the index table and jump to that file position
sub getFilePos {
    my ($locAry, $chr, $qry) = @_;
    my $fIn = file($DIR_Variant, "$chr\_C.txt");
    my $fInH = new IO::File $fIn, "r";
    my $locPos = {};
    my ($chrPos_prev, $filePos_prev) = (0, 0);
    for my $loc (sort {$a<=>$b} @$locAry) {
        while(my $line = readline($fInH)) {
            next if $line =~ /^Reference/i;
            my @eleAry = split("\t", $line);
            die("Invalid position: ".$eleAry[3]."\n") unless $eleAry[3] =~ /^[\d\.]+$/;
            my ($chrPos, $filePos) = ($eleAry[3], tell($fInH));
            if($loc > $chrPos_prev && $loc <= $chrPos) {
                die("$loc exists\n") if exists $locPos->{$loc};
                $locPos->{$loc} = $filePos_prev;
                backOneLine($fInH);
                last;
            }
            ($chrPos_prev, $filePos_prev) = ($chrPos, $filePos);
        }
    }
    #for my $gene (sort keys %$geneLoc) {
    #  print join("\t", $gene, $geneLoc->{$gene}, $locPos->{$geneLoc->{$gene}})."\n";
    #}
    $qry->storeLocIdx($chr, $locPos);
}
sub indexLocation {
    my ($dbs, $chrAry, $types) = rearrange(['db', 'chr', 'types'], @_);
    my $qry = Qry->new(-table=>'geneidx', -cols1=>['chr', 'loc'], -cols2=>'pos1');
    for my $chr (@$chrAry) {
        my @loci;
        my $i = 0;
        for my $db (@$dbs) {
            my $ld = Localdb->new(-db=>$db);
            my @fe = $ld->{biodb}->features(-seq_id=>$chr, -types=>$types);
            print join("\t", $db, $chr, scalar(@fe))."\n";
            for my $fe (sort {$a->start <=> $b->start} @fe) {
                push @loci, $fe->start;
                last if ++$i > 1000000;
            }
        }
        @loci = uniq @loci;
        getFilePos(\@loci, $chr, $qry);
    }
}
my @dbs = qw/mt_30 mt_defl/;
my @chrAry = map {"chr".$_} (0);
my @types  = qw/mRNA tRNA rRNA/;
indexLocation(-db=>\@dbs, -chr=>\@chrAry, -types=>\@types);

