package Blast;
use strict; 
use Common;
use Bio::SeqIO;
use File::Basename;
use Data::Dumper;
use File::Path qw/make_path remove_tree/;
use List::Util qw/min max sum/;
use POSIX qw/ceil floor/;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;
use URI::Escape;
use LWP::UserAgent;
use HTTP::Request::Common qw/POST/;

use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw/split_blast_output blast_tiling blast_nr blast_nr_batch/;
@EXPORT_OK = qw//;

sub write_blast {
    my ($id, $rows, $fo) = @_;
    open(my $fho, ">$fo") or die "cannot open $fo for writing\n";
    print $fho join("\t", qw/qId qBeg qEnd strand hId hBeg hEnd pct e score/)."\n";
    print $fho join("\n", map {join("\t", @$_)} @$rows)."\n";
    close $fho;
}
sub split_blast_output {
    my ($fi, $dirO) = @_;
    make_path($dirO) unless -d $dirO;
    
    open(my $fhi, "<$fi") or die "cannot open $fi for reading\n";
    my ($idP, $lineP) = ("", "");
    my @rows;
    while(my $line = <$fhi>) {
        chop($line);
        my ($qId, $hId, $pct, $alnLen, $mm, $gap, $qBeg, $qEnd, $hBeg, $hEnd, $e, $score) = split("\t", $line);
        my $strand = $hBeg > $hEnd ? "-" : "+";
        ($hBeg, $hEnd) = ($hEnd, $hBeg) if $strand eq "-";
        my $row = [$qId, $qBeg, $qEnd, $strand, $hId, $hBeg, $hEnd, $pct, $e, $score];
        if($qId eq $idP) {
            push @rows, $row;
        } else {
            write_blast($idP, \@rows, "$dirO/$idP.tbl") if $idP;
            @rows = ($row);
            $idP = $qId;
        }
    }
    write_blast($idP, \@rows, "$dirO/$idP.tbl") if @rows > 0;
}


sub blast_tiling {
    my ($fl, $dirI, $fo) = @_;
    open(my $fho, ">$fo") or die "cannot open $fo for writing\n";
    print $fho join("\t", qw/qId qBeg qEnd strand hId hBeg hEnd qLen hLen pct e score/)."\n";
   
    my $tl = readTable(-in=>$fl, -header=>1);
    my @ids = $tl->col("id");
    my ($cntB, $cntG) = (0, 0);
    for my $id (sort @ids) {
        my $fi = "$dirI/$id.tbl";
        if( ! -s $fi ) {
            $cntB ++;
            next;
        }
       
        my $ti = readTable(-in=>$fi, -header=>1);
        my @locs = map { [$ti->elm($_, "qBeg"), $ti->elm($_, "qEnd")] } (0..$ti->nofRow-1);
        my @es = $ti->col("e");
        my $refs = tiling(\@locs, \@es, 1);
        for (@$refs) {
            my ($beg, $end, $idx) = @$_;
            next if ($end - $beg + 1) < 100;
            my ($qId, $qBeg, $qEnd, $strd, $hId, $hBeg, $hEnd, $pct, $e, $score) = $ti->row($idx);
            my $begH = sprintf "%d", $hBeg + ($beg-$qBeg) * ($hEnd-$hBeg)/($qEnd-$qBeg);
            my $endH = sprintf "%d", $hEnd - ($qEnd-$end) * ($hEnd-$hBeg)/($qEnd-$qBeg);
            my $qLen = $end - $beg + 1;
            my $hLen = $endH - $begH + 1;
            print $fho join("\t", $id, $beg, $end, $strd, $hId, $begH, $endH, $qLen, $hLen, $pct, $e, $score)."\n";
        }
        printf " %4d: %10d\n", (++$cntG)+$cntB, $ti->nofRow;
    }
    close $fho;
    printf "\n %4d / %4d tiled\n", $cntG, $cntG+$cntB;
}
sub blast_nr {
    my ($id, $seq, $fo) = @_;
    my $ua = LWP::UserAgent->new;

    my $qry = uri_escape(">$id\n$seq\n");

    my ($program, $db) = ("blastn", "nr");
    my $args = "CMD=Put&PROGRAM=$program&DATABASE=$db&QUERY=".$qry;
    my $req = new HTTP::Request POST => 'http://www.ncbi.nlm.nih.gov/blast/Blast.cgi';
    $req->content_type('application/x-www-form-urlencoded');
    $req->content($args);

    my $response = $ua->request($req);
    $response->content =~ /^    RID = (.*$)/m;
    my $rid=$1;
    $response->content =~ /^    RTOE = (.*$)/m;
    my $rtoe=$1;
    printf "  %s: RID [%s] [%5ds]\n", $id, $rid, $rtoe;
    sleep $rtoe;

    while(1) {
        sleep 3;

        $req = new HTTP::Request GET => "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=$rid";
        $response = $ua->request($req);

        if ($response->content =~ /\s+Status=WAITING/m) {
            # print STDERR "Searching...\n";
            next;
        }
        if ($response->content =~ /\s+Status=FAILED/m) {
            print STDERR "Search $rid failed; please report to blast-help\@ncbi.nlm.nih.gov.\n";
            exit 4;
        }
        if ($response->content =~ /\s+Status=UNKNOWN/m) {
            print STDERR "Search $rid expired.\n";
            exit 3;
        }
        if ($response->content =~ /\s+Status=READY/m) {
            if ($response->content =~ /\s+ThereAreHits=yes/m) {
                #  print STDERR "Search complete, retrieving results...\n";
                last;
            } else {
                print STDERR "No hits found.\n";
                exit 2;
            }
        }

        # if we get here, something unexpected happened.
        exit 5;
    }

    $req = new HTTP::Request GET => "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_TYPE=Text&RID=$rid";
    $response = $ua->request($req);

    open (FHO, ">$fo") || die "Can't open file $fo for writing: $!\n";
    print FHO $response->content;
    close FHO;
}
sub blast_nr_batch {
    my ($fi, $do) = @_;
    make_path($do) unless -d $do;

    my $seqH = Bio::SeqIO->new(-file=>"<$fi", -format=>'fasta');
    while(my $seqO = $seqH->next_seq()) {
        my ($id, $seq) = ($seqO->id, $seqO->seq);
        my $fo = "$do/$id.txt";
        if( -s $fo ) {
            print "  $id: exists - skipped\n";
            next;
        }
        blast_nr($id, $seq, $fo);
    }
    $seqH->close();
}

1;
__END__
