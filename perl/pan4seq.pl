#!/usr/bin/perl -w
use strict;
use FindBin;
use lib "$FindBin::Bin";
use Common;
use Bio::SeqIO;
use Bio::DB::Fasta;
use Data::Dumper;
use File::Path qw/make_path remove_tree/;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;

my $dir = "/home/youngn/zhoup/Data/misc3/pan4seq";
chdir $dir || die "cannot chdir to $dir\n";

my @orgs = qw/HM340.APECCA HM034 HM056/;

#merge_seq(\@orgs, '01.fas');
#runCmd("usearch.pl -i 01.fas -o 11");
clu_group("11.clu", "21.tbl", \@orgs);
select_clu_seq("21.tbl", "01.fas", "23.fas", 50);

sub merge_seq {
  my ($orgs, $fo) = @_;
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  my $seqHO = Bio::SeqIO->new(-fh=>$fho, -format=>'fasta');
  for my $org (@$orgs) {
    my $fi = "../$org\_HM101/41_novseq/21.fas";
    open(my $fhi, "<$fi") or die "cannot read $fi\n";
    my $seqHI = Bio::SeqIO->new(-fh=>$fhi, -format=>'fasta');
    while(my $seqo = $seqHI->next_seq()) {
      my $oid = $seqo->id;
      $oid =~ /^([\w\-]+)\-(\d+)\-(\d+)\-(\d+)\-(\d+)$/
        || die "unknown id [$oid] in $fi\n";
      my ($id, $b, $e, $rb, $re) = ($1, $2, $3, $4, $5);
      my ($beg, $end) = ($b + $rb -1, $b + $re - 1);
      my $nid = "$org-$id-$beg-$end";
      $seqo->id($nid);
      $seqHO->write_seq($seqo);
    }
    close $fhi;
  }
  close $fho;
}
sub clu_group {
  my ($fi, $fo, $org_order) = @_;
  open(my $fhi, "<$fi") or die "cannot read $fi\n";
  my $hc;
  while(<$fhi>) {
    chomp;
    my ($id, $beg, $end, $srd, $cid) = split "\t";
    $hc->{$cid} ||= [];
    push @{$hc->{$cid}}, [$id, $beg, $end, $srd];
  }
  close $fhi;

  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho join("\t", qw/cid len org orgs cnts strs/)."\n";
  for my $cid (sort {$a <=> $b} keys(%$hc)) {
    my (@strs, $ho);
    my $len;
    for (@{$hc->{$cid}}) {
      my ($id, $beg, $end, $srd) = @$_;
      push @strs, "$id:$beg:$end:$srd";

      $len ||= $end - $beg + 1;
      $len == $end - $beg + 1 || die "len conflict: cid $cid\n";
      my ($org) = split("-", $id);
      $ho->{$org} ||= 0;
      $ho->{$org} ++;
    }
    my @orgs = grep {exists $ho->{$_}} @$org_order;
    my @cnts = map {$ho->{$_}} @orgs;
    my $org = $orgs[0];
    my $str_org = join(",", @orgs);
    my $str_cnt = join(",", @cnts);
    my $str_loc = join(" ", @strs);
    print $fho join("\t", $cid, $len, $org, $str_org, $str_cnt, 
      $str_loc)."\n";
  }
  close $fho;
}
sub select_clu_seq {
  my ($fi, $fs, $fo, $minlen) = @_;
  my $t = readTable(-in => $fi, -header => 1);
  my $db = Bio::DB::Fasta->new($fs);
  my $seqHO = Bio::SeqIO->new(-file => ">$fo", -format=>'fasta');
  for my $i (0..$t->nofRow-1) {
    my ($cid, $len, $org, $str_org, $str_cnt, $str_loc) = $t->row($i);
    next if $len < $minlen;
    my @strs = split(" ", $str_loc);
    my $idx = first_index {/^\Q$org\E/} @strs;
    my $str = $strs[$idx];
    my ($id, $beg, $end, $srd) = split(":", $str);
    my $seq = $db->seq($id, $beg, $end);
    $seq = Bio::Seq->new(-seq=>$seq)->revcom->seq if $srd eq "-";
    my $nid = "$cid|$str";
    $seqHO->write_seq(Bio::Seq->new(-id=>$nid, -seq=>$seq));
  }
  $seqHO->close();
}

#make_chr("$dir/$org.tbl", "$data/genome/$org/11_genome.fa", $org, "$dir/$org");
#make_chr("$dir/$org\_unc.tbl", "$data/genome/$org/11_genome.fa", "$org\_unc", "$dir/$org\_unc");
sub make_chr {
  my ($fi, $fs, $tag, $fop) = @_;
  my $foa = "$fop.agp";
  
  open(FHO, ">$foa") or die "cannot write to $foa\n";
  print FHO join("\t", qw/id beg end srd sid sbeg send len/)."\n";
  my $t = readTable(-in=>$fi, -header=>1);
  my $seq = "";
  my $l = 0;
  for my $i (0..$t->nofRow-1) {
    my ($id, $beg, $end, $len) = $t->row($i);
    my $seq1 = seqRet([[$beg, $end]], $id, "+", $fs);
    if($i > 0) {
      $seq .= "N" x 50;
      $l += 50;
    }
    print FHO join("\t", $tag, $l+1, $l+$len, "+", $id, $beg, $end, $len)."\n";
    $seq .= $seq1;
    $l += $len;
  }
  close FHO;

  my $seqH = Bio::SeqIO->new(-file=>">$fop.fa", -format=>"fasta");
  $seqH->write_seq(Bio::Seq->new(-id=>$tag, -seq=>$seq));
  $seqH->close();
}



