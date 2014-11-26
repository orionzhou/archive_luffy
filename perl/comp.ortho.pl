#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  comp.ortho.pl - 

=head1 SYNOPSIS
  
  comp.ortho.pl [-help]

  Options:
    -h (--help)   brief help message

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use File::Path qw/make_path remove_tree/;
use File::Basename;
use Tabix;
use Bio::DB::Fasta;
use Data::Dumper;
use Common;
use Location;
use Gtb;
use Gal;
use List::Util qw/min max sum/;

my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
) or pod2usage(2);
pod2usage(1) if $help_flag;

my $dir = "$ENV{'misc3'}/comp.ortho";
-d $dir || make_path($dir);
chdir $dir || die "cannot chdir to $dir\n";

my $tgt = "HM101";
my @qrys = qw/
  HM058 HM125 HM056 HM129 HM060
  HM095 HM185 HM034 HM004 HM050 
  HM023 HM010 HM022 HM324 HM340
/;
my @orgs = ($tgt, @qrys);

#extract_unortho_seq("06.no.ortho.tbl", "07.seq");
#runCmd("qsub.blat.pl -i 11.fas -o 12.blat");
#runCmd("seq 8 15 | xargs -i printf \"%02d\\n\" {} | parallel -j 8 blat -prot -out=blast8 11.fas 12.blat/part.{}.fas 12.blat/part.{}.tbl");
#runCmd("cat 12.blat/part.*.tbl > 12.blat.tbl");
#runCmd("sort -k1,1 -T /scratch/zhoup 14.pairs.tbl -o 15.pairs.sort.tbl");
#get_best_hit('15.pairs.sort.tbl', '21.best.pairs.tbl');
#runCmd("\$soft/mcl/bin/mcl 21.best.pairs.tbl -te 16 -I 5.0 --abc -o 25.txt");
#parse_mcl("25.txt", "26.tbl", \@orgs);
#ortho_group_refine("01.ortho.tbl", "26.tbl", "31.ortho.tbl", \@orgs);
ortho_group_cat("31.ortho.tbl", "33.ortho.cat.tbl", \@orgs);

sub extract_unortho_seq {
  my ($fi, $do) = @_;
  -d $do || make_path($do);
  
  my $hi;
  my $ti = readTable(-in => $fi, -header => 1);
  for my $i (0..$ti->lastRow) {
    my ($org, $id) = $ti->row($i);
    $hi->{$org} ||= [];
    push @{$hi->{$org}}, $id;
  }

  my @orgs = sort(keys(%$hi));

  for my $org (@orgs) {
    my $fs = "$ENV{'genome'}/$org/51.fas";
    my $db = Bio::DB::Fasta->new($fs);
    my $fo = "$do/$org.fasta";
    my $seqHO = Bio::SeqIO->new(-file => ">$fo", -format => 'fasta');
    for my $id (@{$hi->{$org}}) { 
      my $seq = $db->seq($id);
      $seqHO->write_seq(Bio::Seq->new(-id => "$org|$id", -seq => $seq));
    }
    $seqHO->close();
  }
}
sub write_one_hit {
  my ($tid, $h, $fho) = @_;
  for my $org (keys(%$h)) {
    my ($qid, $score) = @{$h->{$org}};
    print $fho join("\t", $tid, $qid, $score)."\n";
  }
}
sub get_best_hit {
  my ($fi, $fo) = @_;
  open(my $fhi, "<$fi") or die "cannot read $fi\n";
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  my $h;
  my $ptid = "";
  while(<$fhi>) {
    chomp;
    my ($tid, $qid, $torg, $qorg, $e1, $e2, $ident, $cov) = split "\t";
    next if $torg eq $qorg || $cov < 0.5;
    my $score = $ident * $cov / 10000;
    if($ptid ne $tid && $ptid ne "") {
      write_one_hit($ptid, $h, $fho);
      $h = {$qorg => [$qid, $score]};
    } else {
      $h->{$qorg} ||= [$qid, $score];
      $h->{$qorg} = [$qid, $score] if $score > $h->{$qorg}->[1];
    }
    $ptid = $tid;
  }
  write_one_hit($ptid, $h, $fho);
  close $fhi;
  close $fho;
}
sub parse_mcl {
  my ($fi, $fo, $orgs) = @_;
  open(my $fhi, "<$fi") or die "cannot read $fi\n";
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho join("\t", @$orgs)."\n";
  while(<$fhi>) {
    chomp;
    my @ps = split "\t";
    my $h = { map {$_ => []} @$orgs };
    my $hn = { map {$_ => 0} @$orgs };
    for (@ps) {
      my ($org, $id) = split /\|/;
      push @{$h->{$org}}, $id;
      $hn->{$org} ++;
    }
    my $n = max(values(%$hn));
    for my $i (1..$n) {
      my @ort;
      for my $org (@$orgs) {
        if($hn->{$org} >= $i) {
          push @ort, $h->{$org}->[$i-1];
        } else {
          push @ort, '';
        }
      }
      print $fho join("\t", @ort)."\n";
    }
  }
  close $fhi;
  close $fho;
}
sub ortho_group_refine {
  my ($fi, $fr, $fo, $orgs) = @_;
  my %hi = map {$_ => $orgs->[$_]} (0..@$orgs-1);
  my $ti = readTable(-in => $fi, -header => 1);
  $ti->delCols(['cat', 'n_org']);

  my $tr = readTable(-in => $fr, -header => 1);
  my $hr;
  for my $i (0..$tr->lastRow) {
    my @ps = $tr->row($i);
    if($ps[0] ne '') {
      exists $hr->{$ps[0]} && die "$ps[0] appeared >1 times in $fr\n";
      $hr->{$ps[0]} = \@ps;
    }
  }

  my $nc = 0;
  for my $i (0..$ti->lastRow) {
    my @ps = $ti->row($i);
    my $id = $ps[0];
    exists $hr->{$id} || next;
    my @psr = @{$hr->{$id}};
    my @psn = ('') x @$orgs;
    for my $j (1..@$orgs-1) {
      if($ps[$j] eq 'NA') {
        $ps[$j] = $psr[$j] if $psr[$j] ne '';
      } else {
        if($psr[$j] ne '' && $psr[$j] ne $ps[$j]) {
          $psn[$j] = $psr[$j];
        }
      }
    }
    my $nd = scalar( grep {$_ ne ''} @psn );
    if($nd > 0) {
      $nc ++;
      $ti->addRow(\@psn);
    }
  }
  print "$nc\n";
  
  for my $i (0..$tr->lastRow) {
    my @ps = $tr->row($i);
    $ti->addRow(\@ps) if $ps[0] eq '';
  }

  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho $ti->tsv(1);
  close $fho;
}
sub ortho_group_cat {
  my ($fi, $fo, $orgs) = @_;
  my $ti = readTable(-in => $fi, -header => 1);
  my $h;
  for my $org (@$orgs) {
    my $fg = "$ENV{'genome'}/$org/51.gtb";
    my $tg = readTable(-in => $fg, -header => 1);
    my %hg = map {$tg->elm($_, 'id') => 
      [$tg->elm($_, 'cat2'), $tg->elm($_, 'cat3')]} (0..$tg->lastRow);
    $h->{$org} = \%hg;
  }

  my (@cats2, @cats3);
  for my $i (0..$ti->lastRow) {
    my @ps = $ti->row($i);
    for my $j (0..@$orgs-1) {
      if($ps[$j] ne '' && $ps[$j] ne 'NA') {
        my $org = $orgs->[$j];
        push @cats2, $h->{$org}->{$ps[$j]}->[0];
        push @cats3, $h->{$org}->{$ps[$j]}->[1];
        last;
      }
    }
  }
  $ti->addCol(\@cats2, 'cat2');
  $ti->addCol(\@cats3, 'cat3');
  open(my $fho, ">$fo") or die "cannot write $fo\n";
  print $fho $ti->tsv(1);
  close $fho;
}


__END__

