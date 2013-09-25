package Eutils;
use strict; 
use Bio::Seq;
use Bio::DB::EUtilities;
use Bio::DB::Taxonomy;
use Data::Dumper;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw/gi2Taxid annotate_taxid/;
@EXPORT_OK = qw//;

our $dir = "/home/youngn/zhoup/Data/db/ncbi_taxon";
our $fnodes = "$dir/nodes.dmp";
our $fnames = "$dir/names.dmp";

sub get_gi_note {
    my @ids = @_;
    my $fac = Bio::DB::EUtilities->new(-eutil=>'esummary', -email=>'zhoupenggeni@gmail.com', -db=>'nuccore', -id=>\@ids);
    my @notes;
    while(my $ds = $fac->next_DocSum) {
        my @items = $ds->get_Items_by_name("Title");
        my $note = $items[0]->get_content;
        push @notes, $note;
    }
    return @notes;
}
sub gi2Taxid {
    my @gis = @_;
    my $h = { map {$_=>''} @gis };
    @gis = keys(%$h);
    my $fac = my $factory = Bio::DB::EUtilities->new(-eutil=>'esummary', 
        -email=>'zhoupenggeni@gmail.com',
        -db=>'nucleotide', -id=>\@gis);
    while( my $ds = $fac->next_DocSum) {
        my $gi;
        while(my $it = $ds->next_Item("flattened")) {
            $gi = $it->get_content if $it->get_name eq "Gi";
            if($it->get_name eq "TaxId") {
                die "no GI error\n" unless $gi;
                $h->{$gi} = $it->get_content;
            }
        }
    }
    return $h;
}
sub gi2Taxon_entrez {
    my @gis = @_;
    my $db = Bio::DB::Taxonomy->new(-source=>'entrez');
    my $h1 = { map {$_=>''} @gis };
    my $h2 = { map {$_=>''} @gis };
    my @ids = keys(%$h1);
    for my $i (0..$#ids) {
        my $gi = $ids[$i];
        my $node = $db->get_Taxonomy_Node(-gi=>$gi, -db=>'nucleotide');
        $h2->{$gi} = $node->scientific_name;
        my %l;
        my @ary;
        for my $j (1..30) {
            my $pa = $node->ancestor();
            last if(!$pa);
            if($pa->rank eq "no rank") {
                push @ary, $pa->scientific_name;
            } else {
                $l{$pa->rank} = $pa->scientific_name;
            }
            $node = $pa;
        }
        my $cat1;
        if(!exists($l{'kingdom'}) && !exists($l{'superkingdom'})) {
            $cat1 = $ary[-1];
        } elsif(exists($l{'kingdom'})) {
            $cat1 = $l{'kingdom'};
            if($l{'kingdom'} =~ /(Viridiplantae)|(Metazoa)/i) {
                $cat1 = $l{"family"} if exists $l{"family"};
            }
        } elsif(exists($l{'superkingdom'})) {
            $cat1 = $l{'superkingdom'};
        }
        $h1->{$gi} = $cat1;
        printf "%4d | %4d\n", $i+1, scalar(@ids) if ($i+1) % 10 == 0;
    }
#    my @cats1 = map {$h1->{$_}} @gis;
#    my @cats2 = map {$h2->{$_}} @gis;
#    return (\@cats1, \@cats2);
    return ($h1, $h2);
}
sub annotate_taxid {
    my @ids = @_;
    my $h = { map {$_=>''} @ids };
    @ids = keys(%$h);

    my $db = Bio::DB::Taxonomy->new(-source=>'flatfile', -nodesfile=>$fnodes, -namesfile=>$fnames, -directory=>$dir);
    for my $i (0..$#ids) {
        my $id = $ids[$i];
        my $node = $db->get_taxon(-taxonid=>$id);
        my %l;
        my @ary;
        for my $j (1..30) {
            my $pa = $node->ancestor();
            last if(!$pa);
            if($pa->rank eq "no rank") {
                push @ary, $pa->scientific_name;
            } else {
                $l{$pa->rank} = $pa->scientific_name;
            }
            $node = $pa;
        }
        my $cat1;
        if(!exists($l{'kingdom'}) && !exists($l{'superkingdom'})) {
            $cat1 = $ary[-1];
        } elsif(exists($l{'kingdom'})) {
            $cat1 = $l{'kingdom'};
            if($l{'kingdom'} =~ /(Viridiplantae)|(Metazoa)/i) {
                $cat1 = $l{"family"} if exists $l{"family"};
            }
        } elsif(exists($l{'superkingdom'})) {
            $cat1 = $l{'superkingdom'};
        }
        $h->{$id} = [$node->scientific_name, $cat1];
        printf "%4d | %4d\r", $i+1, scalar(@ids) if ($i+1) % 10 == 0;
    }
    return $h;
}




1;
__END__

