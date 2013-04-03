package Eutils;
use strict; 
use Bio::Seq;
use Bio::DB::EUtilities;
use Bio::DB::Taxonomy;
use Data::Dumper;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw/get_gi_taxonomy/;
@EXPORT_OK = qw//;

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
sub get_gi_taxonomy {
    my @gis = @_;
    my $db = Bio::DB::Taxonomy->new(-source=>'entrez');
    my $h1 = { map {$_=>''} @gis };
    my $h2 = { map {$_=>''} @gis };
    for my $gi (keys(%$h1)) {
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
    }
    my @cats1 = map {$h1->{$_}} @gis;
    my @cats2 = map {$h2->{$_}} @gis;
    return (\@cats1, \@cats2);
}

1;
__END__

