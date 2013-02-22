package Eutils;
use strict; 
use Bio::Seq;
use Bio::DB::EUtilities;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw/get_gi_note/;
@EXPORT_OK = qw//;

sub get_gi_note {
#my @ids = qw/223556040 62899127 119359623/;
#my @notes = get_gi_note(@ids);
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

1;
__END__

