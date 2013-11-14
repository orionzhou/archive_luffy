package Medicago;
use strict;
use Common;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw/get_mt_ids/;
@EXPORT_OK = qw//;
sub get_mt_ids {
    my ($opt) = @_;
    my $fi = $ENV{"code"}."/conf/acc_ids.tbl";
    my $t = readTable(-in=>$fi, -header=>1);

    my $idx = -1;
    if($opt =~ /^(acc|opt)([\d\.]+)$/) {
        $idx = first_index {$_ eq $opt} $t->header;
    } 
    
    $idx == -1 || die "unknown opt $opt\n";
    $t = $t->match_pattern("\$_->[$idx] eq 1");
    my $ids = $t->colRef("id");
    return $ids;
}


1;
__END__
