#!/usr/bin/env perl6

subset File of Str where *.IO.f;

sub MAIN (File :$db, File :$ids, Str :$out='out.txt') {
    my %swissprot;
    put "Reading $db";
    for $db.IO.lines -> $line {
        my ($, $sp_id, $tax_id, $) = $line.split(/\t/);
        %swissprot{ $sp_id } = $tax_id;
    }

    my $fh = open $out, :w;

    put "Reading $ids";
    for $ids.IO.lines -> $id {
        if (my $tax_id = %swissprot{ $id }) {
            $fh.print("$id = $tax_id\n");
        }
    }

    put "Done.";
}
