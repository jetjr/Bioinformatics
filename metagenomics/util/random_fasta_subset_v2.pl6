#!/usr/bin/env perl6

sub MAIN {
    my $count = 100;
    my $percent = 0.20;
    my @take = (0...$count).BagHash.grab(round($count * $percent)).sort;
    print @take;
}
