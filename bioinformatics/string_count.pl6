#!/usr/bin/env perl6

sub MAIN (Str $input!) {
  my $bag = $input.comb.Bag;
  for $bag.keys.grep(/<[a..zA..Z]>/).sort -> $char {
    put join "\t", $char, $bag{$char};
  }
}
