#!/usr/bin/env perl6

sub MAIN (Str $dna) {
  my $new = $dna.trans(<A C G T> => < T G C A>);
  print $new;
}
