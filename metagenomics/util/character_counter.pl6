#!/usr/bin/env perl6

sub MAIN (Str $file) {
  my %count;
    for $file.IO.lines -> $line {
    for $line.comb -> $letter {
      %count{ $letter }++;
  }
}
  for %count.kv -> $letter, $value {
    put join "\t", $letter, $value;
  }
}
