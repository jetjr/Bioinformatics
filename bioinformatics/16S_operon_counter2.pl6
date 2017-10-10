#!/usr/bin/env perl6

sub MAIN (Str $file!) {
    die "Not a file ($file)" unless $file.IO.f;
    
    my $seq = $file.IO.lines.join;
    my @count = $seq ~~ m:g/AGAGTTTGATC<[AC]>TGGCTCAG/;   
    
    put join " ", "I found:", @count.elems, "operons";
    put @count.chars;
    put $seq.chars;

    my $out-file = open $file ~ '.operons', :w;
    $out-file.print("I found: ", @count.elems, " operons in $file\n");
}

sub USAGE {
    printf
      "Usage:\n %s --fasta=<File>\n",
      $*PROGRAM-NAME;
}
