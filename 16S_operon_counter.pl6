#!/usr/bin/env perl6

sub MAIN (Str :$file!, Str :$fprimer! is copy, Str :$rprimer! is copy) {
    die "Not a file ($file)" unless $file.IO.f;
    
    my %xform = ( 
        R => '<[AG]>',
        Y => '<[CT]>',
        S => '<[GC]>',
        W => '<[AT]>',
        K => '<[GT]>',
        M => '<[AC]>',
        B => '<[CGT]>',
        D => '<[AGT]>',
        H => '<[ACT]>',
        V => '<[ACG]>',
        N => '<[ACTG]>'
    );

    for %xform.kv -> $base, $new {
        $fprimer.subst-mutate($base, $new);
        $rprimer.subst-mutate($base, $new);
    }
    
    my $seq = $file.IO.lines.join;

    my @count = $seq ~~ m:g/$fprimer.**? 5..*$rprimer/;   
    
    put join " ", "I found:", @count.elems, "operons";
    put @count.chars;
    put $seq.chars;

    my $out-file = open $file ~ '.operons', :w;
    $out-file.print("I found: ", @count.elems, " operons in $file\n");
    $out-file.print("Forward primer: ", $fprimer, "Reverse primer: ", $rprimer, "\n");
    $out-file.print("Operon sequences: ", @count,"\n");
    $out-file.print("Operon total size: ", @count.chars, "\n");
    $out-file.print("Total sequence size: ", $seq.chars, "\n");
}

sub USAGE {
    printf
      "Usage:\n %s --fasta=<File> --fprimer=<forward primer sequence> --rprimer=<reverse primer sequence>\n",
      $*PROGRAM-NAME;
}
