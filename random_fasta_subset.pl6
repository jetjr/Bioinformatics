#!/usr/bin/env perl6

subset File of Str where *.IO.f;

my $num_lines=0;
my $i=1;
my $num_sequences=0;
my @odds;
my @l;
my $num_reads=0;
my $header;
my $sequence;

sub MAIN (File :$fasta, Str :$percent) {
    for $fasta.IO.lines -> $line {
        $num_lines++;
        if $line ~~ /^ \> / {
            $num_sequences++;
        }
    }
    for (1...$num_lines) -> $line_num { 
        if $line_num %% 2 { next; }
        else {@odds.push($line_num);}
}
    $num_reads=($num_sequences * $percent);
    $num_reads=$num_reads.round;
    
    repeat while $i <= $num_reads {
        $header = @odds.pick;
        $sequence = ($header + 1);
        if $header eq any(@l) { next; }
        else { @l.push($header, $sequence), $i++; }
        
    }   

    my $count = 1;
    for $fasta.IO.lines -> $seq {
        if $count eq any(@l) { print "$seq\n"; }
        $count++;
    }
}

