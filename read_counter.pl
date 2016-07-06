#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use feature 'say';
use Getopt::Long;
use Pod::Usage;
use Bio::SeqIO;
use File::Basename;

my @files_list;

main();

# --------------------------------------------------
sub main {

    my %args = get_args();

    if ($args{'help'} || $args{'man_page'}) {
        pod2usage({
            -exitval => 0,
            -verbose => $args{'man_page'} ? 2 : 1
        });
    }; 

my $out_file= $args{'out_file'} || "read_counts.txt";
    
open my $fh, '>', $out_file;
    print $fh join("\t", "Sample Name", "Read Count"), "\n";
    printf ("%-30s %10s\n", "Sample Name", "Read Count");
    foreach my $file (@files_list) {
        my $in = Bio::SeqIO -> new(-file => $file, -format => "fasta");
        my $count = 0;
        while ($in -> next_seq ) {
            $count++;
        }
        printf ("%-30s %-10s\n", basename($file), $count);
        print $fh join("\t", basename($file), $count), "\n"; 
    }
    say "Read counts sucessfully stored in $out_file";
}

# --------------------------------------------------
sub get_args {
    my %args;
    GetOptions(
        \%args,
        'files|f=s@{1,}' => \@files_list,
        'out_file=s',
        'help',
        'man',
    ) or pod2usage(2);

    return %args;
}

__END__

# --------------------------------------------------

=pod

=head1 NAME

read_counter.pl - a script to count reads in a FASTA file

=head1 SYNOPSIS

  read_counter.pl -f [FASTA file or files] 

Options:

  --help   Show brief help and exit
  --man    Show full documentation

=head1 DESCRIPTION

Counts the number of reads in a FASTA file. Input can be a single or multiple FASTA files.

=head1 SEE ALSO

perl.

=head1 AUTHOR

James Eric Thornton E<lt>jamesthornton@email.arizona.eduE<gt>.

=head1 COPYRIGHT

Copyright (c) 2016 jamesthornton

This module is free software; you can redistribute it and/or
modify it under the terms of the GPL (either version 1, or at
your option, any later version) or the Artistic License 2.0.
Refer to LICENSE for the full license text and to DISCLAIMER for
additional warranty disclaimers.

=cut
