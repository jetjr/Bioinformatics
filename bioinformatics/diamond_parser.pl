#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use feature 'say';
use Getopt::Long;
use Pod::Usage;


if (!@ARGV) {
    pod2usage(-verbose => 0, -message => "$0: FASTA file required\n")
}

my @matches;

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

    my $database= $args{'database'};
    my $out_file= $args{'out'};

    open my $fh, '<', @ARGV;
    open my $db, '<', $database;
    open my $out, '>', $out_file;

    my %swissprot;
    while (my $line = <$db>) {
        chomp($line);
        my ($md5, $sp_id, $tax_id, $src) = split /\t/, $line;
        $swissprot{ $sp_id } = $tax_id;
    }

    while (my $id = <$fh>) {
        chomp($id);
        if (my $tax_id = $swissprot{ $id }) {
            print $out "$id = $tax_id\n";
        }
    }
}

# --------------------------------------------------
sub get_args {
    my %args;
    GetOptions(
        \%args,
        'help',
        'man',
        'database=s',
        'out=s',
    ) or pod2usage(2);

    return %args;
}

__END__

# --------------------------------------------------

=pod

=head1 NAME

diamond_parser.pl - a script

=head1 SYNOPSIS

  diamond_parser.pl 

Options:

  --help   Show brief help and exit
  --man    Show full documentation

=head1 DESCRIPTION

Describe what the script does, what input it expects, what output it
creates, etc.

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
