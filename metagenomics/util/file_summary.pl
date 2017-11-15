#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use feature 'say';
use Getopt::Long;
use Pod::Usage;

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

    say "OK";
}

my $file = shift;

open my $fh, '<', $file;

my ($num_lines, $num_char) = (0, 0);

while (my $line = <$fh>) {
    $num_lines++;
    $num_char += length($line);
}

say "Reads = $num_lines, chars = $num_char;
# --------------------------------------------------
sub get_args {
    my %args;
    GetOptions(
        \%args,
        'help',
        'man',
    ) or pod2usage(2);

    return %args;
}

__END__

# --------------------------------------------------

=pod

=head1 NAME

file_summary - a script

=head1 SYNOPSIS

  file_summary 

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

Copyright (c) 2015 jamesthornton

This module is free software; you can redistribute it and/or
modify it under the terms of the GPL (either version 1, or at
your option, any later version) or the Artistic License 2.0.
Refer to LICENSE for the full license text and to DISCLAIMER for
additional warranty disclaimers.

=cut
