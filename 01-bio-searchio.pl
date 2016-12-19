#!/usr/bin/env perl

use warnings;
use strict;
use feature 'say';
use File::Basename;
use Bio::SearchIO;

main();
#---------------------------------------

sub main {

my $blast_out = shift or die sprintf("Usage: %s blast.out\n", basename($0));

my %args = get_args();

    if ($args{'help'} || $args{'man_page'}) {
        pod2usage ({
            -exitval => 0,
            -verbose => $args{'man_page'} ? 2 : 1
    });
};

my $min       = 1e-10 || $args{'e_value'};
my $searchio  = Bio::SearchIO->new( 
    -format   => 'blast',
    -file     => $blast_out,
);

say join "\t", qw[query hit evalue];
while (my $result = $searchio->next_result ) {
    while (my $hit = $result->next_hit) {
        while (my $hsp = $hit->next_hsp) {
            if ($hsp->significance <= $min) {
                say join "\t", 
                    $result->query_name,
		    $hit->name,
	        $hit->description,
		    $hsp->percent_identity,
		    $hsp->length('total'),
 		    $hsp->evalue;
            }
        }
    }
}
}
#--------------------------------------
sub get_args {
    my %args;
    GetOptions (
        \%args,
        'help',
        'man',
        'e_value=i'
    ) or pod2usage(2);
    
    return %args;
}
__END__

#----------------------------------------

=pod 

=head1 NAME

biosearch.pl - a script

=head1 SYNOPSIS

  ./biosearch.pl [blast.out] -options

Options:

    --e_value   Provides the maximum e_value (default 1e-10)
    --help      Show brief help and exit
    --man       Show full documentation

=head1 DESCRIPTION

Under construction

=head1 SEE ALSO

perl.

=head1 AUTHOR

James Eric Thornton E<lt>jamesthornton@email.arizona.eduE<gt>.

=head1 COPYRIGHT

Copyright (c) 2016 jamesthornton

=cut
