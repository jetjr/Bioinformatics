#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use feature 'say';
use Cwd 'cwd';
use Getopt::Long;
use Pod::Usage;
use Bio::SeqIO;
use Bio::Tools::SeqStats;
use File::Basename;
use File::Spec::Functions;
use File::Path 'make_path';

main();

# --------------------------------------------------
sub main {
     
     unless (@ARGV) {
        die "Please provide a FASTA file.\n";
}
     my %args = get_args();

    if ($args{'help'} || $args{'man_page'}) {
        pod2usage({
            -exitval => 0,
            -verbose => $args{'man_page'} ? 2 : 1
        });
    }; 

    

my $out_dir = $args{'out_dir'} || cwd();

unless (-d $out_dir) {
    make_path($out_dir);
}

for my $file (@ARGV) {
    my $count = 0;
    my $seq_in = Bio::SeqIO->new( -file => $file, -format => 'Fasta');
    my ($total_base, $total_gc);
    my $out_name = catfile($out_dir, basename($file) . '.' . "gc");
    open (my $fh, '>', $out_name) or die "Could not open file '$out_name'";
    while (my $seq = $seq_in->next_seq) {
        $count++;
        my $seq_stats = Bio::Tools::SeqStats->new('-seq'=>$seq);
        my $hash_ref = $seq_stats->count_monomers();   
        say $fh "---------------------------------";    
        say $fh "ID$count: ", $seq->id;
        say $fh "Sequence: ", $seq->seq();
        say $fh "Length: ", $seq->length;
        $total_base += $seq->length;
        $total_gc += $hash_ref->{'G'} + $hash_ref->{'C'};   
        printf $fh "GC: %.2f%%\n", ($hash_ref->{'G'} + $hash_ref->{'C'}) / $seq->length() * 100;
        
        }

    print "File successfully created: ($out_name)\n";
    close $fh;
    }
}
    
# --------------------------------------------------
sub get_args {
    my %args;
    GetOptions(
        \%args,
        'help',
        'man',
        'out_dir=s'
    ) or pod2usage(2);

    return %args;
}

__END__

# --------------------------------------------------

=pod

=head1 NAME

gc_content.pl - a script

=head1 SYNOPSIS

  ./gc_content.pl [File(s)] -options

Options:

  --out    Names the out directory for file (Default current directory)
  --help   Show brief help and exit
  --man    Show full documentation

=head1 DESCRIPTION

Calculates the GC content of EACH sequence found in a FASTA file. 

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
