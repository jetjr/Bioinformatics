#!/usr/bin/env perl

use strict;
use warnings;

my $string = shift;
$string =~ tr/ATCG/TAGC/g ;
print $string;

