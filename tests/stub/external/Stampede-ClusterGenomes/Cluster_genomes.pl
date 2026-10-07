#!/usr/bin/env perl
# Stub Stampede Cluster_genomes.pl -f <fasta> -i ID -c COV: every sequence is its own cluster.
use strict; use warnings;
my %a = @ARGV; my $f = $a{'-f'} or die "need -f\n";
open(my $fh, '<', $f) or die "$f: $!\n";
my $n = 0;
while (<$fh>) { if (/^>(\S+)/) { $n++; print "Cluster $n\n$1\n"; } }
