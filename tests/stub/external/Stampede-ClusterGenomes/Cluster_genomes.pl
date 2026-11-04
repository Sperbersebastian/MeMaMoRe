#!/usr/bin/env perl
# Stub Stampede Cluster_genomes.pl -f <fasta> -i ID -c COV: every sequence is its own
# cluster. Writes like the real tool: verbose log on stdout, <stem>_<ID>-<COV>.fna
# (representatives) and <stem>_<ID>-<COV>.clstr (">Cluster_N<TAB>rep<TAB>len").
use strict; use warnings;
my %a = @ARGV; my $f = $a{'-f'} or die "need -f\n";
my ($id, $cov) = ($a{'-i'} // 95, $a{'-c'} // 80);
(my $stem = $f) =~ s/\.(fasta|fa|fna)$//;
open(my $fh, '<', $f) or die "$f: $!\n";
my (@ids, %seq, $cur);
while (<$fh>) { chomp; if (/^>(\S+)/) { $cur = $1; push @ids, $cur; $seq{$cur} = ''; } elsif (defined $cur) { $seq{$cur} .= $_; } }
open(my $fo, '>', "${stem}_${id}-${cov}.fna") or die;
open(my $fc, '>', "${stem}_${id}-${cov}.clstr") or die;
print "comparing the contigs with mummer -> stub\n";
for my $i (0 .. $#ids) {
  print $fo ">$ids[$i]\n$seq{$ids[$i]}\n";
  printf $fc ">Cluster_%d\t%s\t%d\n", $i, $ids[$i], length $seq{$ids[$i]};
}
