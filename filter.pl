#!/usr/bin/perl -w
use strict;

my $input = shift;
my $seqtot = shift;
my $seqvar = shift;
my $rnatot = shift;
my $rnavar = shift;
my $valtot = shift; # validation - solid coverage
my $valvar = shift;

open IN, "$input" or die "Can't read";
while(<IN>) {
  my $line = $_;
  if ($line =~ /^#|^contig/) {
    print $line;
    next;
  }

  my @ln = split(/\t/, $line);
  
  if ($ln[39] >= $seqtot && $ln[40] >= $seqvar) {     # sequencing = illumina coverage
    if ($ln[35] >= $valtot && $ln[36] >= $valvar) {   # validation = solid coverage
      print $line;
    } elsif ($ln[37] >= $rnatot && $ln[38] >= $rnavar) { # rna coverage
      print $line;
    }
  } 
}
close IN; 

