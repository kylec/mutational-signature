#!/usr/bin/perl -w

use Analysis;

my $fin = shift;
my $fasta = "/home/kchang3/references/hsap_36.1_hg18.fa";

open IN, "$fin" or die "Can't read $fin";
my $header = <IN>;
chomp $header;
print "$header\tSeqContext\n";
while(<IN>) {
  chomp;
  my $line = $_;
  my @ln = split(/\t/, $_);  
  my $chr = "chr$ln[4]";
  my $start = $ln[5]-3;
  my $end = $ln[6]+3;
  my $ref = $ln[10];
  my $refbase = Analysis->getRefBases($chr, $start, $end, $fasta);
  my @bases = split(//, $refbase);
  my $seqcontext = ".";
  #print "ref = $refbase\n";
  #print "array = @bases\n";
  #print "$bases[3],$ref\n";
  if ($bases[3] eq $ref) {
    $seqcontext = "$bases[0]$bases[1]$bases[2]x$bases[4]$bases[5]$bases[6]";
  } else {
    die "error\t$line\n";
  }
  print "$line\t$seqcontext\n";
}
close IN;
