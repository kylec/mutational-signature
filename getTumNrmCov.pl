#!/usr/bin/perl -w
use Analysis;
use strict;
use Getopt::Long;

=usage
	take a maf file and add normal and tumor coverage to it
=cut

my $fin;
my $fout;
my $bamdir; 	# ~jenn/datalinks/coca/illumina/we/bam/full
my $covLabel; # wgs, solid
my $shortSample;

GetOptions(
	'i|input=s' => \$fin,
	'bamdir=s'=>\$bamdir,
	'label=s'=>\$covLabel,
	'short'=>\$shortSample
);

$fout = "$fin.cov";
my $fa = "/home/kchang3/references/hsap_36.1_hg18.fa";

# input/output
open OUT, ">$fout" or die "Can't write output";
open IN, $fin or die "Can't read input";
#skip comment
<IN>;
my $header = <IN>;
chomp $header;
print "$header\tTTot\tTVar\tNTot\tNVar\n";
while(<IN>) {
  chomp $_;
  my $ln = $_;
  my @ln = split(/\t/, $ln);
  my $chr = $ln[0];
  my $start = $ln[1];
  my $ref = $ln[3];
  my $var = $ln[4];
  my $tsam = $ln[5];
  my $nsam = $ln[6];

  if ($shortSample) {
    $tsam = Analysis->processSampleName($tsam, "a");
    $nsam = Analysis->processSampleName($nsam, "a");
  }

  ### Init cov
  my $ttot = my $tvar = my $ntot = my $nvar = 0;
	
  ###  BAMS
  my $tbam = `ls $bamdir/$tsam*.bam | head -1`;
  my $nbam = `ls $bamdir/$nsam*.bam | head -1`;
  chomp $tbam; chomp $nbam;
	 
  # Extract coverage
  ($ttot, $tvar) = Analysis->getPileupCovByBamAndCoord($chr, $start, $ref, $var, $tbam, $fa) if ($tbam);
  ($ntot, $nvar) = Analysis->getPileupCovByBamAndCoord($chr, $start, $ref, $var, $nbam, $fa) if ($nbam);
	
  print OUT "$ln\t$ttot\t$tvar\t$ntot\t$nvar\n";
}
close OUT;
close IN;
