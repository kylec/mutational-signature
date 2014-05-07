package Analysis;

=usage
1 create snp id
2. parseFile
3. connectOracle
4. getAmpliconBySnpCoord
5. processSampleName
6. getPatientName
7. getPileupGenotype 

=cut

use strict;
#require DBI;
#require Tie::Handle::CSV;
#require Text::CSV_XS;
use List::Util qw(sum);
my $CODEDIR="/stornext/snfs6/cancer-analysis/code";

sub new {
	my $class = shift;
}

############## create snp id #################
sub getSnpidByCoordsAndType {
        my ($self, $projInit, $ampId, $global_start, $amp_start, $class ) = @_;
	
	my $run = 1;
        my $method = "SD3";
   
        my $mutationDesc;
        
        my $pos = $global_start - $amp_start + 1;

        
        if ($class =~ /SNP/i) {
                $mutationDesc = $projInit.'_'.$ampId.'_'.$pos.'_'.$method.'_'.$run;
        } elsif ($class =~ /INS/i) {
                $mutationDesc = $projInit.'_'.$ampId.'_'.$pos.'i_'.$method.'_'.$run;
        } elsif ($class =~ /DEL/i) {
                $mutationDesc = $projInit.'_'.$ampId.'_'.$pos.'g_'.$method.'_'.$run;
        } else {
                $mutationDesc = $projInit.'_'.$ampId.'_'.$pos.'_'.$method.'_'.$run;
                die "error: wrong class = $class, $global_start, $mutationDesc\n";
                
        }
        return $mutationDesc;

}

############## return structured file hash ##############
sub parseFile {
	my ($self, $input) = @_;
	# create input and output handle
	my $fileParser = Text::CSV_XS->new( {binary=> 1, sep_char => "\t", eol => $/, allow_loose_quotes => 1, quote_char=>'', escape_char=>'',auto_diag=>1} );
	my $inputHandle = Tie::Handle::CSV->new( $input, header => 1, csv_parser => $fileParser );

	return $inputHandle;
}

############## connection to Oracle #####################
sub connectOracle {
	my ($self, $id ) = @_;
	my $host; my $port; my $name; my $pwd;	
	my $dbh;
	if ($id eq "gscdevel") {
		$host = "hondo.hgsc.bcm.tmc.edu";
		$port = 1521;
		$name = "mutations";
		$pwd = "mutations";
		$dbh = DBI->connect("dbi:Oracle:sid=$id;host=$host; port=$port", $name,$pwd,
                {RaiseError => 1,
                 PrintError => 0
                }) or die "Can't connect to Oracle database: $DBI::errstr\n";


	} elsif( $id eq "gsc") {
		$host = "gsc-scan.hgsc.bcm.edu";
		$port = 1521;
		$name = "mutations";
		$pwd = "mutations";
		# gsc connection string is contructed differntly because it's on 11g instead of 10
		$dbh = DBI->connect("dbi:Oracle:$id", $name,$pwd,
                {RaiseError => 1,
                 PrintError => 0
                }) or die "Can't connect to Oracle database: $DBI::errstr\n";

	} elsif ($id eq "gscres") {
    	$host = "hgsc-ora2.hgsc.bcm.tmc.edu";
        $name = "analysis";
        $pwd = "analysis";
        $port= 1527;
    	$dbh = DBI->connect("dbi:Oracle:sid=$id;host=$host; port=$port", $name,$pwd,
                {RaiseError => 1,
                 PrintError => 0
                }) or die "Can't connect to Oracle database: $DBI::errstr\n";	
	} else {
    die "Invalid oracle SID=$id.";
  }

	$dbh->{LongReadLen}=4000;
	return $dbh;
}


######### query amplicon designs for a given site ##########
sub getAmpliconBySnpCoord {
	my ($self, $dbh, $chr, $begin, $end, $project)  = @_;
	 
	my $getAmp = "select distinct a.dna_seq_id, dp.chromosome, dp.begin, dp.end
	from star.amplicons a, star.dna_positions dp, star.dna_sequences s , star.amplicon_statuses st,star.projects p, star.project_amplicons pa
	where a.dna_seq_id=dp.dna_seq_id 
	and a.status_id = st.id
	and dp.build_id=21 
	and s.id = a.dna_seq_id
	and dp.chromosome = ?
	and dp.begin < ?
	and dp.end > ?
	and st.id in (3,4)
	and a.id = pa.amplicon_id
	and p.id = pa.project_id
	and p.name = ?";
	
	my $sth = $dbh->prepare($getAmp);
	$sth->execute($chr, $begin, $end, $project);
	
	my @results = ();
	while (my ($ampId, $ampChr, $ampStart, $ampEnd) = $sth->fetchrow_array) {
		my $ampHash = ();
		$ampHash->{id} = $ampId;
		$ampHash->{chr} = $ampChr;
		$ampHash->{begin} = $ampStart;
		$ampHash->{end} = $ampEnd;
		push(@results, $ampHash);
	}
		
	return \@results;
}

################# truncate sample name
sub processSampleName {
	my ($self, $sample, $px) = @_;
	
	my $sampleShort = "";
	if ( $px eq "triplets") {
		if ($sample =~ /(^\w+\-\w+\-\w+\-[0-9]{1})\w+\-/) {
			$sampleShort = $1;
		}elsif($sample =~ /(EPNCL\-\w+)\-/) {
			$sampleShort = $1;
		}
	} elsif ($sample =~ /(^TCGA\-\w+\-\w+\-[0-9]{2})\w+\-/) { 
		$sampleShort = $1;
	} elsif ($sample =~ /(^TARGET\-\w+\-\w+\-[0-9]{2})\w+/) {
		$sampleShort = $1;
	} elsif ( ($px eq "HC") or ($px eq "HCC") ) {
		$sampleShort = $1 if ($sample =~ /(\d+\-\d+\-[TNBR])/);
	} else {
		$sampleShort = $sample;
	}
	
	return $sampleShort;
}

sub getPatientName {
	my ($self, $sample ) = @_;
	my $patient = "";
	 
    my $numfields = `. $CODEDIR/svn/scripts/analysis/init-parms.sh -v | grep NUMFIELDS | sed 's/^NUMFIELDS *//g'`;
    chomp $numfields;
    $patient = `echo $sample | cut -d- -f1-$numfields`;
    chomp $patient;
	return $patient;
}

sub getPileupGenotype {
		# ref, var, varstring
	my ($self, $a1, $a2, $var_str)= @_;
	my %hash = ();
$hash{"M"}{"A"} = "C";
$hash{"S"}{"G"} = "C";
$hash{"W"}{"A"} = "T";
$hash{"B"}{"G"} = "T";
$hash{"R"}{"G"} = "A";
$hash{"Y"}{"T"} = "C";
$hash{"K"}{"G"} = "T";
$hash{"M"}{"C"} = "A";
$hash{"S"}{"C"} = "G";
$hash{"W"}{"T"} = "A";
$hash{"B"}{"T"} = "G";
$hash{"R"}{"A"} = "G";
$hash{"Y"}{"C"} = "T";
$hash{"N"}{"T"} = "N";
$hash{"N"}{"C"} = "N";
$hash{"N"}{"G"} = "N";
$hash{"N"}{"A"} = "N";

my %iupac = (
	"M" => "A\tC",
	"R" => "A\tG",
	"W" => "A\tT",
	"S" => "C\tG",
	"Y" => "C\tT",
	"K" => "G\tT"
);
	#look up iupac variant

	$var_str =~ s/(\-|\+){1}\w+//g; # remove indels in variant string
	$var_str =~ s/\^.{1}//g; 	# remove ascii character after ^ that designates a quality value 
	 
	my @count1;
	my @count2;
	# iupac coded variant
	unless ($a2 =~ /[AGTC]/) {
		my $tmp= $hash{$a2}{$a1};
			
		#multi allelic site
		if (!$tmp) {
			my ($b1, $b2) = split/\t/,$iupac{$a2};
			(@count1) = $var_str =~ /$b1/ig;
			(@count2) = $var_str =~ /$b2/ig;
			#$c1 = $var_str =~ s/$b1/$1/ig;
			#$c2 = $var_str =~ s/$b2/$1/ig;
			my $c1 = scalar @count1;
			my $c2 = scalar @count2;
			if ($c1 > $c2) {
				$a2 = $b1;
			} elsif ( $c1 < $c2) {
				$a2 = $b2;
			} else {
				#print "same no. of allles!\n";
				#exit;
				$a2 = $b1; #pick a random allele
			}
				
		} else {
			$a2 = $tmp;
		}
	} else {
		$a1 = $a2;
	}

	return ($a1, $a2);
}

sub getPileupVarCov {
		# ref, var, varstring
	my ($self, $a1, $a2, $var_str)= @_;
	my %hash = ();
$hash{"M"}{"A"} = "C";
$hash{"S"}{"G"} = "C";
$hash{"W"}{"A"} = "T";
$hash{"B"}{"G"} = "T";
$hash{"R"}{"G"} = "A";
$hash{"Y"}{"T"} = "C";
$hash{"K"}{"G"} = "T";
$hash{"M"}{"C"} = "A";
$hash{"S"}{"C"} = "G";
$hash{"W"}{"T"} = "A";
$hash{"B"}{"T"} = "G";
$hash{"R"}{"A"} = "G";
$hash{"Y"}{"C"} = "T";
$hash{"N"}{"T"} = "N";
$hash{"N"}{"C"} = "N";
$hash{"N"}{"G"} = "N";
$hash{"N"}{"A"} = "N";

my %iupac = (
	"M" => "A\tC",
	"R" => "A\tG",
	"W" => "A\tT",
	"S" => "C\tG",
	"Y" => "C\tT",
	"K" => "G\tT"
);
	#look up iupac variant

	$var_str =~ s/(\-|\+){1}\w+//g; # remove indels in variant string
	$var_str =~ s/\^.{1}//g; 	# remove ascii character after ^ that designates a quality value 
	 
	my @count1;
	my @count2;
	my $cov = 0;
	
	# iupac coded variant
	unless ($a2 =~ /[AGTC]/) {
		my $tmp= $hash{$a2}{$a1};
			
		#multi allelic site
		if (!$tmp) {
			my ($b1, $b2) = split/\t/,$iupac{$a2};
			(@count1) = $var_str =~ /$b1/ig;
			(@count2) = $var_str =~ /$b2/ig;
			#$c1 = $var_str =~ s/$b1/$1/ig;
			#$c2 = $var_str =~ s/$b2/$1/ig;
			my $c1 = scalar @count1;
			my $c2 = scalar @count2;
			if ($c1 > $c2) {
				$a2 = $b1;
			} elsif ( $c1 < $c2) {
				$a2 = $b2;
			} else {
				#print "same no. of allles!\n";
				#exit;
				$a2 = $b1; #pick a random allele
			}
				
		} else {
			$a2 = $tmp;
		}
	} else {
		$a1 = $a2;
	}
	
	(my @tmp) = $var_str =~ /$a2/ig;
	$cov = scalar @tmp; 
	return $cov;	
}

sub getBamCov {
	my ($self, $chr, $position, $bam) = @_;
	my $cov = `samtools view -F 0X0404 $bam $chr:$position-$position | wc -l`;
	chomp $cov;
	return $cov;
}

sub getPileupRefCov {
	# ref, var, varstring
	my ($self, $var_str)= @_;
	$var_str =~ s/(\-|\+){1}\w+//g; # remove indels in variant string
	$var_str =~ s/\^.{1}//g; 	# remove ascii character after ^ that designates a quality value 
	my @count1;
	my $cov = 0;
	(@count1) = $var_str =~ /([.,])/ig;		
	$cov = scalar @count1;
	return $cov;	
}
sub getPileupStringAlleleCov {
  
  my ($self, $var, $var_str)= @_;

  $var_str =~ s/(\-|\+){1}\w+//g; # remove indels in variant string
  $var_str =~ s/\^.{1}//g;  # remove ascii character after ^ that designates a quality value
  my @count1;
  my $cov = 0;
  (@count1) = $var_str =~ /([$var])/ig;
  $cov = scalar @count1;
  return $cov;	
}

sub getATCGCovByBamAndCoord {
  my ($self, $chr, $coord, $ref, $var, $bam, $fa) = @_;
  my $varcov = my $totcov = my $af = 0;
  my $a_cov = my $t_cov = my $c_cov = my $g_cov = 0;
  my $coord2 = $coord;
  my $samtools = "/stornext/snfs6/cancer-analysis/code/src/samtools-0.1.8/samtools";
  my $varString; my $baseString;

  my $pileup = `$samtools view -F 0X0404 -u $bam $chr:$coord-$coord2 | $samtools pileup -c -f $fa - | grep -P "\t$coord\t" | cut -f8-| tail -1`;
  #print "$samtools view -F 0X0404 -u $bam $chr:$coord-$coord2 | $samtools pileup -c -f $fa - | grep -P \"\t$coord\t\" | cut -f8-| tail -1\n";
  if ($pileup) {
    #print "pileup = $pileup\n";
    ($totcov, $varString, $baseString) = split(/\t/, $pileup);
    #print "$totcov, $varString, $baseString\n"; 
    $varcov = $self->getPileupStringAlleleCov($var, $varString);
    $a_cov = $self->getPileupStringAlleleCov("A", $varString);
    $t_cov = $self->getPileupStringAlleleCov("T", $varString);
    $c_cov = $self->getPileupStringAlleleCov("C", $varString);
    $g_cov = $self->getPileupStringAlleleCov("G", $varString);
    $af = $varcov/$totcov if ($totcov > 0);
    $af = sprintf("%.3f", $af);
  }
  return "$totcov\t$varcov\t$af\t$a_cov\t$t_cov\t$c_cov\t$g_cov";
}
sub getATCGQualByBamAndCoord {
  my ($self, $chr, $coord, $ref, $var, $bam, $fa) = @_;
  my $varqual = my $a_qual = my $t_qual = my $c_qual = my $g_qual = 0;
  my $coord2 = $coord;
  my $samtools = "/stornext/snfs6/cancer-analysis/code/src/samtools-0.1.8/samtools";
  my $varString; my $baseString;
  my $pileup = `$samtools view -F 0X0404 -u $bam $chr:$coord-$coord2 | $samtools pileup -c -f $fa - | grep -P "\t$coord\t" | cut -f9-| tail -1`;
  #print "$samtools view -F 0X0404 -u $bam $chr:$coord-$coord2 | $samtools pileup -c -f $fa - | grep -P \"\t$coord\t\" | cut -f8-| tail -1\n";
  if ($pileup) {
    chomp $pileup;
    #print "pileup = $pileup\n";
    ($varString, $baseString) = split(/\t/, $pileup);
    #print "$varString, $baseString\n";
    $varqual = $self->getPileupStringBaseQual($var, $varString, $baseString);
    $a_qual = $self->getPileupStringBaseQual("A", $varString, $baseString);
    $t_qual = $self->getPileupStringBaseQual("T", $varString, $baseString);
    $c_qual = $self->getPileupStringBaseQual("C", $varString, $baseString);
    $g_qual = $self->getPileupStringBaseQual("G", $varString, $baseString);
  }
  return "$varqual\t$a_qual\t$t_qual\t$c_qual\t$g_qual";
}

sub getPileupStringBaseQual {
  my ($self, $var, $var_str, $base_str) = @_;
  $var_str =~ s/(\-|\+){1}\w+//g; # remove indels in variant string
  $var_str =~ s/\^.{1}//g;  # remove ascii character after ^ that designates a quality value
  $var_str =~ s/\$//g; # remove end of character marker
  my @varIndices;
  my $avgVarQual = -1;

  # indices of varbases if any
  while( $var_str =~ /([$var])/ig) {
    push(@varIndices, $-[0]);
  }

  if (@varIndices) {
    my @baseQuals = split(//, $base_str);
    # get base qual base on indices
    my @varQuals = @baseQuals[@varIndices];
    my @sangerVarQuals = map { ord($_)-33 } @varQuals;
    $avgVarQual = sprintf("%.3f", sum(@sangerVarQuals)/@sangerVarQuals);
  }
  return $avgVarQual;
}

sub getPileupCovByBamAndCoord {
	my ($self, $chr, $coord, $ref, $var, $bam, $fa) = @_;
	my $varcov = my $totcov = 0;
	my $coord2 = $coord;
my $samtools = "samtools";	
	# indel
	# samtools report indels 1 bp before
	
	my $varstr = '';
	($totcov, $varstr) =split(/\t/, `samtools view -F 0X0404 -u $bam $chr:$coord-$coord | samtools mpileup -f $fa - | grep -P "$chr\t$coord\t" | cut -f4,5`);  
	
	if ($varstr) {
		chomp $varstr;
		$varcov = $self->getPileupVarCov($ref, $var, $varstr);
	} else {
		$totcov = 0;
	}

	return ($totcov, $varcov);
}
	
sub revcomp {
	my ($self, $string) = @_;
	my $reverse = reverse($string);
	$reverse =~ tr/ACGTNacgtn-/TGCANtgcan-/;
	return $reverse;

}

sub getRefBases {
	my ($self, $chrom, $start, $end, $fasta) = @_;
	
	# reference fasta paths
	unless (-e $fasta) {
		die "ERROR: Couldn't locate reference fasta - $fasta";
	}
	
	#look up the reference base before the indel event
	my @tmp = `samtools faidx $fasta $chrom:$start-$end`;

	if(exists $tmp[1]) {
		my $bases;
      		shift @tmp; # skip the header (>chrX:xxx-xxx)
     	 	foreach my $line (@tmp) {
                	chomp $line;
        		$bases .= $line;
      		}
		return $bases;
	} else {
		die "ERROR: Out of range - samtools faidx $fasta $chrom:$start-$end\n";
	}
}

# get tcga uuid
sub getTcgaUuid {
  my ($self, $sample) = @_;
  my $string = `python $CODEDIR/submission-tools/fetch-tcga-uuids $sample | cut -f1`;
  chomp $string;
  return $string;
}

