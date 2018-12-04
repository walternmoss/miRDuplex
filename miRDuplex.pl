#!/usr/bin/perl

# miRDuplex takes a set of targt RNAs in FASTA format (all sequences on single lines) 
# and a list of miRNAs. It compares all potential interactions between the RNA and miRNA 
# seed complementarity and uses the program RNAup to predict potential duplexes. To run
# this program you need to have the Vienna RNA Structure package installed. 
#
# Usage: 
#
# $ perl HTP_miRDuplex.pl target.fasta miRNAs.fasta -options > output.txt
#
# the options right now are to use loose (-l; allows GU pairs) matching, otherwise
# strict (no GU pairs) matching is enforced. 

$targetfile = $ARGV[0];
open (INFILE1, "$targetfile") || die "Can't open the infile1!\n";
my @TargetFasta = <INFILE1>;
close(INFILE1) || die "can't close file";

$miRFile = $ARGV[1];
open (INFILE2, "$miRFile") || die "Can't open the infile2!\n";
my @MirFasta = <INFILE2>;
close(INFILE2) || die "can't close file";

my $Option = $ARGV[2];

#Put Target seqs into a Hash
my %FASTA = ();

#Loop through target file and put in hash
for (my $i = 0; $i < @TargetFasta; $i += 2) {

    $SeqName = $TargetFasta[$i];
    $SeqName =~ s/\s/_/g;
    $SeqName =~ s/>//g;
    chomp $SeqName;
	
    $Seq = $TargetFasta[$i + 1];
    chomp $Seq;
    $Seq =~ s/\R//g;	

    $FASTA{$SeqName} = $Seq;
}

my @TargetNames = keys(%FASTA);

#Loop through miRNA fasta
for (my $j = 0; $j < @MirFasta; $j += 2) { 

    my $miRNAname = $MirFasta[$j];
    chomp $miRNAname;
    $miRNAname =~ s/\s+/_/g;
    $miRNAname =~ s/>//g;
    chomp $miRNAname;

    my $miRNAseq = $MirFasta[$j + 1];
    chomp $miRNAseq;
    $miRNAseq =~ s/T/U/g;

    #Remove newlines from all systems
    $miRNAseq =~ s/\R//g;
	
    #Get seed sequence for miRNA
    my $Seed = substr($miRNAseq, 1, 6);
    my $RCSeed = RevComp ($Seed);

	
	#Test for loose option
    if (defined $Option) {$RCSeed =~ s/C/(C|U)/g;}
	
    #Loop through all targets
    foreach my $TargetName (@TargetNames) {
        chomp $TargetName;   
        my $TargetSeq = $FASTA{$TargetName};
        chomp $TargetSeq;
        $TargetSeq =~ s/\s+//g;
        $TargetSeq =~ s/-//g;
        $TargetSeq = uc $TargetSeq;
		
		#Get the location of seed matches			
		my @MatchCoords = ();			
		while ($TargetSeq =~ m/$RCSeed/g) {push (@MatchCoords,$-[0]);}
		my $Hits = @MatchCoords;
        if ($Hits == 0) {next;}		
		#if ($TargetSeq !~ m/$RCSeed/g) {next;}
	
	    print "$TargetName\t$TargetSeq\t$miRNAname\t$miRNAseq\t";
	    my $Input = "$miRNAseq&$TargetSeq";
	
	    my ($Out1, $Out2) = `echo "$Input" | RNAup -o -b`;
		chomp $Out1;
		chomp $Out2;
		$Out1 =~ s/\R//g;
			
		#Print the RNAup results
	    print "$Out2\t$Out1\t";
						
		#Get the location of seed matches			
		my @MatchCoords = ();			
		while ($TargetSeq =~ m/$RCSeed/g) {push (@MatchCoords,$-[0]);}
		my $Hits = @MatchCoords;
			
		#Print out the number and locations of all seed matches.
		print "$Hits\t";
		foreach my $MatchCoord (@MatchCoords) {print "$MatchCoord\t";}
		print "\n";       	
	}
}


################################################################################

sub RevComp {

   my $InSeq = $_[0];
   
   my $RevSeq = reverse $InSeq;
   $RevSeq =~ tr/ACGU/UGCA/;

   return $RevSeq;
   
   }
   

 

