#!/bin/perl
#Regulation_plot.pl is a script that permits to obtain a plot in R
#associated to a particular differential expression analysis you did 
#with Annocript utility DE_analysis_generic.R. Specifically it will show the
#number of GO terms up and down regulated in the samples comparison
#
# Author: Francesco Musacchia PhD (Copyright 2015)


#REQUIRES: it needs the that the script DE_analysis_generic.R is located in the same
#				folder
#HOW TO USE: Put this script in a folder you want together with DE_analysis_generic.R
#						and open a terminal. Then execute:
# 				perl regulation_plot.pl ann_out de_table out_type organism_name go_term_association
#WHERE: 
#	ann_out is the filtered Annocript output
# de_table is the output table from the script DE_analysis.R 
# out_type [go|pathways] specifies which kind of ontology you want to plot
# organism_name name of the organism you are using. Will be used for names of files
# go_terms_association: needed if out_type=go, you can choose go terms associated 
#												to proteins or domains [proteins|domains]
#OUTPUT: you will find both text tables and R plots in a folder named 
#				 expression_plots in the working folder you are using


use strict;
use warnings;
use Data::Dumper;
use Cwd;
#for each row of the Annocript output table I read the
#GO terms associated to the SwissProt id and concatenate the transcript name to a list
#of transcript names in a hash.
#No transcript can be repeated in the lists, thus, any control is unnecessary

#INPUT
#Annocript filtered output table
my $annOut = $ARGV[0];
#Complete table with the t-test executed. This script will decide which genes
#are differentially expressed
my $de_table = $ARGV[1];

#GO terms or pathways have to be plot?
my $dataType = $ARGV[2];

#Organism name
my $organism = $ARGV[3];

#GO terms should be associated to 
my $goTermAss = $ARGV[4];

#Output folder
my $outputFolder = "expression_plots";
unless (-e $outputFolder) { mkdir $outputFolder or die "Cannot create folder $outputFolder check permissions...\n";}


#INDEXES to acceed and SUFFIXES TO add
#Decide which indexes to use in the table depending by the GO terms wanted
my $indexes;
#Suffixes to save and search files
my $suffixes;

if ($dataType  eq 'go'){
	$suffixes->{'index1'} = "bp";
	$suffixes->{'index2'} = "mf";
	$suffixes->{'index3'} = "cc";
	if ( $goTermAss  eq 'proteins'){
		$indexes->{'index1'} = 24;
		$indexes->{'index2'} = 26;
		$indexes->{'index3'} = 28;
	}elsif ( $goTermAss  eq 'domains'){
			$indexes->{'index1'} = 35;
			$indexes->{'index2'} = 37;
			$indexes->{'index3'} = 39;
		}
}elsif ($dataType  eq 'pathways'){
		$indexes->{'index1'} = 12;
		$indexes->{'index2'} = 14;
		$indexes->{'index3'} = 16;
		$suffixes->{'index1'} = 'pwl1';
		$suffixes->{'index2'} = 'pwl2';
		$suffixes->{'index3'} = 'pwl3';
	}
$indexes->{'tr_id_index'} = 1;
		

#Hashes to be used to store the informations belonging to 
#the GO terms
my $bpHash;
my $mfHash;
my $ccHash;

##Where should I find the information? Which columns?
#my $bpIdIndex = 22;
#my $mfIdIndex = 24;
#my $ccIdIndex = 26;

#Indexes in the table of differentially expressed
my $trNameIndex = 0;
my $fcIndex = 1;
my $fdrIndex = 4;
#What is used to separate the GO terms in the Annocript output?
my $separator = "]---[";

#I will use the same parameters always used to decide if 
#a transcript is differentially expressed or not
my $minFDR = 0.05;
my $logFCThr = 2;

#Top to show in the table
my $topN = 40;

#Max length of the GO descriptions
my $maxLengthDescs = 50;


#Filling the hash with all the differential expression of the transcripts 
print "Filling the hash with all the differential expression of the transcripts \n";
open (DE,"<$de_table");
my $de_hash;
while (my $row_de = <DE>){
		next if $row_de =~ /logFC/;#Jumps the header
		chomp ($row_de);
		my @parts = split("\t",$row_de);	
		$de_hash->{$parts[$trNameIndex]}->{'fc'}= $parts[$fcIndex];#Get the fold change
		$de_hash->{$parts[$trNameIndex]}->{'fdr'}= $parts[$fdrIndex];#Get the FDR
}
#print Dumper $de_hash;
close(DE);

#Function to compute the log in base 2
sub log2 {
my $n = shift;
return log($n)/log(2);
}

print "For each row of the Annocript output, takes the transcript expression and increases\n";

#For each row of the Annocript output, takes the transcript expression and increases
#a counter in the hash for GO classesBP, MF and CC for the associated upregulated, downregulated
#or not differentially expressed genes
open(ANN, "<$annOut") or die "Cannot open $annOut";
while (my $row = <ANN>){
	chomp($row);
	next if $row =~ /TranscriptName/;
  my @parts = split("\t",$row);
  my $deType = '';#sring for up,down or not changing regulation
    
  #Define the change or not in expression of this transcript
  #To be defined as upregulated or downregulated the FDR should be 
  #very low and the Fold Change will decide if it is up-down
  #otherwise it is NDE
  my $transcript = $parts[0];
  
  #To be up or downregulated the transcript should have
  #a given minimum False Discovery Ratea and
  #depending by the fold change sign, if it is greater or lower of a given
  #threshold it is up or down regulated
  if ($de_hash->{$transcript}->{'fdr'} <= $minFDR){
		if ( $de_hash->{$transcript}->{'fc'} >= log2($logFCThr) ){
			$deType = 'up';
		}elsif ($de_hash->{$transcript}->{'fc'} <= -(log2($logFCThr))){
			$deType = 'down';
		 }else{
			 $deType = 'nde';
			} 
	}else{
		#print "Analyzing: $transcript\n";
		#print "FDR: ".$de_hash->{$transcript}->{'fdr'}."\tFC: ".$de_hash->{$transcript}->{'fc'}."\n";
		$deType = 'nde';
	}
	
	
	# Use these transcripts with the hash of the differentially expressed 
	# and get the fold change. Depending by the fold change, sum 1 in the hash 
	# of the GO term assigned to upregulated, down regulated or NDE (if it is 0) 
	# (here we should see if this value is very little and retain zero.  
	#print "GOs: ".$parts[$bpIdIndex]."\n";
	if ( $parts[$indexes->{'index1'}] ne '-'){
		my @goTermsBP = split (/\Q$separator\E/,$parts[$indexes->{'index1'}]);
		my @goDescs = split (/\Q$separator\E/,$parts[$indexes->{'index1'}+1]);
		#print $goDescs[0];
		my $GOpos = 0;
		foreach my $goTerm (@goTermsBP){
			#This variable is needed to count the number of transcripts which
			#have this GO term but it will not keep in count the 
			#transcripts which are not differentially expressed (NDE)
			#If them all are NDE, this count will be zero
			if ( not defined ($bpHash->{$goTerm}->{'count'}) ){
				if ($deType ne 'nde'){ 
					$bpHash->{$goTerm}->{'count'} = 1;
				}else{
					$bpHash->{$goTerm}->{'count'} = 0;
				}			
			}else{
				#To count only the total of differentially expressed genes
				$bpHash->{$goTerm}->{'count'}++ unless ($deType eq 'nde');
		  }
		  #reducing the description to be used in a plot in R
		  if (length($goDescs[$GOpos]) > $maxLengthDescs+10){
				my $reducedDesc = substr $goDescs[$GOpos], 0, $maxLengthDescs-10;
				$bpHash->{$goTerm}->{'desc'} = $reducedDesc." ($goTerm)";
			}else{
				my $reducedDesc = substr $goDescs[$GOpos], 0, $maxLengthDescs;
				$bpHash->{$goTerm}->{'desc'} = $reducedDesc;
		  }   
			#Sum this transcript to the count for this GO term
			if (not defined $bpHash->{$goTerm}->{$deType}){ 
				$bpHash->{$goTerm}->{$deType} = 1;
			}else{
				$bpHash->{$goTerm}->{$deType}++;
			}
			$GOpos++;
		}
  }
  if ( $parts[$indexes->{'index2'}] ne '-'){
		my @goTermsMF = split (/\Q$separator\E/,$parts[$indexes->{'index2'}]);
		my @goDescs = split (/\Q$separator\E/,$parts[$indexes->{'index2'}+1]);
		my $GOpos = 0;
		foreach my $goTerm (@goTermsMF){
			#Counting only differentially expressed terms
			if ( not defined $mfHash->{$goTerm}->{'count'} ){
				if ($deType ne 'nde'){ 
					$mfHash->{$goTerm}->{'count'} = 1;
				}else{
					$mfHash->{$goTerm}->{'count'} = 0;
				}
			}else{
				$mfHash->{$goTerm}->{'count'}++ unless ($deType eq 'nde');
		  }
		  #reducing the description to be used in a plot in R
		  if (length($goDescs[$GOpos]) > $maxLengthDescs+10){
				my $reducedDesc = substr $goDescs[$GOpos], 0, $maxLengthDescs-10;
				$mfHash->{$goTerm}->{'desc'} = $reducedDesc." ($goTerm)";
			}else{
				my $reducedDesc = substr $goDescs[$GOpos], 0, $maxLengthDescs;
				$mfHash->{$goTerm}->{'desc'} = $reducedDesc;
		  } 
			#Sum this transcript to the count for this GO term
			if (not defined $mfHash->{$goTerm}->{$deType}){ 
				$mfHash->{$goTerm}->{$deType} = 1;
			}else{
				$mfHash->{$goTerm}->{$deType}++;
			}
			$GOpos++;
		}
  }
  if ( $parts[$indexes->{'index3'}] ne '-'){
		#Splits the GO terms and descriptions using the separator
		my @goTermsCC = split (/\Q$separator\E/,$parts[$indexes->{'index3'}]);
		my @goDescs = split (/\Q$separator\E/,$parts[$indexes->{'index3'}+1]);
		
		#This is used to determine which is the position to fetch the description		
		my $GOpos = 0;
		#For each GO term sum this transcript to their count of transcripts
		#associated
		foreach my $goTerm (@goTermsCC){
			if ( not defined $ccHash->{$goTerm}->{'count'}){
				if ($deType ne 'nde'){ 
					$ccHash->{$goTerm}->{'count'} = 1;
				}else{
					$ccHash->{$goTerm}->{'count'} = 0;
				}
			}else{
				$ccHash->{$goTerm}->{'count'}++ unless ($deType eq 'nde');
		  }
		  #reducing the description to be used in a plot in R
		  if (length($goDescs[$GOpos]) > $maxLengthDescs+10){
				my $reducedDesc = substr $goDescs[$GOpos], 0, $maxLengthDescs-10;
				$ccHash->{$goTerm}->{'desc'} = $reducedDesc." ($goTerm)";
			}else{
				my $reducedDesc = substr $goDescs[$GOpos], 0, $maxLengthDescs;
				$ccHash->{$goTerm}->{'desc'} = $reducedDesc;
		 } 
			#Sum this transcript to the count for this GO term
			if (not defined $ccHash->{$goTerm}->{$deType}){ 
				$ccHash->{$goTerm}->{$deType} = 1;
			}else{
				$ccHash->{$goTerm}->{$deType}++;
			}
			$GOpos++;
		}
	}
} 
close(ANN);


print "Printing tables sorted..\n";
printTableSorted($outputFolder."/expressionTable_".$suffixes->{'index1'},"count",$bpHash);
printTableSorted($outputFolder."/expressionTable_".$suffixes->{'index2'},"count",$ccHash);
printTableSorted($outputFolder."/expressionTable_".$suffixes->{'index3'},"count",$mfHash);


#The previous code created new tables with the abundances of expression
#relatively to bp, mf and cc by using the same name as the expression tables and adding
#a suffix. This suffix is searched here to plot all these tables
print "Printing plots for biological processes...\n";
my $parameters =  " $outputFolder expressionTable_".$suffixes->{'index1'}." $organism none  $minFDR 3 $annOut ";
my $command = "R CMD BATCH  --no-save --no-restore '--args ".$parameters." ' DE_analysis_generic.R";		
print "Executing command: $command\n";
system($command);# or die "Unable to execute command: $command\n";

$parameters =  " $outputFolder expressionTable_".$suffixes->{'index2'}." $organism none  $minFDR 3 $annOut ";
$command = "R CMD BATCH  --no-save --no-restore '--args ".$parameters." ' DE_analysis_generic.R";		
print "Executing command: $command\n";
system($command);# or die "Unable to execute command: $command\n";
		
$parameters =  " $outputFolder expressionTable_".$suffixes->{'index3'}." $organism none  $minFDR 3 $annOut ";
$command = "R CMD BATCH  --no-save --no-restore '--args ".$parameters." ' DE_analysis_generic.R";		
print "Executing command: $command\n";
system($command);# or die "Unable to execute command: $command\n";		


#Prints a table with for each row a GO term and three columns represent
#the number of upregulated, downregulated, not differentially expressed
#transcripts associated to such term.
sub printTableSorted {
	my $expTable = shift;
	my $field = shift;
	my $hash = shift;

	open (EXP_T,">$expTable") or die "Unable to open $expTable";
	my $numKeys = 0;
	
	 # The hash is sorted putting before the GO terms which have more
	 # differentially expressed transcripts: the 'count' is given as the 
	 # sum of both the up and down regulated.
	 # GO terms build a final table with 
	 # each column being a GO term and three rows represent the 
	 # number of UpRegulated, DownRegulated and NDE genes
	 foreach my $key ( #
		 sort { $hash->{$b}->{$field} <=> $hash->{$a}->{$field} } #
		 keys %{$hash}
	 )
	 {
			#These variables could be not initialized if the particular differential expression
			#is not found. Thus the following is initialization
			if (not defined $hash->{$key}->{'up'}){ 
				$hash->{$key}->{'up'} = 0;
			}
			if (not defined $hash->{$key}->{'down'}){ 
				$hash->{$key}->{'down'} = 0;
			}
			if (not defined $hash->{$key}->{'nde'}){ 
				$hash->{$key}->{'nde'} = 0;
			}
			
			my $totCountDE = $hash->{$key}->{'count'};#Count of only the differentially expressed
			my $overallCount = $totCountDE + $hash->{$key}->{'nde'};#Overall count with also NDE
			#print "GO: $key Overall count: $overallCount - Total count diff expr: $totCountDE \n";
			#print $hash->{$key}->{'up'}." ".$hash->{$key}->{'down'}." ".$hash->{$key}->{'nde'}."\n";
			my $upPerc = 0;
			my $downPerc = 0;
			my $ndePerc = 0;
			
			#To avoid the division by zero I used this commented check
			#But this should never be zero
			#if ( $overallCount != 0){#DEBUGCODE
			 $upPerc = ($hash->{$key}->{'up'}/$overallCount)*100;
			 $downPerc = ($hash->{$key}->{'down'}/$overallCount)*100;
			 $ndePerc = ($hash->{$key}->{'nde'}/$overallCount)*100;
			#}
			
			my $desc = $hash->{$key}->{'desc'}." [$overallCount]";
			#$desc =~ s/ /_/g;
			
			print EXP_T join("\t",$desc,$upPerc,$ndePerc,$downPerc);
			#print EXP_T join("\t",$key,$upPerc,$ndePerc,$downPerc);
			print EXP_T "\n";
			$numKeys ++;
			last unless $numKeys < $topN;
	 } 
	close (EXP_T); 
}



#Using the final table to build the stacked plot in R
