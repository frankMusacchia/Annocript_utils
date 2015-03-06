#Differential Expression Analysis scripts 
Author: Francesco Musacchia (2015)

This scripts has the purpose to both get a differential expression analysis and then perform enrichment of pathways and go terms by using the output from Annocript

Please create a folder with the following files:

- a file with the raw counts after the mapping of the reads against the transcriptome;
- the filtered output from Annocript
- a target file tab separated which describes the experiments you carried on.
	All the files in the same condition will be taken together and compared with all the other conditions.

```
	name	condition
	sample1A	A
	sample1B	B
	sample2A	A
	sample2B	B
```

 Where sample1A means the sample 1 in the condition A.

###Differential Expression Analysis (DE_analysis_generic.R)
Takes as input:
- the sample_file with the grouping of replicates
- the table with transcripts raw counts
- OPTIONALLY the filtered output from *Annocript* (if given it will print also the annotation)

Execute this script with the following command in linux/unix systems:

**R CMD BATCH --no-save --no-restore '--args input_folder  count_filtered organism\_name target_file minFDR annocript\_out' DE\_analysis\_generic.R**

Where:
- input_folder -> Folder where all is taken and output goes
- counts_filtered -> Raw reads counts file
- organism_name -> name of the organism you are analyzing
- target_file -> file tab separated as described before
- minFDR -> minimum False Discovery Rate for the difference among two transcripts to be significant
- annocript_out -> filtered Annocript output

It returns tabular files for each combination of samples.

###GO terms and Pathways enrichments (GO_analysis.R and PW_analysis.R)

Execute this script with the following command in linux/unix systems:

**R CMD BATCH --no-save --no-restore '--args input\_folder annocript\_filt\_out significant\_transcripts\_file organism\_name p-value> min\_transcr min\_respect\_to[ann|sign] go mapping' GO\_analysis_4.R**

Where:
- input_folder -> input folder where your differential analysis table and Annocript output are present
- annocript\_filt\_out -> name of the output from Annocript
- significant\_transcripts\_file -> the file with the significant transcripts (de analysis executed before)
- organism\_name -> the name of the organism you are using
- p-value -> Adjusted p-value cutoff to consider a class as significant (def: 0.1)
- min\_transcr -> The minimum number of transcripts associated to a GO class to take it into consideration for the analysis
- min\_respect\_to -> a string that: if it is "sign", the minimum number of transcripts is from the significant table  while if it is "ann" it is from the annocript output table
- go_map -> the file with the GO terms mapping (included in the folder)

Run the script for the enrichments of the pathways with the following command identic to the previous but without the go_map:

**R CMD BATCH --no-save --no-restore '--args  input_folder annocript_filt_out significant_transcripts_file organism_name p-value min_transcr min_respect_to[ann|sign]' PW\_analysis\_1.R**

###Regulation plots (Regulation_plot.pl)

Regulation_plot.pl is a script that permits to obtain a plot in R associated to a particular differential expression analysis you did with Annocript utility DE_analysis_generic.R. Specifically it will show the number of GO terms up and down regulated in the samples comparison


REQUIRES: it needs the that the script DE_analysis_generic.R is located in the same folder
HOW TO USE: Put this script in a folder you want together with DE_analysis_generic.R and open a terminal. Then execute:

**perl regulation_plot.pl ann_out de_table out_type organism_name go_term_association**

WHERE: 
- ann_out is the filtered Annocript output
- de_table is the output table from the script DE_analysis.R 
- out_type [go|pathways] specifies which kind of ontology you want to plot
- organism_name name of the organism you are using. Will be used for names of files
- go_terms_association: needed if out_type=go, you can choose go terms associated to proteins or domains [proteins|domains]
OUTPUT: you will find both text tables and R plots in a folder named expression_plots in the working folder you are using
