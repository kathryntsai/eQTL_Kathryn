# ==========================================================
# QUESTION 2_3
# ==========================================================
#Search for motif in sequence with swapped alleles (according to cis eQTL ref/alt).: swap these? franke_cis_data[,"AssessedAllele"], franke_cis_data[,"OtherAllele"]

# http://www.eqtlgen.org/cis-eqtls.html
# This README accompanies the files with cis-eQTL results from eQTLGen
# 
# Files
# -----
#   File with full cis-eQTL results: cis-eQTLs_full_20180905.txt.gz
# File with significant (FDR<0.05) cis-eQTL results: cis-eQTL_significant_20181017.txt.gz
# 
# Column Names
# ------------
# Pvalue - P-value
# SNP - SNP rs ID
# SNPChr - SNP chromosome
# SNPPos - SNP position
# Zscore - Z-score
# AssessedAllele - Assessed allele, the Z-score refers to this allele
# OtherAllele - Not assessed allele
# Gene - ENSG name (Ensembl v71) of the eQTL gene
# GeneSymbol - HGNC name of the gene
# GeneChr - Gene chromosome
# GenePos - Centre of gene position
# NrCohorts - Total number of cohorts where this SNP-gene combination was tested
# NrSamples - Total number of samples where this SNP-gene combination was tested
# FDR - False discovery rate estimated based on permutations
# 
# Additional information
# ----------------------
#   These files contain all cis-eQTL results from eQTLGen, accompanying the article.
# 19,960 genes that showed expression in blood were tested.
# Every SNP-gene combination with a distance <1Mb from the center of the gene and  tested in at least 2 cohorts was included.
# Associations where SNP/proxy positioned in Illumina probe were not removed from combined analysis.


# Motif Finding with HOMER from FASTA files
# Most of HOMER's functionality is built around either promoter or genomic position based analysis, and aims to manage the sequence manipulation, hiding it from the user.  However, if you have some sequences that you would like HOMER to analyze, the program findMotifs.pl accepts FASTA formatted files for analysis.  Alternatively you could use the homer2 executable which also accepts FASTA files as input.
# 
# HOMER is designed to analyze high-throughput data using differential motif discovery, which means that it is HIGHLY recommended that you have both target and background sequences, and in each case you should have several (preferably thousands) of sequences in each set that are roughly the same length.  If you absolutely can't think of the proper background, homer will scramble your input sequences for you (starting v4.3, you can also call the scrambling script directly: scrambleFasta.pl).
# 
# A quick note about FASTA files - Each sequence should have a unique identifier.  In theory, HOMER should be flexible with what is in the header line, but if you're having trouble please just keep it simple with minimal quite-space, especially tabs.  For example:
# 
# >NM_003456
# AAGGCCTGAGATAGCTAGAGCTGAGAGTTTTCCACACG
# 
# Running findMotifs.pl with FASTA files:
# To find motifs from FASTA files, run findMotifs.pl with the target sequence FASTA file as the first command-line argument, and use the option "-fasta <file>" to specify the background FASTA file.  You should make every attempt to get sequences that represent a thoughtful background file - it would defeat the purpose of differential motif finding not to have it!
# 
# findMotifs.pl <targetSequences.fa> fasta <output directory> -fasta <background.fa> [options]
# 
# NOTE: you must choose an "organism" for the 2nd argument to keep with the structure of the command, even though this isn't actually relevant for FASTA based analysis.  Organism doesn't have to match the data in the FASTA files.  You can use a valid organism or just put "fasta" as a place holder. i.e.:
# 
# findMotifs.pl chuckNorrisGenes.fa human analysis_output/ -fasta normalHumanGenes.fa
# 
#  Many other options are available to control motif finding parameters.  findMotifs.pl will perform GC normalization and autonormalization be default (see here for more details).
# 
# Finding instances of motifs with FASTA files:
# To find instance of a motif, run the same command used for motif discovery above but add the option "-find <motif file>".  Motif results will be sent to stdout, so to capture the results in a file Add "> outputfile" to the end of the command.
# 
# findMotifs.pl <targetSequences.fa> fasta <output directory> -fasta <background.fa> [options] -find motif1.motif > outputfile.txt
# 
# For more information on the output file format, see here.
# 
# Using homer2 directly with FASTA files:
# homer2 is the motif finding executable, and it can choke down FASTA files if you want to avoid all the nonsense above.  Running the homer2 command will also give you access to other options for optimizing the motif finding process.  homer2 works by first specifying a command, and then the appropriate options:
# 
# homer2 <command> [options]
# i.e. homer2 denovo -i input.fa -b background.fa > outputfile.txt
# 
# To find instances of the output motifs, use "homer2 find".  To see other commands, just type "homer2".

franke_cis_data[,"AssessedAllele"]
franke_cis_data[,"OtherAllele"]
franke_cis_data[,"GeneSymbol"]
franke_cis_data[,"Gene"]

# https://www.biostars.org/p/75700/
tmp = tempfile()
efetch(franke_cis_data[,"Gene"], db="nucleotide", retmode="text", rettype="fasta", destfile=tmp)
readDNAStringSet(tmp)

# https://www.google.com/search?q=match+gene+name+to+sequence+r&oq=match+gene+name+to+sequence+r&aqs=chrome..69i57j33l5.5769j0j1&sourceid=chrome&ie=UTF-8
# https://www.bioconductor.org/packages/release/workflows/vignettes/annotation/inst/doc/Annotation_Resources.html
library(rentrez)
library(seqinr)
COI <- entrez_fetch(db = "nucleotide", id = 167843256, file_format = "fasta")
coi.fa <- read.fasta(file = textConnection(COI), as.string = T)
