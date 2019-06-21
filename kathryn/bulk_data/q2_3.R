# UNFINISHED: Can't retrieve info for FASTA sequence by NCBI, need to do this manually
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
# homer2 is the motif finding executable, and it can choke down FASTA files if you want tomail avoid all the nonsense above.  Running the homer2 command will also give you access to other options for optimizing the motif finding process.  homer2 works by first specifying a command, and then the appropriate options:
# 
# homer2 <command> [options]
# i.e. homer2 denovo -i input.fa -b background.fa > outputfile.txt
# 
# To find instances of the output motifs, use "homer2 find".  To see other commands, just type "homer2".

## Apply functions
# https://bianalystblog.wordpress.com/2013/07/13/r-grouping-functions-sapply-vs-lapply-vs-apply-vs-tapply-vs-by-vs-aggregate-vs-plyr/
# http://www.datasciencemadesimple.com/apply-function-r/

## Important functions
franke_cis_data[,"AssessedAllele"]
franke_cis_data[,"OtherAllele"]
franke_cis_data[,"GeneSymbol"]
franke_cis_data[,"Gene"]
franke_cis_data[,c("GeneSymbol", "GeneChr", "GenePos")]

## NCBI Efetch References
# https://www.ncbi.nlm.nih.gov/books/NBK25499/ 
# https://dataguide.nlm.nih.gov/edirect/efetch.html
# https://dataguide.nlm.nih.gov/eutilities/utilities.html#efetch 
# https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=cancer&reldate=60&datetype=edat&retmax=100&usehistory=y

## Google Searches
# https://www.google.com/search?ei=BdIHXfPnH-vB_QaymL_wCA&q=match+ensg+id+to+pubmed+id+in+r&oq=match+ensg+id+to+pubmed+id+in+r&gs_l=psy-ab.3..0i71l8.2036.2455..2523...0.0..0.0.0.......0....1..gws-wiz.4mJaRAHm_sQ
# https://www.google.com/search?ei=sc8HXfCkHqW4ggeOvr_oBQ&q=entrez+eutils+efetch+fcgi&oq=entrez+eutils+efetch&gs_l=psy-ab.1.0.0i71l8.0.0..21540...0.0..0.0.0.......0......gws-wiz.6ddBfuOMt70 
# https://cran.r-project.org/web/packages/easyPubMed/easyPubMed.pdf
# https://www.google.com/search?q=link+ensg+ids+to+pubmed+id&oq=link+ensg+ids+to+pubmed+id&aqs=chrome..69i57j33.5203j0j1&sourceid=chrome&ie=UTF-8

## Figuring out NCBI Efetch
# https://askvoprosy.com/tegi/bioinformatics/golosov/stranitsa/14
# https://www.biostars.org/p/293965/
# https://combine-australia.github.io/2017-05-19-bioconductor-melbourne/data_structures.html
# http://ogee.medgenius.info/cancer/Burkitt's%20lymphoma

## NCBI Texts
# https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?report=objectonly
# https://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T._entrez_unique_identifiers_ui/?report=objectonly 

# Efetch Information
# https://www.biostars.org/p/75700/
# https://github.com/gschofl/reutils/issues/1 
# https://rdrr.io/cran/reutils/man/efetch.html
library(reutils)

## Example: this uses HG38 - how to find hg19?  is it worth it?  try manually searching only a few of them, because the loop yields weird information ($gene_nuccore_pos[1] is very specific and may not be accurate.)
# https://www.ncbi.nlm.nih.gov/nuccore/NC_000020.11?report=fasta&from=50934855&to=50958564&strand=true
# efetch(uid="ENSG00000000419",db="gene") --> not that useful but works https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?id=ENSG00000000419&db=gene
# esearch(term="ENSG00000000419", db="gene",rettype="uilist") --> get uid  
# Not useful: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=ENSG00000000419&usehistory=y&WebEnv=web1&retmode=text&rettype=ID
franke_cis_data_unique_genes <-  franke_cis_data[!duplicated(franke_cis_data[,"Gene"]),"Gene"]
for (i in 1:nrow(franke_cis_data_unique_genes)){
  x <- esearch(term=franke_cis_data_unique_genes[i,"Gene"], db="gene")
  # elink(uid="8813", dbFrom="gene", dbTo="nuccore")
  y <- linkset(elink(uid=x, dbFrom="gene", dbTo="nuccore"))$gene_nuccore_pos[1]
  z <- efetch(uid=y, db="nuccore", retmode="text", rettype="fasta", strand=TRUE, retstart = NULL, retmax = NULL, seqstart = NULL, seqstop = NULL, outfile=paste("sequence",i,"nuccore",franke_cis_data_unique_genes[i,"Gene"], ".fasta", sep=""))
}

# Method 2: THIS USES HG19
# look up transcription start site
# if it's on the minus strand: +1000
# if it's on the plus strand: -1000

library("Biostrings")
# https://stackoverflow.com/questions/21263636/read-fasta-into-a-dataframe-and-extract-subsequences-of-fasta-file

# Other method for matching
# https://www.google.com/search?q=match+gene+name+to+sequence+r&oq=match+gene+name+to+sequence+r&aqs=chrome..69i57j33l5.5769j0j1&sourceid=chrome&ie=UTF-8
# https://www.bioconductor.org/packages/release/workflows/vignettes/annotation/inst/doc/Annotation_Resources.html
# library(rentrez)
# library(seqinr)
# COI <- entrez_fetch(db = "nucleotide", id = 167843256, file_format = "fasta")
# coi.fa <- read.fasta(file = textConnection(COI), as.string = T)

# https://www.nature.com/articles/nmeth.3799
# file://localhost/Users/kathryntsai/Zotero/storage/AC35NIIW/nrg.2017.html

# Links 6/17/19
# https://www.ncbi.nlm.nih.gov/nuccore?term=ENSG00000000419
# https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=ENSG00000000419
# https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?id=ENSG00000000419&db=gene
# https://www.ncbi.nlm.nih.gov/gene/8813
# https://www.ncbi.nlm.nih.gov/nuccore/NC_000020.11?report=fasta 
# https://github.com/gschofl/reutils/issues/1 - Repeat
# https://www.ncbi.nlm.nih.gov/books/NBK25499/ - Repeat

# Links 6/17/19 #2
# https://github.com/immunogenomics/IMPACT/blob/master/Features/Tcell_features/readme.txt.gz
# https://www.ncbi.nlm.nih.gov/gene/?term=2519
# https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.25/
# https://www.ncbi.nlm.nih.gov/nuccore?LinkName=assembly_nuccore_insdc&from_uid=37871
# https://www.ncbi.nlm.nih.gov/nuccore?LinkName=assembly_nuccore_refseq&from_uid=37871
# https://www.ncbi.nlm.nih.gov/nuccore/NC_000006.12?report=fasta&from=143494812&to=143511720&strand=true
# https://www.ncbi.nlm.nih.gov/nuccore/NC_000006.12?report=fasta&from=143494812&to=143511720&strand=true