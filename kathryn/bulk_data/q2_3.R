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

library(reutils)


# Method 2: THIS USES HG19
# look up transcription start site
# if it's on the minus strand: +1000
# if it's on the plus strand: -1000

# http://hgdownload.cse.ucsc.edu/downloads.html
# http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/ 
# http://genome.ucsc.edu/admin/git.html
# http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr21%3A33031597-33041570&hgsid=731188783_MLKAzbhB9nJUnx64kxdfc05oomHe
# franke_cis_data_unique_genes_subset <- franke_cis_data_unique_genes %>% sample_n(10)

a <- franke_cis_data[which(Gene==c('ENSG00000127191', 
                              'ENSG00000214279', 
                              'ENSG00000168765', 
                              'ENSG00000125966',
                              'ENSG00000258531',
                              'ENSG00000237437'
                              )), c('Gene', 'GeneSymbol', 'GeneChr', 'GenePos')]
a <- a[!duplicated(a),]


# http://grch37.ensembl.org/Homo_sapiens/Gene/Sequence?db=core;g=ENSG00000127191
# http://grch37.ensembl.org/Homo_sapiens/Gene/Sequence?db=core;g=ENSG00000214279
# http://grch37.ensembl.org/Homo_sapiens/Gene/Sequence?db=core;g=ENSG00000168765
# http://grch37.ensembl.org/Homo_sapiens/Gene/Sequence?db=core;g=ENSG00000125966
# http://grch37.ensembl.org/Homo_sapiens/Gene/Sequence?db=core;g=ENSG00000258531
# http://grch37.ensembl.org/Homo_sapiens/Gene/Sequence?db=core;g=ENSG00000237437


aa <- data.frame(a, rbind(
  # TRAF2
  '>chromosome:GRCh37:9:139775364:139821659:1
GTGACTGGTCCTGGCGGGTGCTTAGGAGGCGGGGGCAGAGACATGGCTCACTCCAGCCTCAACCTCCTGGGCTCAAGCAATTCTCCCACCTCAGCCACTCAGGTAGCTGGGGCTACAGGCGTGCACCATCACGCCCGGCTATTTTTTTTTTTTTTTTGAGACAGAATTTGGCTCTTATTGCCCAGGCTGGAGTGCAATGGTGCAATCTCAGCCACAACCTCTGCCTCCCAGGGTCAAGTGATTCTCCTGCCTCAGCCTCCCAAGTAGCTAGGATTACAGGTGTGCACCACCACGCCGGGCTAATTTTTGTATTTTTAGTAGAAACGGGGTTTCTCCATGTTGGTCAGGTTGGTCTCGAACTCCCAACCTCAGGTGATCTGCCTGCCTCAGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCACACCCGGCCTAATTTTTGCATTTTTTGGAGAGACGGAGGTTTCACTATATTGCCCAGGCTGGTCTTGAACTCCCGAGCTCAAGTGATCCAGCCTTGGCCTCCCAAAGTGCTGGGATTACAGGTGTGAATCTTCTATTACTTTCCTGTGTGTCCTTCCCAAGTTCACAGAAAACACAAGCTGACAGAAATCCATTCTTTCCCACGGTTGCTGGCACGGATGCTGCTGTGTGACAGCACAAGTTCATGTGCTTGAGGGTCTCCTGGCCACAGCCAGCCCCAGGCTTGCCCTGGACAGGCGGGCTTCAGATTCCCCCCAAGGTTCTATGCCATTGAGTTGGAGGACAAGAGCTTCCCCACTGTCCAGCGTTGAAGCAAACCTCCCTCCAGGTCCCACCCCTCCCCAGGCCTTTGTCCTTTTATCAGGGATTGCCCTGGGGGAGCCCCCCTGGGGAGCTGATGAGGAGGGAGCCTCTGTGCCCCGTGCACTGGCCTGAGGGGGTGAGGCCAGGGCTGGGGCTGCTGCAGGGGACTCTGGCCCGGAGCTGCGGCTGCCTTCCTGGAAACACTCCTGTC',

# RP11-108K14.4
'>chromosome:GRCh37:10:135266432:135337662:1
TGCTCTGTCACCCAGGCTGGAGTGCAGTGGTGCAATCTCAGCTCACTGCAACCTCCGCCTCTTGGGCTCACGCCATTCTCATGCCTCAGCCTCCCGAGTAGCTGGGACTACAGGTGCATGTCACCCTGCCCGGCTAATTTTTGTATTTTTAGTAGAGACAGGGTTTCACCATGTTGGCCAGGGTGGTCTCGAACTCCTGACCTCAGGTATCTGCCCACCTCAGCCTCCCAAAGTGCTAGGATTACAGGCGTGAGCCACTGCGCCTAGCCTGTTCCATTTTTATAACTTCTATTTCTTTGCTGAGATTTCAATTTTTCCACTTATTTCAAGAGCATTTGCAATTTGTTGTGGAAACATTTTTATGATGCCTGTTTTAAAATCCTTGTCAGATAATTCTGACACCTGACTCTTCCTGGTGCCGCCCTTAGTTGGTTGTCTCTTACTGAAATGGTCGTGTCCCTGGGTCTTCATATGATGAGTGATTCTTAACTGTATCTTCTCTCAAGGGTTTTCATTCCAACGACAGTTCAGTTCTCAGAGTTCTTGCAAAGCTCTTCCGTGTCATTCTTCTGGTGCTGCTGGGACTCCCACTCAAACCCTAATCCTACCCCTACCCCTAACTCTGCTAGTGCCATCTGCTGGCACAGAGGACGCCTCCTTGGCTGTCTGGGGCCACTGCCTTCTGTGGTGGTGGCATGGCTGAAGCCGACGGGCCCGGGTGGAGGGCAGAGACGCAGGGCTTCACTGATGCTGTGCTGTGGGCAGGTTGGACCATCTCTGCCTGATGGTTGGGTCTAGTCCCATCAGCAGCTTGATCAGTGTGGGAGGAAACCCACAGCCTGCACCCTCCCCATGTCTGGTGATGCTCAGGCCAGCAAGGCATTGGCTGCTCTGCTCACGGCTGAGAGGCAGAGCTGCTGATTCGCTCTGGAAAGGCTGTGGCTGTCCCCATTGGACCCTGGAGGTCCTTGGCAGCTCCCACTTCTGTATATGTGACTCT',

  # GSTM4
'>chromosome:GRCh37:1:110197703:110208718:1
GGGTGAATACTAATGTAACAATATCAATTCTGTTTCTTTCATAATTATTCTCATAAGTGTTAGCTTTCTTAAAAAATCATCTTCACAATTGGACCTAGACATTTTATCTTTCCCAGCTCAAGATGGATAAACTCAGAATGACAGGCCAGGCACCGTGGCTCATGCCTGTAATCCCAGCACTTTGGGAGGCGGAGGCGGGTGGATCACTTGAGACCAGGAGTTCAAGACCAGTCTGGTCAACACGGTGAAATCCCATCTCTACAGAAAATACAAAAAAGTTAGCAGGGCGTCGTGGCGCTCGTCTGTAATCCCAGCTGCTCGGGAGGCTGAGGCAGGAGAATCGCTTGAGCCCGGGAGGCGGAGGTTGCAATGAACCGAGATCGCCCTACTGCACTCCAGCCTGGGCAACACAGCGAGTCCCCGTTCCAAAAAAAAAAAAAAAAAAAAGACTCAGAATGACACATTCTTAATGTCTTAATGTATAATGAGAAGTATAACTATTCAGAAGTCCAAGCACTTGTGTGTGCATGGTGCTATCGTCTTGCTGTAGGCAAACCTACATGGGAATCCACCTTGACACATAGGCCACTTCCTGAGCCTGGACCCAGTCTCAGAGCTGGGGAACTGGCCCAATGCAAAAGGGTCGGGAGCATCTGCAACAGAGACTGAGCTCTATCAGCTTCGGTGACATAGCCTCCATTCACGCTCCCCAACTCAGCAGAGAGAGCACACCATCAGACTTCTAAGACTTAGTAGCCAAGAAGTGTTGAATTAAACTCTCTGAGACCTCTCTTTAGTCTGACCCTGGCAGCCTCAGTCTCCCAGAGCCTGTGGGAACTCGGCAGCCGAGAGGCAGAAGGCTGGGCGACGTCCGGAGAAGAAGAAACGGGGGAAGAACTTTTCTCTTACGATCTGGCTTTACTCTCACGCGCACAGCCGAGTCCCTGGGGACCCAGCAGAGGTCCGAAGCGGAGCGGGGCGGGGCGGGGCTACGGAAGCT',  
  
  # MMP24
  '>chromosome:GRCh37:20:33813457:33865404:1
TCAGGCACATCCTGTACACCCTCAACTTCCTCCCTGAAATGCCCTCTCCTCTCCCCTCCCCTGGTTTAAATCCAGCCTTTCCTTCAGGTCCTAGCTAGAGTCCCACCCTCCTCCATTACACTTTTTCCTAACTTCCCCAGCCCCTCATCAGGTCTCCCTCGACCCCTGGAGCAGCCGCCTGTATCCTTGTTGACACGTGGTGTTTGCCATTCTGTTTTTCTTTGCGTTTATGTTTTCTCATCAAGCTGTAAACCTCTGCACTCCAGGAAATCGGCCTTCTCTTTTTTGTTTTTGTTTTTGGGTGTTTTTTTTTTTCCTTTCAGCATCCCAAGTAGCCCAGAGCCCTACAAATAGAAAGCGCTTGATAAATATGGGTTGTTTTGGATGAAAGATGATCTTCCCAGCTGGATGAGCTCTTTCCCTCTAGGGGTGAGGGACTGAGCTTCACCTCAACTTCACCTAATGTCTCACAGGAATGCCCCCTCTTAAACGCTCCTCAGGTCTTCTCTACCCCTTTCCCCCACCACCCCGGCCCAGGAGAGCCACTCAGAGCTCCCTCTCCGGAGTCACCCTCTCAGACCCAGGGCTATTTCTTTAGCTTGCTCTCGAACCCCGCCTTTACCCCAGAGCCGCTCCTCAGTCTCCTCCCCCAGGCTGCCTCTCAGGCCCCCTCTCCTCCAATGCCACTCCCTAGGACCCTCTCCCCAGGGTTCCGGACCCCCAACCCCGGATATGTTCCCCCAGGGCCATCCAGTTGTTCTACCTCCAGGATTACTGATTTAGCCTCTTCCCCCTCCTGCTCAGGGCTACCCCCAGACCCCGACCTCCAGTTCCACCCCTTCGGACCCTGCCTCTAGAACTGTTCTTTTAGCCTCTGCAACCCCCTCCCCAGCGCCGCTCCGCGGCCCCTGCCCCAGGGCTGGTCGGGGAGCCACTGCCCTCGCAGGTCCCGCCCGGGTGCTGCCCCCGCGCCCCCGGCCGCCCCCCACCGGGGCGGGGCGCGCGGAGGCGGGGGCGCGC',  
  
# BANF1P1
  '>chromosome:GRCh37:14:69402319:69404181:1
TCTCAAACTCCTGACCTCAAGCGATCTGTCCACCTTGGCTTCCCAAAGTGCTGGGATTACAGGCATGAGCCACCACACCCGGCTAGAGCGCTCTATTCTGATGTTTAAGCCTCAACACACAATCCAGGCCATCCACTCTGGCTGACCTCCCCAGAGAAGCCTCCATTGACAAAGGTGCTTTTTCCATTTGTTGCACAGTGGAATAAACAAGATAACCCACCTGCAGCCAAACAAGACTAGCTGGGCAGTGGTGGGACAGCCCAGGCTCCAAACTCACGTGCTGGCCACACCTCTTCCAAGCACATCCTCCCACTCTCCCGGAATCAAAAATACCAGCCACCATTCAGGCCAGGGGAGTTAAGTGTTCAGAGGATGGGATCAGGAAACCTGGGACCTGCCCTGGCTCAGCCAGGCGTGTGACCTCAGATAAGTCACTTCAATCCCTTCCAGCTCTGCTAGCGCAAACAAGTCCCTCCTACCTTGCAGGCTGGAGCCCAGAAACTTTTCCCATCAGAACACGCCTAATGAGAAGACACAACCAAACAGCACCCGGCAACAATCTTTTGTTCCTGTGTAGGGTTTATTCTGATTGGTGGGCTCTCAAAGAAAGCCAAAAAATTTCCCAGGGACAATTTCATCTCCATGAACTCAAAACATACTCAGTGTTTCTGTCTCAAGACAGCCAGCTCCTAACCCTTCTCAGCAAGCAGAATCCTGACGTTAGTGTCCCGGACAGCATCTAAAAGCTTTATGCTAGAACATTCCTAGCAGGACACAGAAGGACCAGGATCATCGGCTGCCTCCAGGGAGAGAAATGGGAAGGCAGGGGGTAGAGAACACTTGCTTTTCACAGAATATGCTTTTGTTTGAGTTTTGTACCAAAATGTGCTGTCATTACCTTTTTTTTTTTTTTAAGTATGCAAATAAAAGAACAGTCTGTCTAGGGGGTGGCGGGAGGAACTGTTACGGGAATTGAAGCTGCCGATTAGGCCTAATCAAG',
  
  # ASS1P12
  '>chromosome:GRCh37:9:32944994:32947820:1
TCACTTGTTTGTCTTTTTTTGACAAATGTCTATTCAAATCCTTTGACCATTTTTTAATTGCATTATTTGTTTTCTTACAATTGAGTTGTTTCAGTTCCTTATTATTTTAGATATTAATCCCTTATCAGATGTATGGTTTGCAAATATTTTCTCCCATTCTGTAGGTTGTCTCTTAACTGTGTAGATTATTTCCTTGGCTGTGCAGAAGCTTTTTAGTTTGATACAATCCCATTTGTCTATTATTACTTTTGTTGCCTGTGCTTTTGGGTTTATATTTAAAAATTCAGGCCTGGGTACCGTGGCTCATGCCTCTAATCCTAACACTTCAGGAGGCCAAGGCAGGAGGATCTCTTGAGGCCAGGAGTTTGAGACCATCCTGGGCAACATAGCAAGGCCTCGTCTTTAAAAAAAATTAGCTGGGTATAGTGGTGCGTGCCTGTATTCCTTATTACTTGGGAAACTAAGGCAGGAGGATCCCTTGAGCCCAGGAGTATGAGGCTGCAGTACAGCTATGCTCCTGCCACTGCACTCCAGCCTGGGCAACAGAGCAAGACCCTGTCTCTGAAAATAAATAAATAAATAAATTTTTAAAAATCTTTAAATAGATTTTATGTTTCCAAGTCATTAAAGAACATGAGATGTCTTTCCATTTGTTCAGATAACATTTATAATTTTAATTGAGCTTTTATAACTTTAAATAGGTCTTATATCTTATTAATTTTACTGTCAAGGTGAAACTATTTTGATTTCTGAAATTTGCTTTAAAGTAGTACTACAGTAGGTAAGGTGGAGAGGCAGGATGAAACAAGAAGGGCAGAATATTGACGATTGTTGAAGCTGGTTATGGGTACTATTCTCTGGGTTTTTGTGGATGGTTTGGTATTTTGAAATTATACTCCCTGCTCTGTCGCCTGCCACCGCTCCCTGAGCCCGAGTGGTTCACCGTACCATGAAGACAGATGGCAGATGCCAGGAACTCGAGCCTCCAATCCCAGATGCT'  
))

colnames(aa)[5] <- "Fasta Sequence"

write.table(aa$"Fasta Sequence", 'q2_3_output/fasta_file1.fa', sep='\t', quote=F, row.names = F, col.names = F)

# Try with bioconductor

#===================

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

# findMotifs.pl fasta_file1.fa hg19 analysis_output/ -fasta normalHumanGenes.fa > output.txt
# 
# Selected Options:
#   Input file = fasta_file1.fa
# Promoter Set = human
# Output Directory = analysis_output/
#   Will use FASTA files for motif finding
# Target Sequences = fasta_file1.fa
# Background Sequences = normalHumanGenes.fa
# !!!!
#   Could not open file fasta_file1.fa
# !!!!
#   kat-juno:Github kathryntsai$ 
#   kat-juno:Github kathryntsai$ findMotifs.pl fasta_file1.fa human analysis_output/ -fasta normalHumanGenes.fa > output.txt
# 
# Selected Options:
#   Input file = fasta_file1.fa
# Promoter Set = human
# Output Directory = analysis_output/
#   Will use FASTA files for motif finding
# Target Sequences = fasta_file1.fa
# Background Sequences = normalHumanGenes.fa
# !!!!
#   Could not open file fasta_file1.fa
# !!!!
#   kat-juno:Github kathryntsai$ cd .. ..
# kat-juno:Documents kathryntsai$ cd ..
# kat-juno:~ kathryntsai$ cd homer
# kat-juno:homer kathryntsai$ findMotifs.pl fasta_file1.fa human analysis_output/ -fasta normalHumanGenes.fa > output.txt
# 
# Selected Options:
#   Input file = fasta_file1.fa
# Promoter Set = human
# Output Directory = analysis_output/
#   Will use FASTA files for motif finding
# Target Sequences = fasta_file1.fa
# Background Sequences = normalHumanGenes.fa
# Found mset for "human", will check against vertebrates motifs
# Parsing FASTA format files...
# Found 6 sequences
# !! 6 of 6 contained bad nucleotide characters [not ACGTN], replaced with N
# readline() on closed filehandle IN at /Users/kathryntsai/homer/bin/fasta2tab.pl line 16.
# Found 0 sequences
# Use of uninitialized value $len in numeric lt (<) at /Users/kathryntsai/homer/bin/cleanUpSequences.pl line 36, <IN> line 1.
# Use of uninitialized value $len in numeric gt (>) at /Users/kathryntsai/homer/bin/cleanUpSequences.pl line 40, <IN> line 1.
# Use of uninitialized value $line[1] in substitution (s///) at /Users/kathryntsai/homer/bin/cleanUpSequences.pl line 44, <IN> line 1.
# Use of uninitialized value $line[1] in substitution (s///) at /Users/kathryntsai/homer/bin/cleanUpSequences.pl line 45, <IN> line 1.
# Use of uninitialized value $line[1] in substitution (s///) at /Users/kathryntsai/homer/bin/cleanUpSequences.pl line 46, <IN> line 1.
# Use of uninitialized value $line[1] in substitution (s///) at /Users/kathryntsai/homer/bin/cleanUpSequences.pl line 47, <IN> line 1.
# Use of uninitialized value $line[1] in substitution (s///) at /Users/kathryntsai/homer/bin/cleanUpSequences.pl line 48, <IN> line 1.
# Use of uninitialized value $line[1] in pattern match (m//) at /Users/kathryntsai/homer/bin/cleanUpSequences.pl line 50, <IN> line 1.
# Use of uninitialized value $line[0] in concatenation (.) or string at /Users/kathryntsai/homer/bin/cleanUpSequences.pl line 55, <IN> line 1.
# Use of uninitialized value $line[1] in concatenation (.) or string at /Users/kathryntsai/homer/bin/cleanUpSequences.pl line 55, <IN> line 1.
# Use of uninitialized value in concatenation (.) or string at /Users/kathryntsai/homer/bin/cleanUpPeakFile.pl line 26, <IN> line 1.
# 
# Progress: Step4 - removing redundant promoters
# 
# Progress: Step5 - adjusting background sequences for GC/CpG content...
# 
# Sequences processed:
#   Auto detected maximum sequence length of 1034 bp
# 6 total
# 
# Frequency Bins: 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.6 0.7 0.8
# Freq	Bin	Count
# 0.4	4	1
# 0.5	6	1
# 0.6	7	4
# Bin	# Targets	# Background	Background Weight
# 
# Normalizing lower order oligos using homer2
# 
# Reading input files...
# 0 total sequences read
# Autonormalization: 1-mers (4 total)
# A	inf%	inf%	nan
# C	inf%	inf%	nan
# G	inf%	inf%	nan
# T	inf%	inf%	nan
# Autonormalization: 2-mers (16 total)
# AA	inf%	inf%	nan
# CA	inf%	inf%	nan
# GA	inf%	inf%	nan
# TA	inf%	inf%	nan
# AC	inf%	inf%	nan
# CC	inf%	inf%	nan
# GC	inf%	inf%	nan
# TC	inf%	inf%	nan
# AG	inf%	inf%	nan
# CG	inf%	inf%	nan
# GG	inf%	inf%	nan
# TG	inf%	inf%	nan
# AT	inf%	inf%	nan
# CT	inf%	inf%	nan
# GT	inf%	inf%	nan
# TT	inf%	inf%	nan
# Autonormalization: 3-mers (64 total)
# Normalization weights can be found in file: analysis_output//seq.autonorm.tsv
# Converging on autonormalization solution:
#   ...............................................................................
# Final normalization:	Autonormalization: 1-mers (4 total)
# A	inf%	inf%	nan
# C	inf%	inf%	nan
# G	inf%	inf%	nan
# T	inf%	inf%	nan
# Autonormalization: 2-mers (16 total)
# AA	inf%	inf%	nan
# CA	inf%	inf%	nan
# GA	inf%	inf%	nan
# TA	inf%	inf%	nan
# AC	inf%	inf%	nan
# CC	inf%	inf%	nan
# GC	inf%	inf%	nan
# TC	inf%	inf%	nan
# AG	inf%	inf%	nan
# CG	inf%	inf%	nan
# GG	inf%	inf%	nan
# TG	inf%	inf%	nan
# AT	inf%	inf%	nan
# CT	inf%	inf%	nan
# GT	inf%	inf%	nan
# TT	inf%	inf%	nan
# Autonormalization: 3-mers (64 total)
# 
# Progress: Step6 - Gene Ontology Enrichment Analysis
# Skipping...
# 
# Progress: Step7 - Known motif enrichment
# 
# Reading input files...
# 0 total sequences read
# 428 motifs loaded
# Cache length = 11180
# Using hypergeometric scoring
# Checking enrichment of 428 motif(s)
# |0%                                    50%                                  100%|
#   =================================================================================
#   Illegal division by zero at /Users/kathryntsai/homer/bin/findKnownMotifs.pl line 152.
# 
# Progress: Step8 - De novo motif finding (HOMER)
# 
# Scanning input files...
# !!! Something is wrong... are you sure you chose the right length for motif finding?
#   !!! i.e. also check your sequence file!!!
#   
#   Scanning input files...
# !!! Something is wrong... are you sure you chose the right length for motif finding?
#   !!! i.e. also check your sequence file!!!
#   
#   -blen automatically set to 2
# Scanning input files...
# !!! Something is wrong... are you sure you chose the right length for motif finding?
#   !!! i.e. also check your sequence file!!!
#   Use of uninitialized value in numeric gt (>) at /Users/kathryntsai/homer/bin/compareMotifs.pl line 1389.
# !!! Filtered out all motifs!!!
#   Job finished


# trying to figure out fasta sequences:
d <- read.fasta('q2_3_output/fasta_file1.fa')

# kat-juno:homer kathryntsai$ findMotifs.pl fasta_file1.fa hg19 analysis_output/ > output.txt
# Selected Options:
#   Input file = fasta_file1.fa
# Promoter Set = hg19
# Output Directory = analysis_output/
#   
#   !!! hg19 not found in /Users/kathryntsai/homer/.//config.txt
# Try typing "perl /Users/kathryntsai/homer/.//configureHomer.pl -list" to see available promoter sets
# If avaliable, type "perl /Users/kathryntsai/homer/.//configureHomer.pl -install hg19" to install

# kat-juno:homer kathryntsai$ findMotifs.pl fasta_file1.fa human analysis_output/ > output.txt
# 
# Selected Options:
#   Input file = fasta_file1.fa
# Promoter Set = human
# Output Directory = analysis_output/
#   Found mset for "human", will check against vertebrates motifs
# 
# Progress: Step1 - Convert input file to refseq IDs
# Percentage of IDs converted into refseq: 0.0% (0 out of 12)
# !!!! Homer converted less than 5% of the rows in your file.
# !!!! Check to be sure the input file has valid gene identifiers
# !!!! Check to be sure the input file has valid gene identifiers
# kat-juno:homer kathryntsai$ 

####

# d[1][["TRAF2"]]
# [1] " " " " "g" "t" "g" "a" "c" "t" "g" "g" "t" "c" "c" "t" "g" "g" "c" "g" "g" "g" "t" "g" "c" "t" "t"
# [26] "a" "g" "g" "a" "g" "g" "c" "g" "g" "g" "g" "g" "c" "a" "g" "a" "g" "a" "c" "a" "t" "g" "g" "c" "t"
# [51] "c" "a" "c" "t" "c" "c" "a" "g" "c" "c" "t" "c" " " " " "a" "a" "c" "c" "t" "c" "c" "t" "g" "g" "g"
# [76] "c" "t" "c" "a" "a" "g" "c" "a" "a" "t" "t" "c" "t" "c" "c" "c" "a" "c" "c" "t" "c" "a" "g" "c" "c"
# [101] "a" "c" "t" "c" "a" "g" "g" "t" "a" "g" "c" "t" "g" "g" "g" "g" "c" "t" "a" "c" "a" "g" "g" "c" " "
# [126] " " "g" "t" "g" "c" "a" "c" "c" "a" "t" "c" "a" "c" "g" "c" "c" "c" "g" "g" "c" "t" "a" "t" "t" "t"
# [151] "t" "t" "t" "t" "t" "t" "t" "t" "t" "t" "t" "t" "t" "g" "a" "g" "a" "c" "a" "g" "a" "a" "t" "t" "t"
# [176] "g" "g" "c" "t" "c" "t" "t" "a" "t" "t" "g" " " " " "c" "c" "c" "a" "g" "g" "c" "t" "g" "g" "a" "g"
# [201] "t" "g" "c" "a" "a" "t" "g" "g" "t" "g" "c" "a" "a" "t" "c" "t" "c" "a" "g" "c" "c" "a" "c" "a" "a"
# [226] "c" "c" "t" "c" "t" "g" "c" "c" "t" "c" "c" "c" "a" "g" "g" "g" "t" "c" "a" "a" "g" "t" "g" " " " "
# [251] "a" "t" "t" "c" "t" "c" "c" "t" "g" "c" "c" "t" "c" "a" "g" "c" "c" "t" "c" "c" "c" "a" "a" "g" "t"
# [276] "a" "g" "c" "t" "a" "g" "g" "a" "t" "t" "a" "c" "a" "g" "g" "t" "g" "t" "g" "c" "a" "c" "c" "a" "c"
# [301] "c" "a" "c" "g" "c" "c" "g" "g" "g" "c" " " " " "t" "a" "a" "t" "t" "t" "t" "t" "g" "t" "a" "t" "t"
# [326] "t" "t" "t" "a" "g" "t" "a" "g" "a" "a" "a" "c" "g" "g" "g" "g" "t" "t" "t" "c" "t" "c" "c" "a" "t"
# [351] "g" "t" "t" "g" "g" "t" "c" "a" "g" "g" "t" "t" "g" "g" "t" "c" "t" "c" "g" "a" "a" "c" " " " " "t"
# [376] "c" "c" "c" "a" "a" "c" "c" "t" "c" "a" "g" "g" "t" "g" "a" "t" "c" "t" "g" "c" "c" "t" "g" "c" "c"
# [401] "t" "c" "a" "g" "c" "c" "t" "c" "c" "c" "a" "a" "a" "g" "t" "g" "c" "t" "g" "g" "g" "a" "t" "t" "a"
# [426] "c" "a" "g" "g" "c" "g" "t" "g" "a" " " " " "g" "c" "c" "a" "c" "c" "a" "c" "a" "c" "c" "c" "g" "g"
# [451] "c" "c" "t" "a" "a" "t" "t" "t" "t" "t" "g" "c" "a" "t" "t" "t" "t" "t" "t" "g" "g" "a" "g" "a" "g"
# [476] "a" "c" "g" "g" "a" "g" "g" "t" "t" "t" "c" "a" "c" "t" "a" "t" "a" "t" "t" "g" "c" " " " " "c" "c"
# [501] "a" "g" "g" "c" "t" "g" "g" "t" "c" "t" "t" "g" "a" "a" "c" "t" "c" "c" "c" "g" "a" "g" "c" "t" "c"
# [526] "a" "a" "g" "t" "g" "a" "t" "c" "c" "a" "g" "c" "c" "t" "t" "g" "g" "c" "c" "t" "c" "c" "c" "a" "a"
# [551] "a" "g" "t" "g" "c" "t" "g" "g" " " " " "g" "a" "t" "t" "a" "c" "a" "g" "g" "t" "g" "t" "g" "a" "a"
# [576] "t" "c" "t" "t" "c" "t" "a" "t" "t" "a" "c" "t" "t" "t" "c" "c" "t" "g" "t" "g" "t" "g" "t" "c" "c"
# [601] "t" "t" "c" "c" "c" "a" "a" "g" "t" "t" "c" "a" "c" "a" "g" "a" "a" "a" "a" "c" " " " " "a" "c" "a"
# [626] "a" "g" "c" "t" "g" "a" "c" "a" "g" "a" "a" "a" "t" "c" "c" "a" "t" "t" "c" "t" "t" "t" "c" "c" "c"
# [651] "a" "c" "g" "g" "t" "t" "g" "c" "t" "g" "g" "c" "a" "c" "g" "g" "a" "t" "g" "c" "t" "g" "c" "t" "g"
# [676] "t" "g" "t" "g" "a" "c" "a" " " " " "g" "c" "a" "c" "a" "a" "g" "t" "t" "c" "a" "t" "g" "t" "g" "c"
# [701] "t" "t" "g" "a" "g" "g" "g" "t" "c" "t" "c" "c" "t" "g" "g" "c" "c" "a" "c" "a" "g" "c" "c" "a" "g"
# [726] "c" "c" "c" "c" "a" "g" "g" "c" "t" "t" "g" "c" "c" "c" "t" "g" "g" "a" "c" " " " " "a" "g" "g" "c"
# [751] "g" "g" "g" "c" "t" "t" "c" "a" "g" "a" "t" "t" "c" "c" "c" "c" "c" "c" "a" "a" "g" "g" "t" "t" "c"
# [776] "t" "a" "t" "g" "c" "c" "a" "t" "t" "g" "a" "g" "t" "t" "g" "g" "a" "g" "g" "a" "c" "a" "a" "g" "a"
# [801] "g" "c" "t" "t" "c" "c" " " " " "c" "c" "a" "c" "t" "g" "t" "c" "c" "a" "g" "c" "g" "t" "t" "g" "a"
# [826] "a" "g" "c" "a" "a" "a" "c" "c" "t" "c" "c" "c" "t" "c" "c" "a" "g" "g" "t" "c" "c" "c" "a" "c" "c"
# [851] "c" "c" "t" "c" "c" "c" "c" "a" "g" "g" "c" "c" "t" "t" "t" "g" "t" "c" " " " " "c" "t" "t" "t" "t"
# [876] "a" "t" "c" "a" "g" "g" "g" "a" "t" "t" "g" "c" "c" "c" "t" "g" "g" "g" "g" "g" "a" "g" "c" "c" "c"
# [901] "c" "c" "c" "t" "g" "g" "g" "g" "a" "g" "c" "t" "g" "a" "t" "g" "a" "g" "g" "a" "g" "g" "g" "a" "g"
# [926] "c" "c" "t" "c" "t" " " " " "g" "t" "g" "c" "c" "c" "c" "g" "t" "g" "c" "a" "c" "t" "g" "g" "c" "c"
# [951] "t" "g" "a" "g" "g" "g" "g" "g" "t" "g" "a" "g" "g" "c" "c" "a" "g" "g" "g" "c" "t" "g" "g" "g" "g"
# [976] "c" "t" "g" "c" "t" "g" "c" "a" "g" "g" "g" "g" "a" "c" "t" "c" "t" " " " " "g" "g" "c" "c" "c" "g"
# [ reached getOption("max.print") -- omitted 34 entries ]
# attr(,"name")
# [1] "TRAF2"
# attr(,"Annot")
# [1] ">TRAF2"
# attr(,"class")
# [1] "SeqFastadna"
