# UNFINISHED: Attempt Bioconductor packages to retrieve fasta promoter sequences, along with 
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

a_gene_snp <- franke_cis_data[which(Gene==c('ENSG00000127191', 
                              'ENSG00000214279', 
                              'ENSG00000168765', 
                              'ENSG00000125966',
                              'ENSG00000258531',
                              'ENSG00000237437'
              )), c('Gene', 'GeneSymbol', 'GeneChr', 'GenePos', 'SNP', 'SNPChr', 'SNPPos')]
a_gene <- a_gene_snp[!duplicated(a_gene_snp$Gene),]


# http://grch37.ensembl.org/Homo_sapiens/Gene/Sequence?db=core;g=ENSG00000127191
# http://grch37.ensembl.org/Homo_sapiens/Gene/Sequence?db=core;g=ENSG00000214279
# http://grch37.ensembl.org/Homo_sapiens/Gene/Sequence?db=core;g=ENSG00000168765
# http://grch37.ensembl.org/Homo_sapiens/Gene/Sequence?db=core;g=ENSG00000125966
# http://grch37.ensembl.org/Homo_sapiens/Gene/Sequence?db=core;g=ENSG00000258531
# http://grch37.ensembl.org/Homo_sapiens/Gene/Sequence?db=core;g=ENSG00000237437


aa <- data.frame(a, rbind(
  '>chromosome:GRCh37:9:139775364:139821659:1
GTGACTGGTCCTGGCGGGTGCTTAGGAGGCGGGGGCAGAGACATGGCTCACTCCAGCCTCAACCTCCTGGGCTCAAGCAATTCTCCCACCTCAGCCACTCAGGTAGCTGGGGCTACAGGCGTGCACCATCACGCCCGGCTATTTTTTTTTTTTTTTTGAGACAGAATTTGGCTCTTATTGCCCAGGCTGGAGTGCAATGGTGCAATCTCAGCCACAACCTCTGCCTCCCAGGGTCAAGTGATTCTCCTGCCTCAGCCTCCCAAGTAGCTAGGATTACAGGTGTGCACCACCACGCCGGGCTAATTTTTGTATTTTTAGTAGAAACGGGGTTTCTCCATGTTGGTCAGGTTGGTCTCGAACTCCCAACCTCAGGTGATCTGCCTGCCTCAGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCACACCCGGCCTAATTTTTGCATTTTTTGGAGAGACGGAGGTTTCACTATATTGCCCAGGCTGGTCTTGAACTCCCGAGCTCAAGTGATCCAGCCTTGGCCTCCCAAAGTGCTGGGATTACAGGTGTGAATCTTCTATTACTTTCCTGTGTGTCCTTCCCAAGTTCACAGAAAACACAAGCTGACAGAAATCCATTCTTTCCCACGGTTGCTGGCACGGATGCTGCTGTGTGACAGCACAAGTTCATGTGCTTGAGGGTCTCCTGGCCACAGCCAGCCCCAGGCTTGCCCTGGACAGGCGGGCTTCAGATTCCCCCCAAGGTTCTATGCCATTGAGTTGGAGGACAAGAGCTTCCCCACTGTCCAGCGTTGAAGCAAACCTCCCTCCAGGTCCCACCCCTCCCCAGGCCTTTGTCCTTTTATCAGGGATTGCCCTGGGGGAGCCCCCCTGGGGAGCTGATGAGGAGGGAGCCTCTGTGCCCCGTGCACTGGCCTGAGGGGGTGAGGCCAGGGCTGGGGCTGCTGCAGGGGACTCTGGCCCGGAGCTGCGGCTGCCTTCCTGGAAACACTCCTGTC',

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

# https://www.ncbi.nlm.nih.gov/snp/rs1033797
a_gene <- data.frame(a_gene, rbind(
  '>chromosome:GRCh37:9:139775364:139821659:1 #139798711
  GTGACTGGTCCTGGCGGGTGCTTAGGAGGCGGGGGCAGAGACATGGCTCACTCCAGCCTCAACCTCCTGGGCTCAAGCAATTCTCCCACCTCAGCCACTCAGGTAGCTGGGGCTACAGGCGTGCACCATCACGCCCGGCTATTTTTTTTTTTTTTTTGAGACAGAATTTGGCTCTTATTGCCCAGGCTGGAGTGCAATGGTGCAATCTCAGCCACAACCTCTGCCTCCCAGGGTCAAGTGATTCTCCTGCCTCAGCCTCCCAAGTAGCTAGGATTACAGGTGTGCACCACCACGCCGGGCTAATTTTTGTATTTTTAGTAGAAACGGGGTTTCTCCATGTTGGTCAGGTTGGTCTCGAACTCCCAACCTCAGGTGATCTGCCTGCCTCAGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCACACCCGGCCTAATTTTTGCATTTTTTGGAGAGACGGAGGTTTCACTATATTGCCCAGGCTGGTCTTGAACTCCCGAGCTCAAGTGATCCAGCCTTGGCCTCCCAAAGTGCTGGGATTACAGGTGTGAATCTTCTATTACTTTCCTGTGTGTCCTTCCCAAGTTCACAGAAAACACAAGCTGACAGAAATCCATTCTTTCCCACGGTTGCTGGCACGGATGCTGCTGTGTGACAGCACAAGTTCATGTGCTTGAGGGTCTCCTGGCCACAGCCAGCCCCAGGCTTGCCCTGGACAGGCGGGCTTCAGATTCCCCCCAAGGTTCTATGCCATTGAGTTGGAGGACAAGAGCTTCCCCACTGTCCAGCGTTGAAGCAAACCTCCCTCCAGGTCCCACCCCTCCCCAGGCCTTTGTCCTTTTATCAGGGATTGCCCTGGGGGAGCCCCCCTGGGGAGCTGATGAGGAGGGAGCCTCTGTGCCCCGTGCACTGGCCTGAGGGGGTGAGGCCAGGGCTGGGGCTGCTGCAGGGGACTCTGGCCCGGAGCTGCGGCTGCCTTCCTGGAAACACTCCTGTC',
  
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

colnames(a_gene)[8] <- "Fasta Sequence Changed"

a_gene$SNP_Subtraction <- a_gene$GenePos - a_gene$SNPPos

write.table(aa$"Fasta Sequence Changed", 'q2_3_output/fasta_file1.fa', sep='\t', quote=F, row.names = F, col.names = F)

# Try with bioconductor

##  kathryntsai$ findMotifs.pl fasta_file1.fa fasta analysis_output/ > output.txt --> couldn't open revised fasta file?

## kathryntsai$ findMotifs.pl /Users/kathryntsai/OneDrive\ -\ Villanova\ University/College/2018-2019/Summer\ 2019/TFs_eQTLs_Research/GitHub/eQTL_Kathryn/kathryn/bulk_data/q2_3_output/fasta_file1_old2.fa fasta analysis_output/ > output.txt

# Selected Options:
#   Input file = /Users/kathryntsai/OneDrive - Villanova University/College/2018-2019/Summer 2019/TFs_eQTLs_Research/GitHub/eQTL_Kathryn/kathryn/bulk_data/q2_3_output/fasta_file1_old2.fa
# Promoter Set = fasta
# Output Directory = analysis_output/
#   
#   !Warning - no background FASTA file specified (Highly recommended)
# !Your input sequences will be randomized to serve as a background instead.
# 
# Found 6 sequences
# Using custom gene IDs for GO analysis
# Parsing FASTA format files...
# Found 6 sequences
# Found 30 sequences
# 
# Progress: Step4 - removing redundant promoters
# 
# Progress: Step5 - adjusting background sequences for GC/CpG content...
# 
# Sequences processed:
#   Auto detected maximum sequence length of 1000 bp
# 36 total
# 
# Frequency Bins: 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.6 0.7 0.8
# Freq	Bin	Count
# 0.4	4	5
# 0.45	5	1
# 0.5	6	6
# 0.6	7	21
# 0.7	8	3
# Bin	# Targets	# Background	Background Weight
# 4	1	4	1.083
# 6	1	5	0.867
# 7	4	17	1.020
# 
# Normalizing lower order oligos using homer2
# 
# Reading input files...
# 32 total sequences read
# Autonormalization: 1-mers (4 total)
# A	24.47%	24.67%	0.992
# C	25.53%	25.33%	1.008
# G	25.53%	25.33%	1.008
# T	24.47%	24.67%	0.992
# Autonormalization: 2-mers (16 total)
# AA	7.19%	6.36%	1.130
# CA	7.25%	6.11%	1.188
# GA	6.25%	6.16%	1.014
# TA	3.80%	6.04%	0.629
# AC	4.70%	6.09%	0.771
# CC	8.18%	6.65%	1.230
# GC	6.42%	6.43%	0.999
# TC	6.25%	6.16%	1.014
# AG	7.86%	6.21%	1.265
# CG	2.23%	6.38%	0.349
# GG	8.18%	6.65%	1.230
# TG	7.25%	6.11%	1.188
# AT	4.72%	6.00%	0.787
# CT	7.86%	6.21%	1.265
# GT	4.70%	6.09%	0.771
# TT	7.19%	6.36%	1.130
# Autonormalization: 3-mers (64 total)
# Normalization weights can be found in file: analysis_output//seq.autonorm.tsv
# Converging on autonormalization solution:
#   ...............................................................................
# Final normalization:	Autonormalization: 1-mers (4 total)
# A	24.47%	24.41%	1.003
# C	25.53%	25.59%	0.998
# G	25.53%	25.59%	0.998
# T	24.47%	24.41%	1.003
# Autonormalization: 2-mers (16 total)
# AA	7.19%	6.28%	1.145
# CA	7.25%	6.16%	1.179
# GA	6.25%	6.18%	1.010
# TA	3.80%	5.80%	0.655
# AC	4.70%	6.05%	0.777
# CC	8.18%	6.82%	1.199
# GC	6.42%	6.53%	0.983
# TC	6.25%	6.18%	1.010
# AG	7.86%	6.27%	1.253
# CG	2.23%	6.35%	0.350
# GG	8.18%	6.82%	1.199
# TG	7.25%	6.16%	1.179
# AT	4.72%	5.81%	0.812
# CT	7.86%	6.27%	1.253
# GT	4.70%	6.05%	0.777
# TT	7.19%	6.28%	1.145
# Autonormalization: 3-mers (64 total)
# 
# Progress: Step6 - Gene Ontology Enrichment Analysis
# Skipping...
# 
# Progress: Step7 - Known motif enrichment
# 
# Reading input files...
# 32 total sequences read
# 994 motifs loaded
# Cache length = 11180
# Using hypergeometric scoring
# Checking enrichment of 994 motif(s)
# |0%                                    50%                                  100%|
#   =================================================================================
#   Illegal division by zero at /Users/kathryntsai/homer/bin/findKnownMotifs.pl line 153.
# 
# Progress: Step8 - De novo motif finding (HOMER)
# 
# Scanning input files...
# Parsing sequences...
# |0%                                   50%                                  100%|
#   ================================
#   Total number of Oligos: 19956
# Autoadjustment for sequence coverage in background: 1.00x
# 
# Oligos: 19956 of 31856 max
# Tree  : 44444 of 159280 max
# Optimizing memory usage...
# Cache length = 11180
# Using hypergeometric scoring
# 
# Global Optimization Phase: Looking for enriched oligos with up to 1 mismatches...
# 
# Screening oligos 19956 (allowing 0 mismatches):
#   |0%                                   50%                                  100%|
#   ================================================================================
#   76.18% skipped, 23.82% checked (4754 of 19956), of those checked:
#   76.18% not in target, 0.00% increased p-value, 0.00% high p-value
# 
# Screening oligos 19956 (allowing 1 mismatches):
#   |0%                                   50%                                  100%|
#   ================================================================================
#   76.18% skipped, 23.82% checked (4754 of 19956), of those checked:
#   0.00% not in target, 18.27% increased p-value, 6.95% high p-value
# Reading input files...
# 32 total sequences read
# Cache length = 11180
# Using hypergeometric scoring
# 
# Local Optimization Phase:
#   1 of 25 Initial Sequence: TCCAGCCT... (-13.313)
# Round 1: -17.40 TCCAGCCT T:13.0(88.78%),B:0.0(0.00%),P:1e-7
# Round 2: -17.40 TCCAGCCT T:13.0(88.78%),B:0.0(0.00%),P:1e-7
# =Final=: -22.17 TCCAGCCT T:6.0(100.00%),B:0.0(0.00%),P:1e-9
# Performing exhaustive masking of motif...
# Reprioritizing potential motifs...
# 2 of 25 Initial Sequence: TTTTTTTT... (-10.297)
# Round 1: -10.30 TTTTTTTT T:34.0(99.76%),B:9.4(7.40%),P:1e-4
# Round 2: -10.30 TTTTTTTT T:34.0(99.76%),B:9.4(7.40%),P:1e-4
# =Final=: -9.86 TTTTTTTT T:5.0(83.33%),B:9.4(7.99%),P:1e-4
# Performing exhaustive masking of motif...
# Reprioritizing potential motifs...
# 3 of 25 Initial Sequence: GGATTACA... (-9.799)
# Round 1: -11.31 GGATTACA T:10.0(83.85%),B:7.7(5.81%),P:1e-4
# Round 2: -11.31 GGATTACA T:10.0(83.85%),B:7.7(5.81%),P:1e-4
# =Final=: -12.59 GGATTACA T:5.0(83.33%),B:4.2(3.60%),P:1e-5
# Performing exhaustive masking of motif...
# Reprioritizing potential motifs...
# 4 of 25 Initial Sequence: CAGCCTCC... (-9.625)
# Round 1: -17.40 CAGCCTCC T:16.0(93.51%),B:0.0(0.00%),P:1e-7
# Round 2: -17.40 CAGCCTCC T:16.0(93.51%),B:0.0(0.00%),P:1e-7
# =Final=: -22.17 CAGCCTCC T:6.0(100.00%),B:0.0(0.00%),P:1e-9
# Performing exhaustive masking of motif...
# Reprioritizing potential motifs...
# 5 of 25 Initial Sequence: TCCTGCCT... (-9.625)
# Round 1: -9.62 TCCTGCCT T:5.0(59.81%),B:0.0(0.00%),P:1e-4
# Round 2: -9.62 TCCTGCCT T:5.0(59.81%),B:0.0(0.00%),P:1e-4
# =Final=: -9.62 TCCTGCCT T:3.0(50.00%),B:0.0(0.00%),P:1e-4
# Performing exhaustive masking of motif...
# Reprioritizing potential motifs...
# 6 of 25 Initial Sequence: AACTCCTG... (-9.625)
# Round 1: -9.62 AACTCCTG T:4.0(51.77%),B:0.0(0.00%),P:1e-4
# Round 2: -9.62 AACTCCTG T:4.0(51.77%),B:0.0(0.00%),P:1e-4
# =Final=: -13.31 AACTCCTG T:4.0(66.67%),B:0.0(0.00%),P:1e-5
# Performing exhaustive masking of motif...
# Reprioritizing potential motifs...
# 7 of 25 Initial Sequence: TGAGCCAC... (-9.625)
# Round 1: -15.61 TGAGCCAC T:10.0(83.85%),B:2.6(1.69%),P:1e-6
# Round 2: -15.61 TGAGCCAC T:10.0(83.85%),B:2.6(1.69%),P:1e-6
# =Final=: -17.40 TGAGCCAC T:5.0(83.33%),B:0.0(0.00%),P:1e-7
# Performing exhaustive masking of motif...
# Reprioritizing potential motifs...
# 8 of 25 Initial Sequence: AAAGTGCT... (-9.625)
# Round 1: -11.72 AAAGTGCT T:7.0(72.09%),B:2.6(1.69%),P:1e-5
# Round 2: -13.40 CAARTGYT T:10.0(83.85%),B:5.0(3.36%),P:1e-5
# Round 3: -13.40 CAARTGYT T:10.0(83.85%),B:5.0(3.36%),P:1e-5
# =Final=: -9.12 CAARTGYT T:4.0(66.67%),B:5.0(4.27%),P:1e-3
# Performing exhaustive masking of motif...
# Reprioritizing potential motifs...
# Remaining seeds don't look promising (After initial 5 motifs, logp -6.215 > -6.778)
# 
# Finalizing Enrichment Statistics (new in v3.4)
# Reading input files...
# 32 total sequences read
# Cache length = 11180
# Using hypergeometric scoring
# Checking enrichment of 8 motif(s)
# |0%                                    50%                                  100%|
# =================================================================================
# Output in file: analysis_output//homerMotifs.motifs8
# 
# 
# Scanning input files...
# Parsing sequences...
# |0%                                   50%                                  100%|
# ================================
# Total number of Oligos: 30482
# Autoadjustment for sequence coverage in background: 1.00x
# 
# Oligos: 30482 of 31792 max
# Tree  : 83388 of 158960 max
# Optimizing memory usage...
# Cache length = 11180
# Using hypergeometric scoring
# 
# Global Optimization Phase: Looking for enriched oligos with up to 1 mismatches...
# 
# Screening oligos 30482 (allowing 0 mismatches):
# |0%                                   50%                                  100%|
# ================================================================================
# 81.70% skipped, 18.30% checked (5577 of 30482), of those checked:
# 81.70% not in target, 0.00% increased p-value, 0.00% high p-value
# 
# Screening oligos 30482 (allowing 1 mismatches):
# |0%                                   50%                                  100%|
# ================================================================================
# 81.70% skipped, 18.30% checked (5577 of 30482), of those checked:
# 0.00% not in target, 12.01% increased p-value, 2.02% high p-value
# Reading input files...
# 32 total sequences read
# Cache length = 11180
# Using hypergeometric scoring
# 
# Local Optimization Phase:
# 1 of 25 Initial Sequence: CTCCTGCCTC... (-17.399)
# Round 1: -17.40 CTCCTGCCTC T:10.0(83.85%),B:0.0(0.00%),P:1e-7
# Round 2: -17.40 CTCCTGCCTC T:10.0(83.85%),B:0.0(0.00%),P:1e-7
# =Final=: -13.31 CTCCTGCCTC T:4.0(66.67%),B:0.0(0.00%),P:1e-5
# Performing exhaustive masking of motif...
# Reprioritizing potential motifs...
# 2 of 25 Initial Sequence: TTTTTTTTTT... (-17.399)
# Round 1: -17.40 TTTTTTTTTT T:25.0(98.74%),B:0.0(0.00%),P:1e-7
# Round 2: -17.40 TTTTTTTTTT T:25.0(98.74%),B:0.0(0.00%),P:1e-7
# =Final=: -13.31 TTTTTTTTTT T:4.0(66.67%),B:0.0(0.00%),P:1e-5
# Performing exhaustive masking of motif...
# Reprioritizing potential motifs...
# 3 of 25 Initial Sequence: ACTCCAGCCT... (-13.313)
# Round 1: -17.40 ACTCCAGCCT T:11.0(83.85%),B:0.0(0.00%),P:1e-7
# Round 2: -17.40 ACTCCAGCCT T:11.0(83.85%),B:0.0(0.00%),P:1e-7
# =Final=: -22.17 ACTCCAGCCT T:6.0(100.00%),B:0.0(0.00%),P:1e-9
# Performing exhaustive masking of motif...
# Reprioritizing potential motifs...
# 4 of 25 Initial Sequence: CCTGGGCTCA... (-13.313)
# Round 1: -17.40 CCTGGGCTCA T:10.0(83.85%),B:0.0(0.00%),P:1e-7
# Round 2: -17.40 CCTGGGCTCA T:10.0(83.85%),B:0.0(0.00%),P:1e-7
# =Final=: -17.40 CCTGGGCTCA T:5.0(83.33%),B:0.0(0.00%),P:1e-7
# Performing exhaustive masking of motif...
# Reprioritizing potential motifs...
# 5 of 25 Initial Sequence: ATTACAGGCG... (-13.313)
# Round 1: -13.31 ATTACAGGCG T:7.0(72.09%),B:0.0(0.00%),P:1e-5
# Round 2: -13.31 ATTACAGGCG T:7.0(72.09%),B:0.0(0.00%),P:1e-5
# =Final=: -13.31 ATTACAGGCG T:4.0(66.67%),B:0.0(0.00%),P:1e-5
# Performing exhaustive masking of motif...
# Reprioritizing potential motifs...
# 6 of 25 Initial Sequence: GCCTCCCAAG... (-13.313)
# Round 1: -17.40 SCYTCCCAAG T:11.0(83.85%),B:0.0(0.00%),P:1e-7
# Round 2: -17.40 SCYTCCCAAG T:11.0(83.85%),B:0.0(0.00%),P:1e-7
# =Final=: -22.17 SCYTCCCAAG T:6.0(100.00%),B:0.0(0.00%),P:1e-9
# Performing exhaustive masking of motif...
# Reprioritizing potential motifs...
# 7 of 25 Initial Sequence: GACGGAGGTT... (-9.625)
# Round 1: -9.62 GACGGAGGTT T:4.0(51.77%),B:0.0(0.00%),P:1e-4
# Round 2: -9.62 GACGGAGGTT T:4.0(51.77%),B:0.0(0.00%),P:1e-4
# =Final=: -9.62 GACGGAGGTT T:3.0(50.00%),B:0.0(0.00%),P:1e-4
# Performing exhaustive masking of motif...
# Reprioritizing potential motifs...
# 8 of 25 Initial Sequence: CCACCACGCC... (-9.625)
# Round 1: -9.62 CCACCACGCC T:6.0(59.81%),B:0.0(0.00%),P:1e-4
# Round 2: -9.62 CCACCACGCC T:6.0(59.81%),B:0.0(0.00%),P:1e-4
# =Final=: -13.31 CCACCACGCC T:4.0(66.67%),B:0.0(0.00%),P:1e-5
# Performing exhaustive masking of motif...
# Reprioritizing potential motifs...
# 9 of 25 Initial Sequence: TTTTTCCATT... (-9.625)
# Round 1: -13.31 TTTTTCCATT T:7.0(72.09%),B:0.0(0.00%),P:1e-5
# Round 2: -17.40 TTTTTYTATT T:10.0(83.85%),B:0.0(0.00%),P:1e-7
# Round 3: -17.40 TTTTTYTATT T:10.0(83.85%),B:0.0(0.00%),P:1e-7
# =Final=: -17.74 TTTTTYTATT T:6.0(100.00%),B:3.4(2.92%),P:1e-7
# Performing exhaustive masking of motif...
# Reprioritizing potential motifs...
# 10 of 25 Initial Sequence: CCCAGAGCCC... (-9.625)
# Round 1: -9.62 CCCAGAGCCC T:5.0(51.77%),B:0.0(0.00%),P:1e-4
# Round 2: -9.62 CCCAGAGCCC T:5.0(51.77%),B:0.0(0.00%),P:1e-4
# =Final=: -9.62 CCCAGAGCCC T:3.0(50.00%),B:0.0(0.00%),P:1e-4
# Performing exhaustive masking of motif...
# Reprioritizing potential motifs...
# 11 of 25 Initial Sequence: GTCTTGAACT... (-9.625)
# Round 1: -9.62 GTCTTGAACT T:4.0(51.77%),B:0.0(0.00%),P:1e-4
# Round 2: -9.62 GTCTTGAACT T:4.0(51.77%),B:0.0(0.00%),P:1e-4
# =Final=: -9.62 GTCTTGAACT T:3.0(50.00%),B:0.0(0.00%),P:1e-4
# Performing exhaustive masking of motif...
# Reprioritizing potential motifs...
# 12 of 25 Initial Sequence: GACACAACCA... (-9.625)
# Round 1: -13.31 GMCACAACCA T:7.0(72.09%),B:0.0(0.00%),P:1e-5
# Round 2: -13.31 GMCACAACCA T:7.0(72.09%),B:0.0(0.00%),P:1e-5
# =Final=: -13.31 GMCACAACCA T:4.0(66.67%),B:0.0(0.00%),P:1e-5
# Performing exhaustive masking of motif...
# Reprioritizing potential motifs...
# 13 of 25 Initial Sequence: TCCCAGGGAC... (-9.625)
# Round 1: -17.40 TCCCAGGGAC T:11.0(83.85%),B:0.0(0.00%),P:1e-7
# Round 2: -17.40 TCCCAGGGAC T:11.0(83.85%),B:0.0(0.00%),P:1e-7
# =Final=: -22.17 TCCCAGGGAC T:6.0(100.00%),B:0.0(0.00%),P:1e-9
# Performing exhaustive masking of motif...
# Reprioritizing potential motifs...
# 14 of 25 Initial Sequence: TTAAAAAATC... (-9.625)
# Round 1: -13.31 TTAAAAAATC T:7.0(72.09%),B:0.0(0.00%),P:1e-5
# Round 2: -13.31 TTAAAAAATC T:7.0(72.09%),B:0.0(0.00%),P:1e-5
# =Final=: -9.62 TTAAAAAATC T:3.0(50.00%),B:0.0(0.00%),P:1e-4
# Performing exhaustive masking of motif...
# Reprioritizing potential motifs...
# 15 of 25 Initial Sequence: GGCAGCCTCA... (-9.625)
# Round 1: -9.62 GGCAGCCTCA T:4.0(51.77%),B:0.0(0.00%),P:1e-4
# Round 2: -9.62 GGCAGCCTCA T:4.0(51.77%),B:0.0(0.00%),P:1e-4
# =Final=: -13.31 GGCAGCCTCA T:4.0(66.67%),B:0.0(0.00%),P:1e-5
# Performing exhaustive masking of motif...
# Reprioritizing potential motifs...
# 16 of 25 Initial Sequence: CTGGTGCTGC... (-8.257)
# Round 1: -11.72 CTGGTGCTGY T:8.0(72.09%),B:2.6(1.69%),P:1e-5
# Round 2: -13.31 CTGGTSCTGC T:8.0(72.09%),B:0.0(0.00%),P:1e-5
# Round 3: -13.31 CTGGTSCTGC T:8.0(72.09%),B:0.0(0.00%),P:1e-5
# =Final=: -13.31 CTGGTSCTGC T:4.0(66.67%),B:0.0(0.00%),P:1e-5
# Performing exhaustive masking of motif...
# Reprioritizing potential motifs...
# 17 of 25 Initial Sequence: ACAGCCCAGG... (-8.257)
# Round 1: -10.78 ACAGCCCAGG T:10.0(83.85%),B:8.6(6.61%),P:1e-4
# Round 2: -13.31 ACAGAMCAGG T:7.0(72.09%),B:0.0(0.00%),P:1e-5
# Round 3: -14.37 WCAGMVCAGG T:10.0(83.85%),B:3.1(2.53%),P:1e-6
# Round 4: -14.37 WCAGMVCAGG T:10.0(83.85%),B:3.1(2.53%),P:1e-6
# =Final=: -17.40 WCAGMVCAGG T:5.0(83.33%),B:0.0(0.00%),P:1e-7
# Performing exhaustive masking of motif...
# Reprioritizing potential motifs...
# 18 of 25 Initial Sequence: GCAACAGAGC... (-7.360)
# Round 1: -9.62 GCAACAGAGC T:4.0(51.77%),B:0.0(0.00%),P:1e-4
# Round 2: -11.72 GCAATATAGT T:7.0(72.09%),B:2.4(1.69%),P:1e-5
# Round 3: -11.72 GCAATATAGT T:7.0(72.09%),B:2.4(1.69%),P:1e-5
# =Final=: -10.63 GCAATATAGT T:4.0(66.67%),B:2.4(2.04%),P:1e-4
# Performing exhaustive masking of motif...
# Reprioritizing potential motifs...
# Remaining seeds don't look promising (After initial 5 motifs, logp -6.686 > -6.778)
# 
# Finalizing Enrichment Statistics (new in v3.4)
# Reading input files...
# 32 total sequences read
# Cache length = 11180
# Using hypergeometric scoring
# Checking enrichment of 18 motif(s)
# |0%                                    50%                                  100%|
#   =================================================================================
#   Output in file: analysis_output//homerMotifs.motifs10
# 
# 
# -blen automatically set to 2
# Scanning input files...
# Parsing sequences...
# |0%                                   50%                                  100%|
#   ================================
#   Total number of Oligos: 31432
# Autoadjustment for sequence coverage in background: 1.00x
# 
# Oligos: 31432 of 31728 max
# Tree  : 90540 of 158640 max
# Optimizing memory usage...
# Cache length = 11180
# Using hypergeometric scoring
# 
# Global Optimization Phase: Looking for enriched oligos with up to 1 mismatches...
# 
# Screening oligos 31432 (allowing 0 mismatches):
#   |0%                                   50%                                  100%|
#   ================================================================================
#   81.78% skipped, 18.22% checked (5727 of 31432), of those checked:
#   81.78% not in target, 0.00% increased p-value, 0.00% high p-value
# 
# Screening oligos 31432 (allowing 1 mismatches):
#   |0%                                   50%                                  100%|
#   ================================================================================
#   81.78% skipped, 18.22% checked (5727 of 31432), of those checked:
#   0.00% not in target, 3.13% increased p-value, 0.36% high p-value
# Reading input files...
# 32 total sequences read
# Cache length = 11180
# Using hypergeometric scoring
# 
# Local Optimization Phase:
#   1 of 25 Initial Sequence: TATTTTTTTTTT... (-17.399)
# Round 1: -17.40 TATTTTTTTTTT T:20.0(96.87%),B:0.0(0.00%),P:1e-7
# Round 2: -17.40 TATTTTTTTTTT T:20.0(96.87%),B:0.0(0.00%),P:1e-7
# =Final=: -13.31 TATTTTTTTTTT T:4.0(66.67%),B:0.0(0.00%),P:1e-5
# Performing exhaustive masking of motif...
# Reprioritizing potential motifs...
# 2 of 25 Initial Sequence: GTAATCCCAGCT... (-13.313)
# Round 1: -17.40 GYAATCCCAGCT T:10.0(83.85%),B:0.0(0.00%),P:1e-7
# Round 2: -17.40 GYAATCCCAGCT T:10.0(83.85%),B:0.0(0.00%),P:1e-7
# =Final=: -17.40 GYAATCCCAGCT T:5.0(83.33%),B:0.0(0.00%),P:1e-7
# Performing exhaustive masking of motif...
# Reprioritizing potential motifs...
# 3 of 25 Initial Sequence: ACCTCAGCCTCC... (-13.313)
# Round 1: -17.40 ACCTCAGCCTCC T:13.0(88.78%),B:0.0(0.00%),P:1e-7
# Round 2: -17.40 ACCTCAGCCTCC T:13.0(88.78%),B:0.0(0.00%),P:1e-7
# =Final=: -13.31 ACCTCAGCCTCC T:4.0(66.67%),B:0.0(0.00%),P:1e-5
# Performing exhaustive masking of motif...
# Reprioritizing potential motifs...
# 4 of 25 Initial Sequence: GAGGCCAGGAGT... (-9.625)
# Round 1: -17.40 GAGGCCAGGAGT T:10.0(83.85%),B:0.0(0.00%),P:1e-7
# Round 2: -17.40 GAGGCCAGGAGT T:10.0(83.85%),B:0.0(0.00%),P:1e-7
# =Final=: -13.31 GAGGCCAGGAGT T:4.0(66.67%),B:0.0(0.00%),P:1e-5
# Performing exhaustive masking of motif...
# Reprioritizing potential motifs...
# 5 of 25 Initial Sequence: TCCAGCCTTGGC... (-9.625)
# Round 1: -15.61 WCCAGCCTGGGY T:10.0(83.85%),B:2.6(1.69%),P:1e-6
# Round 2: -15.61 WCCAGCCTGGGY T:10.0(83.85%),B:2.6(1.69%),P:1e-6
# =Final=: -14.37 WCCAGCCTGGGY T:5.0(83.33%),B:2.6(2.23%),P:1e-6
# Performing exhaustive masking of motif...
# Reprioritizing potential motifs...
# 6 of 25 Initial Sequence: TTTTGTATTTTT... (-9.625)
# Round 1: -11.91 TTTTGTATTTTT T:12.0(86.54%),B:6.2(5.00%),P:1e-5
# Round 2: -17.40 TTTTNTATTTTT T:11.0(83.85%),B:0.0(0.00%),P:1e-7
# Round 3: -17.40 TTTTNTATTTTT T:11.0(83.85%),B:0.0(0.00%),P:1e-7
# =Final=: -22.17 TTTTNTATTTTT T:6.0(100.00%),B:0.0(0.00%),P:1e-9
# Performing exhaustive masking of motif...
# Reprioritizing potential motifs...
# 7 of 25 Initial Sequence: GCATGAGCCACC... (-9.625)
# Round 1: -9.62 GCATGAGCCACC T:4.0(51.77%),B:0.0(0.00%),P:1e-4
# Round 2: -9.62 GCATGAGCCACC T:4.0(51.77%),B:0.0(0.00%),P:1e-4
# =Final=: -13.31 GCATGAGCCACC T:4.0(66.67%),B:0.0(0.00%),P:1e-5
# Performing exhaustive masking of motif...
# Reprioritizing potential motifs...
# Remaining seeds don't look promising (After initial 5 motifs, logp -6.215 > -6.778)
# 
# Finalizing Enrichment Statistics (new in v3.4)
# Reading input files...
# 32 total sequences read
# Cache length = 11180
# Using hypergeometric scoring
# Checking enrichment of 7 motif(s)
# |0%                                    50%                                  100%|
# =================================================================================
# Output in file: analysis_output//homerMotifs.motifs12
# 
# (Motifs in homer2 format)
# Determining similar motifs... 33 reduced to 19 motifs
# Outputing HTML and sequence logos for motif comparison...
# Checking de novo motifs against known motifs...
# Formatting HTML page...
# 1 of 19 (1e-9) similar to CRZ1(MacIsaac)/Yeast(0.704)
# 2 of 19 (1e-9) similar to HNRNPH2(RRM)/Homo_sapiens-RNCMPT00160-PBM/HughesRNA(0.747)
# 3 of 19 (1e-9) similar to EBF1(EBF)/Near-E2A-ChIP-Seq(GSE21512)/Homer(0.724)
# 4 of 19 (1e-9) similar to RLR1?/SacCer-Promoters/Homer(0.770)
# 5 of 19 (1e-7) similar to NCU08034(RRM)/Neurospora_crassa-RNCMPT00209-PBM/HughesRNA(0.670)
# 6 of 19 (1e-7) similar to Gfi1b/MA0483.1/Jaspar(0.784)
# 7 of 19 (1e-7) similar to CRZ1/MA0285.1/Jaspar(0.793)
# 8 of 19 (1e-5) similar to Klf4(Zf)/mES-Klf4-ChIP-Seq(GSE11431)/Homer(0.789)
# 9 of 19 (1e-5) similar to MATALPHA2/MA0328.2/Jaspar(0.720)
# 10 of 19 (1e-5) similar to CRZ1(MacIsaac)/Yeast(0.796)
# 11 of 19 (1e-5) similar to CG33714(RRM)/Drosophila_melanogaster-RNCMPT00009-PBM/HughesRNA(0.739)
# 12 of 19 (1e-5) similar to ttk/dmmpmm(SeSiMCMC)/fly(0.781)
# 13 of 19 (1e-5) similar to SRSF2(RRM)/Homo_sapiens-RNCMPT00072-PBM/HughesRNA(0.749)
# 14 of 19 (1e-4) similar to MOT2/MA0379.1/Jaspar(0.700)
# 15 of 19 (1e-4) similar to ZNF384/MA1125.1/Jaspar(0.960)
# 16 of 19 (1e-4) similar to slbo/dmmpmm(Bergman)/fly(0.744)
# 17 of 19 (1e-4) similar to Ct/dmmpmm(Noyes_hd)/fly(0.764)
# 18 of 19 (1e-4) similar to ESE3(AP2EREBP)/col-ESE3-DAP-Seq(GSE60143)/Homer(0.700)
# 19 of 19 (1e-3) similar to vnd/dmmpmm(Papatsenko)/fly(0.776)
# Job finished