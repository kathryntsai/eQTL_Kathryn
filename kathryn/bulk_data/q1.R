# ==========================================================
# BULK DATA
# ==========================================================

# ==========================================================
# QUESTION 1
# ==========================================================

# How many SNP-gene associations exist in both bulk eqtl datasets? 
library(data.table)
library(dplyr)
library(glue)
# Read Emma Davenport's eQTL R-variables - not used
load("input_data/files_for_pfizer_eqtl.rda")

# Read Emma Davenport's Pfzier Data
davenport_data <- read.csv("input_data/pfizer_eqtl_table_1.csv")
# Read Lude Franke's cis-eQTL data
franke_cis_data <- fread("input_data/cis-eQTL_significant_20181017.txt", header = TRUE, sep = "\t", dec = ".")

# inner_join (only keep observations that are similar) by SNPS+GENE
franke_cis_gene_snps <- unique(data.table(franke_cis_data[,8], franke_cis_data[,2])) #10,527,168x2
davenport_gene_snps <- data.table(sub("*\\.[0-9]", "", davenport_data[,"Ensembl_ID"]), davenport_data[,"SNP"]) #4818x2
colnames(davenport_gene_snps) <- c("Gene", "SNP")
franke_davenport_gene_snps_common <- inner_join(davenport_gene_snps, franke_cis_gene_snps) #2962
# proportions:
nrow(franke_davenport_gene_snps_common)/nrow(davenport_gene_snps) # = 61.47779% proportion!
nrow(franke_davenport_gene_snps_common)/nrow(franke_cis_gene_snps) # = 0.02813672% proportion!

# ==========================================================
# QUESTION 3
# ==========================================================
# What proportion of eQTLgen trans-eQTLs are also cis-eQTLs (with a nearby gene?)

franke_cis_snps <- unique(data.table(franke_cis_data[,2])) #10,527,168x2
# Read Lude Franke's trans-eQTL data
franke_trans_data <- fread("input_data/trans-eQTL_significant_20181017.txt", header = TRUE, sep = "\t", dec = ".")
franke_trans_snps <- unique(data.table(franke_trans_data[,2])) # 59786

# inner_join by SNPs
franke_cis_trans_common_snps <- inner_join(franke_trans_snps, franke_cis_snps, by="SNP") # 3281

nrow(franke_cis_trans_common_snps)/nrow(franke_trans_snps) # 85.15443% proportion!
nrow(franke_cis_trans_common_snps)/nrow(franke_cis_snps) # 0.08867992% proportion!

# NOTE: HOMER PATHS DON'T WORK ANYMORE DUE TO ORGANIZATION 6/11/19
# ==========================================================
# QUESTION 2_1
# ==========================================================
# What % of eQTLs (cis / trans consider separately) overlap 
# a motif or are in LD with a SNP that overlaps a motif? 
# If a SNP is not genotyped or poorly imputed (missing data or incorrect data), 
# we won’t be able to see it’s true association with gene expression. 
# Therefore, a SNP in tight LD with this poorly documented SNP (which is in reality the causal SNP) 
# might show an artificially strong association signal with expression. 

# # use homer to find trans eQTLs that also overlap cis eQTLs
# write.table(franke_cis_trans_common_snps, file="q2_1_written/franke_cis_trans_common.txt",row.names= FALSE,col.names=FALSE,quote=FALSE) # write SNPs to file
# 
# # match cis/trans snps to cis genes, keep only unique genes for Homer findmotifs.pl // do I need more specificity for gene?  i.e. gene.1, gene.2...?  Do i keep only unique genes?
# franke_cis_trans_common_cis_genes <- unique(data.table(right_join(franke_cis_gene_snps, franke_cis_trans_common_snps))[,"Gene"]) # 24490 merged -> 4605 unique
# write.table(franke_cis_trans_common_cis_genes, file="input_data_homer/q2_1_old/franke_cis_trans_common_cis_genes.txt",row.names= FALSE,col.names="Gene_ID",quote=FALSE) # write cis genes to file
# # match cis/trans snps to trans genes, keep only unique genes for Homer findmotifs.pl // same questions as for cis
# franke_cis_trans_common_trans_genes <- unique(data.table(right_join(data.table(franke_trans_data[,"SNP"], franke_trans_data[,"Gene"]), franke_cis_trans_common_snps))[,"Gene"]) #56038 merged -> 6059 unique
# write.table(franke_cis_trans_common_trans_genes, file="input_data_homer/q2_1_old/franke_cis_trans_common_trans_genes.txt",row.names= FALSE,col.names="Gene_ID",quote=FALSE) # write trans genes to file
# 
# # having issues getting homer to find file?  am i putting it in the right directory? FIXED - spell "franke" correctly
# # findMotifs.pl franke_cis_trans_common_trans_genes.txt human  motifResults_trans1/ -find data/knownTFs/vertebrates/known.motifs > /my_dir/trans_output.txt
# 
# trans_homer_results <- fread("input_data_homer/old/trans_output.txt",header = TRUE, sep = "\t", dec = ".") #336,445
#   
# # parse out after parentheses in motif name
# # create table from data + trans_genes
# # print list of factors
# 
# trans_homer_results$MotifNameAbbreviated <-  sub("*\\(.*", "", trans_homer_results$'Motif Name')
# # trans_homer_results_organized <- data.table(trans_homer_results[,"Offset"],trans_homer_results[,"Strand"],trans_homer_results[,"MotifScore"],trans_homer_results[,"Ensembl"],trans_homer_results[,"MotifNameAbbreviated"]) #336,445
# # franke_trans_gene_snps <- data.table(franke_trans_data$SNP, franke_trans_data$Gene) #59,786
# # colnames(franke_trans_gene_snps) <- c("SNP", "Ensembl")
# # trans_homer_results_organized_snps <- left_join(trans_homer_results_organized, franke_trans_gene_snps) #left_join only yields 3,359,095, inner_join only yields 3,359,095 results
# # colnames(trans_homer_results_organized_snps) <- c("Offset", "Strand", "MotifScore", "Trans_Gene", "MotifNameAbbreviated", "SNP")
# # trans_homer_results_organized_snps_cis_gene_snps <- left_join(trans_homer_results_organized_snps, franke_cis_gene_snps) #left_join only yields 30,472,391, inner_join only yields 30,118,485
# # colnames(trans_homer_results_organized_snps_cis_gene_snps) <- c("Offset", "Strand", "MotifScore", "Trans_Gene", "MotifNameAbbreviated", "SNP", "Cis_Gene")
# 
# # WITHOUT STRAND
# trans_homer_results_organized <- data.table(trans_homer_results[,"Offset"],trans_homer_results[,"MotifScore"],trans_homer_results[,"Ensembl"],trans_homer_results[,"MotifNameAbbreviated"]) #336,445
# franke_trans_gene_snps <- data.table(franke_trans_data$SNP, franke_trans_data$Gene) #59,786
# colnames(franke_trans_gene_snps) <- c("SNP", "Ensembl")
# trans_homer_results_organized_snps <- left_join(trans_homer_results_organized, franke_trans_gene_snps) #left_join only yields 3,359,095, inner_join only yields 3,359,095 results
# colnames(trans_homer_results_organized_snps) <- c("Offset", "MotifScore", "Trans_Gene", "MotifNameAbbreviated", "SNP")
# trans_homer_results_organized_snps_cis_gene_snps <- left_join(trans_homer_results_organized_snps, franke_cis_gene_snps) #left_join only yields 30,472,391, inner_join only yields 30,118,485
# colnames(trans_homer_results_organized_snps_cis_gene_snps) <- c("Offset", "MotifScore", "Trans_Gene", "MotifNameAbbreviated", "SNP", "Cis_Gene")
# 
# write.csv(trans_homer_results_organized_snps_cis_gene_snps, file="input_data_homer/q2_1_old/trans_homer_results_organized_snps_cis_gene_snps.csv")
# # WRITE CODE TO GET UNIQUE ONES - can't run on local computer
# 
# franke_cis_trans_common_snps_gene <- inner_join(franke_trans_gene_snps, franke_cis_gene_snps, by="SNP")
