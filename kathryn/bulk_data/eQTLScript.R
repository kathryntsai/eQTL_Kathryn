# ==========================================================
# BULK DATA
# ==========================================================

# ==========================================================
# QUESTION 1
# ==========================================================

# How many SNP-gene associations exist in both bulk eqtl datasets? 
library(data.table)
library(dplyr)
load("files_for_pfizer_eqtl.rda")
davenport_data <- read.csv("pfizer_eqtl_table_1.csv")
franke_cis_data <- fread("cis-eQTL_significant_20181017.txt", header = TRUE, sep = "\t", dec = ".")

# inner_join (only keep observations that are similar) by SNPS+GENE
franke_cis_gene_snps <- unique(data.table(franke_cis_data[,8], franke_cis_data[,2])) #10,527,168x2
davenport_gene_snps <- data.table(sub("*\\.[0-9]", "", davenport_data[,"Ensembl_ID"]), davenport_data[,"SNP"]) #4818x2
colnames(davenport_gene_snps) <- c("Gene", "SNP")
franke_davenport_gene_snps_common <- inner_join(davenport_gene_snps, franke_cis_gene_snps) #2962
# proportions:
nrow(franke_davenport_gene_snps_common)/nrow(davenport_gene_snps) # = 61.47779% proportion!
nrow(franke_davenport_gene_snps_common)/nrow(franke_cis_gene_snps) # = 0.02813672% proportion!

# ==========================================================
# QUESTION 2
# ==========================================================
# What proportion of eQTLgen trans-eQTLs are also cis-eQTLs (with a nearby gene?)

franke_cis_snps <- unique(data.table(franke_cis_data[,2])) #10,527,168x2
franke_trans_data <- fread("trans-eQTL_significant_20181017.txt", header = TRUE, sep = "\t", dec = ".")
franke_trans_snps <- unique(data.table(franke_trans_data[,2])) # 59786

# inner_join by SNPs
franke_cis_trans_common <- inner_join(franke_trans_snps, franke_cis_snps, by="SNP") # 3281

nrow(franke_cis_trans_common)/nrow(franke_trans_snps) # 85.15443% proportion!
nrow(franke_cis_trans_common)/nrow(franke_cis_snps) # 0.08867992% proportion!

# ==========================================================
# QUESTION 3
# ==========================================================
# What % of eQTLs (cis / trans consider separately) overlap 
# a motif or are in LD with a SNP that overlaps a motif? 
# If a SNP is not genotyped or poorly imputed (missing data or incorrect data), 
# we won’t be able to see it’s true association with gene expression. 
# Therefore, a SNP in tight LD with this poorly documented SNP (which is in reality the causal SNP) 
# might show an artificially strong association signal with expression. 

  