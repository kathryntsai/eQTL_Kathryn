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
franke_cis_data <- fread("cis-eQTL_significant_20181017.txt", header = TRUE, sep = "\t", dec = ".")

# inner_join (only keep observations that are similar) by SNPS+GENE
franke_cis_gene_snps <- unique(data.table(franke_cis_data[,8], franke_cis_data[,2])) #10,527,168x2
davenport_gene_snps <- data.table(sub("*\\.[0-9]", "", pairs.int.tss[,1]), pairs.int.tss[,2]) #4976x2
colnames(davenport_gene_snps) <- c("Gene", "SNP")
franke_davenport_gene_snps_common <- inner_join(davenport_gene_snps, franke_cis_gene_snps) #2887
# proportions:
nrow(franke_davenport_gene_snps_common)/nrow(davenport_gene_snps) # = 58.01849% proportion!
nrow(franke_davenport_gene_snps_common)/nrow(franke_cis_gene_snps) # = 0.0002742428 proportion!

# ==========================================================
# QUESTION 2
# ==========================================================
# What proportion of eQTLgen trans-eQTLs are also cis-eQTLs (with a nearby gene?)
# What mechanisms are trans-eQTLs acting through? 
# Do these SNPs affect the expression of a TF (cis eQTL with a TF gene) which goes on to bind maybe targets genome-wide and have a trans-effect? 
# In this case there is no DNA-DNA contact. 
# Or, are trans-eQTLs acting by bringing two distal regions of DNA together by 
# recruiting chromatin remodelers or CTCF (to lock the DNA interaction in place)?

franke_cis_snps <- unique(data.table(franke_cis_data[,2])) #10,527,168x2
franke_trans_data <- fread("trans-eQTL_significant_20181017.txt", header = TRUE, sep = "\t", dec = ".")
franke_trans_snps <- unique(data.table(franke_trans_data[,2])) # 59786

# # inner_join by SNPs+GENE // I don't think it makes sense to do both SNP+Gene - there are 0 matches anyway.
# franke_cis_trans_common <- inner_join(franke_trans_snps, franke_cis_snps) # 0
# 0% proportions for both

# inner_join by SNPs
franke_cis_trans_common <- inner_join(franke_trans_snps, franke_cis_snps, by="SNP") # 601,288

nrow(franke_cis_trans_common)/nrow(franke_trans_snps) # 85.15443% proportion!
# at first: 1005.734% proportion but this doesn't make sense - since only SNPs are being taken into account, we only care about SNPs, not # Gene column.  
nrow(franke_cis_trans_common)/nrow(franke_cis_snps) # 0.08867992% proportion!
# at first: 5.711774% proportion but this doesn't make sense - since only SNPs are being taken into account, we only care about SNPs, not # Gene column.

  