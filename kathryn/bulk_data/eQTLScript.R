#How many SNP-gene associations exist in both bulk eqtl datasets? 
library(data.table)
library(dplyr)
load("files_for_pfizer_eqtl.rda")
franke_cis_data <- fread("cis-eQTL_significant_20181017.txt", header = TRUE, sep = "\t", dec = ".")

# SNPS+GENE
franke_cis_snps <- unique(data.table(franke_cis_data[,8], franke_cis_data[,2])) #10,527,168x2
davenport_snps <- data.table(sub("*\\.[0-9]", "", pairs.int.tss[,1]), pairs.int.tss[,2]) #4976x2
colnames(davenport_snps) <- c("Gene", "SNP")
franke_davenport_common <- inner_join(davenport_snps, franke_cis_snps) #2887
#proportions:
nrow(franke_davenport_common)/nrow(davenport_snps) # = 58.01849% proportion!
nrow(franke_davenport_common)/nrow(franke_cis_snps) # = 0.0002742428 proportion!

# # JUST SNPS
# franke_snps <- unique(franke_cis_data[,2]) #10,527,168 --> 3,699,823
# davenport_snps <- data.table(pairs.int.tss[,2]) #4976
# colnames(davenport_snps) <- c("SNP")
# franke_davenport_common <- inner_join(davenport_snps, franke_snps) #42 # merge gives 24,255 which is not right
# #proportions:
# nrow(franke_davenport_common)/nrow(davenport_snps)
# nrow(franke_davenport_common)/nrow(franke_snps)