# ==========================================================
# Q4
# ==========================================================

# # only use SNPs that are JUST cis or JUST trans
only_cis <- left_join(franke_cis_data, franke_trans_data, by=c("SNP", "SNPChr", "SNPPos", "AssessedAllele", "OtherAllele"))
only_trans <- right_join(franke_cis_data, franke_trans_data, by=c("SNP", "SNPChr", "SNPPos", "AssessedAllele", "OtherAllele"))
# fwrite(only_cis, "q4_output/only_cis.txt", sep="\t")
# fwrite(only_trans, "q4_output/only_trans.txt", sep="\t")
# 
# only_cis <- fread("q4_output/only_cis.txt", sep="\t")
# only_cis_downsampled <- only_cis %>% sample_n(200)
# fwrite(only_cis_downsampled, "q4_output/only_cis_downsampled.txt", sep="\t")
# 
# only_trans <- fread("q4_output/only_trans.txt", sep="\t")
# only_trans_downsampled <- only_trans %>% sample_n(200)
# fwrite(only_trans_downsampled, "q4_output/only_trans_downsampled.txt", sep="\t")

# ==========================================================
# CREATION OF BED FILE WITH EQTLGEN AND PFIZER DATA
# ==========================================================

franke_cis_bed_file <- data.table(paste0("chr",only_cis$SNPChr),
                                  only_cis$SNPPos - 20,
                                  only_cis$SNPPos + 20,
                                  only_cis$SNP, # is this right?
                                  0,
                                  0)
colnames(franke_cis_bed_file) <- c("c1","c2","c3","c4","c5","c6")
franke_cis_bed_file_small <- franke_cis_bed_file[!duplicated(franke_cis_bed_file$c4)]

franke_trans_bed_file <- data.table(paste0("chr",only_trans$SNPChr),
                                    only_trans$SNPPos - 20,
                                    only_trans$SNPPos + 20,
                                    only_trans$SNP, # is this right?
                                  0,
                                  0)
colnames(franke_trans_bed_file) <- c("c1","c2","c3","c4","c5","c6")
franke_trans_bed_file_small <- franke_trans_bed_file[!duplicated(franke_trans_bed_file$c4)]

franke_cis_bed_file_smaller <- franke_cis_bed_file_small %>% sample_n(1000)
write.table(franke_cis_bed_file_smaller, file="q4_output/franke_cis_bed_file_smaller.txt", col.names = F, row.names=F, quote=F, sep='\t')

franke_trans_bed_file_smaller <- franke_trans_bed_file_small %>% sample_n(1000)
write.table(franke_trans_bed_file_smaller, file="q4_output/franke_trans_bed_file_smaller.txt", col.names = F, row.names=F, quote=F, sep='\t')
