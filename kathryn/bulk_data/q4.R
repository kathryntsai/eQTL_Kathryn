# ==========================================================
# Q4
# ==========================================================

# only use SNPs that are JUST cis or JUST trans
only_cis <- left_join(franke_cis_data, franke_trans_data, by=c("SNP", "SNPChr", "SNPPos", "AssessedAllele", "OtherAllele"))
only_trans <- right_join(franke_cis_data, franke_trans_data, by=c("SNP", "SNPChr", "SNPPos", "AssessedAllele", "OtherAllele"))
fwrite(only_cis, "q4_output/only_cis.txt", sep="\t")
fwrite(only_trans, "q4_output/only_trans.txt", sep="\t")

only_cis <- fread("q4_output/only_cis.txt", sep="\t")
only_cis_downsampled <- only_cis %>% sample_n(200)
fwrite(only_cis_downsampled, "q4_output/only_cis_downsampled.txt", sep="\t")

only_trans <- fread("q4_output/only_trans.txt", sep="\t")
only_trans_downsampled <- only_trans %>% sample_n(200)
fwrite(only_trans_downsampled, "q4_output/only_trans_downsampled.txt", sep="\t")
