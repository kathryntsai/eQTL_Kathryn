# ==========================================================
# Q4
# ==========================================================

# only use SNPs that are JUST cis or JUST trans
only_cis <- left_join(franke_cis_data, franke_trans_data, by=c("SNP", "SNPChr", "SNPPos", "AssessedAllele", "OtherAllele"))
fwrite(onlycis, "")
only_trans <- right_join(franke_cis_data, franke_trans_data, by=c("SNP", "SNPChr", "SNPPos", "AssessedAllele", "OtherAllele"))