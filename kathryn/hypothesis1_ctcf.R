# ==========================================================
# HYPOTHESIS 1: HI-C CTCF Modeling
# ==========================================================

HiC_data_2 <- fread("/Users/kathryntsai/OneDrive\ -\ Villanova\ University/College/2018-2019/Summer\ 2019/TFs_eQTLs_Research/Important_Files/JavierreCHICP/DATA_S1/PCHiC_peak_matrix_cutoff5.tsv", sep="\t", header=T) #728,838 x 30
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525

bait <- HiC_data_2[, c("baitChr", "baitStart", "baitEnd", "baitID", "baitName")]
oe <- HiC_data_2[, c("oeChr", "oeStart", "oeEnd", "oeID", "oeName")]