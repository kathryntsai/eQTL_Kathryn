# ==========================================================
# HYPOTHESIS 1: HI-C CTCF Modeling
# ==========================================================

HiC_data_2 <- fread("/Users/kathryntsai/OneDrive\ -\ Villanova\ University/College/2018-2019/Summer\ 2019/TFs_eQTLs_Research/Important_Files/JavierreCHICP/DATA_S1/PCHiC_peak_matrix_cutoff5.tsv", sep="\t", header=T) #728,838 x 30
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525

# ==========================================================
# Load full chr
# ==========================================================

full_chr <- readRDS("q1_input/full_chr.rds")
full_chr_mod <- gsub('c\\(|\\)', "", full_chr)
full_chr_mod <- gsub(')', "", full_chr_mod)

