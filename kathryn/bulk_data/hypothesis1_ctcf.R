# ==========================================================
# HYPOTHESIS 1: HI-C CTCF Modeling
# ==========================================================

# Papers from Friday 6/14/19

# Can't use these because they're HG38:
# main one i was going to use - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3612257
# too small -  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3612258
# Would like to find dataset: hic peripheral blood hg19 - https://www.ncbi.nlm.nih.gov/gds/?term=Hi-C+peripheral+blood
# But eQTLGen (http://www.eqtlgen.org/publications.html) uses data from https://www.chicp.org/chicp/ -- try to download information
# https://www.ncbi.nlm.nih.gov/gds?LinkName=geoprofiles_gds&from_uid=54443598

# actually, just convert it from Hg38 to Hg19 using 

# Find Hi-C data: look at HiGlass or some other dataset// don't use Hi-C data for T-cell or B-cell data specifically
# for closest trans eQTLs: SNP-Gene pairs --> find closest HI-C pair
# of resaonbly close ones (100 bp??): 1 snp 1 gene 1 trans eQTL pair, see where SNp matches in both columns. is the trans gene on the different opposite chromosome
# anti-pair: part of pair that's not the snp, ask if position is close to gene releavnt to snp

# Can StringDB/CTCF mediate Hi-C relations?

# ===================================

# This is the Hg38 Peripheral Blood dataset
HiC_data_1 <- fread("/Users/kathryntsai/OneDrive\ -\ Villanova\ University/College/2018-2019/Summer\ 2019/TFs_eQTLs_Research/Important_Files/PeripheralBloodHiC/GSM3612257_Neutrophil_Ecoli_HiC_rep1_validPairs.txt", sep="\t", header=T)
# # This is the Hg38 ChicP dataset
HiC_data_2 <- fread("/Users/kathryntsai/OneDrive\ -\ Villanova\ University/College/2018-2019/Summer\ 2019/TFs_eQTLs_Research/Important_Files/JavierreCHICP/DATA_S1/PCHiC_peak_matrix_cutoff5.tsv", sep="\t", header=T) #728,838 x 30
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525

# This is the CHICP dataset
