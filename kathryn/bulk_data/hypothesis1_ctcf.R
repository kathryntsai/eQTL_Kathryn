# Ned to find download ChiCP data
# ==========================================================
# HYPOTHESIS 1: HI-C CTCF Modeling
# ==========================================================

# Papers from Friday 6/14/19

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3612258
# Find Hi-C data: look at HiGlass or some other dataset// don't use Hi-C data for T-cell or B-cell data specifically

# for closest trans eQTLs: SNP-Gene pairs --> find closest HI-C pair
# of resaonbly close ones (100 bp??): 1 snp 1 gene 1 trans eQTL pair, see where SNp matches in both columns. is the trans gene on the different opposite chromosome
# anti-pair: part of pair that's not the snp, ask if position is close to gene releavnt to snp

# Can't use these because they're HG38: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3612258
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3612257
# Would like to find dataset: hic peripheral blood hg19
# But eQTLGen (http://www.eqtlgen.org/publications.html) uses data from https://www.chicp.org/chicp/ -- try to download information
# https://www.ncbi.nlm.nih.gov/gds?LinkName=geoprofiles_gds&from_uid=54443598

# Can StringDB/CTCF mediate Hi-C relations?