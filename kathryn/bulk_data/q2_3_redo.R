# ==========================================================
# QUESTION 2_3_Redo
# ==========================================================

a <- franke_cis_data[which(GeneChr==1), c('Gene', 'GeneSymbol', 'GeneChr', 'GenePos')]
a <- a[!duplicated(a),]

a_gene_snp <- franke_cis_data[which(GeneChr==1), c('Gene', 'GeneSymbol', 'GeneChr', 'GenePos', 'SNP', 'SNPChr', 'SNPPos')]
# check
# a_gene <- a_gene_snp[!duplicated(a_gene_snp$Gene),]

# read in chromosomes from hg19 genome, taken from homer/data/genomes/hg19 from q2_3_input folder
chr1 <- fread("q2_3_input/chr1.fa")
chr1_collapsed <- gsub(', |\"|\n', "", paste(chr1,collapse=""))
fwrite(as.list(chr1_collapsed), "chr1_character.txt")

a_gene_snp$SNPWindowStart <- a_gene_snp$SNPPos-20
a_gene_snp$SNPWindowEnd <- a_gene_snp$SNPPos+20

# ryan reynolds should've played gatsby
# jodie switched labs

HiC_data <- fread("/Users/kathryntsai/OneDrive\ -\ Villanova\ University/College/2018-2019/Summer\ 2019/TFs_eQTLs_Research/Important_Files/PeripheralBloodHiC/GSM3612257_Neutrophil_Ecoli_HiC_rep1_validPairs.txt")