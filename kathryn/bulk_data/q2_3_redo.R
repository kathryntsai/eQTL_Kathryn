# ==========================================================
# QUESTION 2_3_Redo
# ==========================================================

a <- franke_cis_data[which(Gene==c('ENSG00000127191', 
                                   'ENSG00000214279', 
                                   'ENSG00000168765', 
                                   'ENSG00000125966',
                                   'ENSG00000258531',
                                   'ENSG00000237437'
)), c('Gene', 'GeneSymbol', 'GeneChr', 'GenePos')]
a <- a[!duplicated(a),]

a_gene_snp <- franke_cis_data[which(Gene==c('ENSG00000127191', 
                                            'ENSG00000214279', 
                                            'ENSG00000168765', 
                                            'ENSG00000125966',
                                            'ENSG00000258531',
                                            'ENSG00000237437'
)), c('Gene', 'GeneSymbol', 'GeneChr', 'GenePos', 'SNP', 'SNPChr', 'SNPPos')]
a_gene <- a_gene_snp[!duplicated(a_gene_snp$Gene),]

# read in chromosomes from hg19 genome, taken from homer/data/genomes/hg19
chr1 <- read.table('q2_3_input/chr1.fa')