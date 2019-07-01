# ==========================================================
# QUESTION 2_3_Redo
# ==========================================================

# ==========================================================
# Create Chr1 Information
# ==========================================================

# read in chromosomes from hg19 genome, taken from homer/data/genomes/hg19 from q2_3_input folder
chr1 <- fread("q2_3_input/chr1.fa")
# all at once is too long so divide it in chunks
# x <- 4985013
# div <- seq_len(x)
# paste(div[x %% div == 0], collapse=", ")
# 1, 3, 11, 29, 33, 87, 319, 957, 5209, 15627, 57299, 151061, 171897, 453183, 1661671, 4985013
# so this didn't work chr1_new <- split(chr1, ceiling(seq_along(chr1)/5209))

# subset chr1
chosen_factor=171897 # need 29 arrays of 171897 rows of 50 chars each
counter=1
total_rows <- nrow(chr1) # 4985013
x <- matrix(0, nrow=total_rows/chosen_factor*2-2, ncol=1)
for (i in 1:nrow(x)){
  #print(i)
  if (i %% 2 == 1){
    x[i,1] <- chosen_factor * counter 
  }
  else{
    x[i,1] <- chosen_factor * counter + 1
    counter = counter + 1
  }
}

# each 
chr1_1 <- gsub(', |\"|\n', "", paste(chr1[1:x[1]],collapse="")) #should be 171897*50=8594850, then nchar(chr1_1) = 8594853
chr1_2 <- gsub(', |\"|\n', "", paste(chr1[x[2]:x[3]],collapse=""))
chr1_3 <- gsub(', |\"|\n', "", paste(chr1[x[4]:x[5]],collapse=""))
chr1_4 <- gsub(', |\"|\n', "", paste(chr1[x[6]:x[7]],collapse=""))
chr1_5 <- gsub(', |\"|\n', "", paste(chr1[x[8]:x[9]],collapse=""))
chr1_6 <- gsub(', |\"|\n', "", paste(chr1[x[10]:x[11]],collapse=""))
chr1_7 <- gsub(', |\"|\n', "", paste(chr1[x[12]:x[13]],collapse=""))
chr1_8 <- gsub(', |\"|\n', "", paste(chr1[x[14]:x[15]],collapse=""))
chr1_9 <- gsub(', |\"|\n', "", paste(chr1[x[16]:x[17]],collapse=""))
chr1_10 <- gsub(', |\"|\n', "", paste(chr1[x[18]:x[19]],collapse=""))
chr1_11 <- gsub(', |\"|\n', "", paste(chr1[x[20]:x[21]],collapse=""))
chr1_12 <- gsub(', |\"|\n', "", paste(chr1[x[22]:x[23]],collapse=""))
chr1_13 <- gsub(', |\"|\n', "", paste(chr1[x[24]:x[25]],collapse=""))
chr1_14 <- gsub(', |\"|\n', "", paste(chr1[x[26]:x[27]],collapse=""))
chr1_15 <- gsub(', |\"|\n', "", paste(chr1[x[28]:x[29]],collapse=""))
chr1_16 <- gsub(', |\"|\n', "", paste(chr1[x[30]:x[31]],collapse=""))
chr1_17 <- gsub(', |\"|\n', "", paste(chr1[x[32]:x[33]],collapse=""))
chr1_18 <- gsub(', |\"|\n', "", paste(chr1[x[34]:x[35]],collapse=""))
chr1_19 <- gsub(', |\"|\n', "", paste(chr1[x[36]:x[37]],collapse=""))
chr1_20 <- gsub(', |\"|\n', "", paste(chr1[x[38]:x[39]],collapse=""))
chr1_21 <- gsub(', |\"|\n', "", paste(chr1[x[40]:x[41]],collapse=""))
chr1_22 <- gsub(', |\"|\n', "", paste(chr1[x[42]:x[43]],collapse=""))
chr1_23 <- gsub(', |\"|\n', "", paste(chr1[x[44]:x[45]],collapse=""))
chr1_24 <- gsub(', |\"|\n', "", paste(chr1[x[46]:x[47]],collapse=""))
chr1_25 <- gsub(', |\"|\n', "", paste(chr1[x[48]:x[49]],collapse=""))
chr1_26 <- gsub(', |\"|\n', "", paste(chr1[x[50]:x[51]],collapse=""))
chr1_27 <- gsub(', |\"|\n', "", paste(chr1[x[52]:x[53]],collapse=""))
chr1_28 <- gsub(', |\"|\n', "", paste(chr1[x[54]:x[55]],collapse=""))
chr1_29 <- gsub(', |\"|\n', "", paste(chr1[x[56]:4985013],collapse=""))

# fwrite(as.list(chr1_1), "q2_3_output/chr1_1.txt")
# fwrite(as.list(chr1_2), "q2_3_output/chr1_2.txt")
# fwrite(as.list(chr1_3), "q2_3_output/chr1_3.txt")
# fwrite(as.list(chr1_4), "q2_3_output/chr1_4.txt")
# fwrite(as.list(chr1_5), "q2_3_output/chr1_5.txt")
# fwrite(as.list(chr1_6), "q2_3_output/chr1_6.txt")
# fwrite(as.list(chr1_7), "q2_3_output/chr1_7.txt")
# fwrite(as.list(chr1_8), "q2_3_output/chr1_8.txt")
# fwrite(as.list(chr1_9), "q2_3_output/chr1_9.txt")
# fwrite(as.list(chr1_10), "q2_3_output/chr1_10.txt")
# fwrite(as.list(chr1_11), "q2_3_output/chr1_11.txt")
# fwrite(as.list(chr1_12), "q2_3_output/chr1_12.txt")
# fwrite(as.list(chr1_13), "q2_3_output/chr1_13.txt")
# fwrite(as.list(chr1_14), "q2_3_output/chr1_14.txt")
# fwrite(as.list(chr1_15), "q2_3_output/chr1_15.txt")
# fwrite(as.list(chr1_16), "q2_3_output/chr1_16.txt")
# fwrite(as.list(chr1_17), "q2_3_output/chr1_17.txt")
# fwrite(as.list(chr1_18), "q2_3_output/chr1_18.txt")
# fwrite(as.list(chr1_19), "q2_3_output/chr1_19.txt")
# fwrite(as.list(chr1_20), "q2_3_output/chr1_20.txt")
# fwrite(as.list(chr1_21), "q2_3_output/chr1_21.txt")
# fwrite(as.list(chr1_22), "q2_3_output/chr1_22.txt")
# fwrite(as.list(chr1_23), "q2_3_output/chr1_23.txt")
# fwrite(as.list(chr1_24), "q2_3_output/chr1_24.txt")
# fwrite(as.list(chr1_25), "q2_3_output/chr1_25.txt")
# fwrite(as.list(chr1_26), "q2_3_output/chr1_26.txt")
# fwrite(as.list(chr1_27), "q2_3_output/chr1_27.txt")
# fwrite(as.list(chr1_28), "q2_3_output/chr1_28.txt")
# fwrite(as.list(chr1_29), "q2_3_output/chr1_29.txt")

saveRDS(chr1_1, 'q2_3_output/chr1_1.rds')  
saveRDS(chr1_2, 'q2_3_output/chr1_2.rds')  
saveRDS(chr1_3, 'q2_3_output/chr1_3.rds')  
saveRDS(chr1_4, 'q2_3_output/chr1_4.rds')  
saveRDS(chr1_5, 'q2_3_output/chr1_5.rds')  
saveRDS(chr1_6, 'q2_3_output/chr1_6.rds')  
saveRDS(chr1_7, 'q2_3_output/chr1_7.rds')  
saveRDS(chr1_8, 'q2_3_output/chr1_8.rds')  
saveRDS(chr1_9, 'q2_3_output/chr1_9.rds')  
saveRDS(chr1_10, 'q2_3_output/chr1_10.rds')  
saveRDS(chr1_11, 'q2_3_output/chr1_11.rds')  
saveRDS(chr1_12, 'q2_3_output/chr1_12.rds')  
saveRDS(chr1_13, 'q2_3_output/chr1_13.rds')  
saveRDS(chr1_14, 'q2_3_output/chr1_14.rds')  
saveRDS(chr1_15, 'q2_3_output/chr1_15.rds')  
saveRDS(chr1_16, 'q2_3_output/chr1_16.rds')  
saveRDS(chr1_17, 'q2_3_output/chr1_17.rds')  
saveRDS(chr1_18, 'q2_3_output/chr1_18.rds')  
saveRDS(chr1_19, 'q2_3_output/chr1_19.rds')  
saveRDS(chr1_20, 'q2_3_output/chr1_20.rds')  
saveRDS(chr1_21, 'q2_3_output/chr1_21.rds')  
saveRDS(chr1_22, 'q2_3_output/chr1_22.rds')  
saveRDS(chr1_23, 'q2_3_output/chr1_23.rds')  
saveRDS(chr1_24, 'q2_3_output/chr1_24.rds')  
saveRDS(chr1_25, 'q2_3_output/chr1_25.rds')  
saveRDS(chr1_26, 'q2_3_output/chr1_26.rds')  
saveRDS(chr1_27, 'q2_3_output/chr1_27.rds')  
saveRDS(chr1_28, 'q2_3_output/chr1_28.rds')  
saveRDS(chr1_29, 'q2_3_output/chr1_29.rds')


# ==========================================================
# Make full Chr1
# ==========================================================

# chr1_1 <- fread("q2_3_output/chr1_1.txt")
# chr1_2 <- fread("q2_3_output/chr1_2.txt")
# chr1_3 <- fread("q2_3_output/chr1_3.txt")
# chr1_4 <- fread("q2_3_output/chr1_4.txt")
# chr1_5 <- fread("q2_3_output/chr1_5.txt")
# chr1_6 <- fread("q2_3_output/chr1_6.txt")
# chr1_7 <- fread("q2_3_output/chr1_7.txt")
# chr1_8 <- fread("q2_3_output/chr1_8.txt")
# chr1_9 <- fread("q2_3_output/chr1_9.txt")
# chr1_10 <- fread("q2_3_output/chr1_10.txt")
# chr1_11 <- fread("q2_3_output/chr1_11.txt")
# chr1_12 <- fread("q2_3_output/chr1_12.txt")
# chr1_13 <- fread("q2_3_output/chr1_13.txt")
# chr1_14 <- fread("q2_3_output/chr1_14.txt")
# chr1_15 <- fread("q2_3_output/chr1_15.txt")
# chr1_16 <- fread("q2_3_output/chr1_16.txt")
# chr1_17 <- fread("q2_3_output/chr1_17.txt")
# chr1_18 <- fread("q2_3_output/chr1_18.txt")
# chr1_19 <- fread("q2_3_output/chr1_19.txt")
# chr1_20 <- fread("q2_3_output/chr1_20.txt")
# chr1_21 <- fread("q2_3_output/chr1_21.txt")
# chr1_22 <- fread("q2_3_output/chr1_22.txt")
# chr1_23 <- fread("q2_3_output/chr1_23.txt")
# chr1_24 <- fread("q2_3_output/chr1_24.txt")
# chr1_25 <- fread("q2_3_output/chr1_25.txt")
# chr1_26 <- fread("q2_3_output/chr1_26.txt")
# chr1_27 <- fread("q2_3_output/chr1_27.txt")
# chr1_28 <- fread("q2_3_output/chr1_28.txt")
# chr1_29 <- fread("q2_3_output/chr1_29.txt")

chr1_1 <- readRDS("q2_3_output/chr1_1.rds")
chr1_2 <- readRDS("q2_3_output/chr1_2.rds")
chr1_3 <- readRDS("q2_3_output/chr1_3.rds")
chr1_4 <- readRDS("q2_3_output/chr1_4.rds")
chr1_5 <- readRDS("q2_3_output/chr1_5.rds")
chr1_6 <- readRDS("q2_3_output/chr1_6.rds")
chr1_7 <- readRDS("q2_3_output/chr1_7.rds")
chr1_8 <- readRDS("q2_3_output/chr1_8.rds")
chr1_9 <- readRDS("q2_3_output/chr1_9.rds")
chr1_10 <- readRDS("q2_3_output/chr1_10.rds")
chr1_11 <- readRDS("q2_3_output/chr1_11.rds")
chr1_12 <- readRDS("q2_3_output/chr1_12.rds")
chr1_13 <- readRDS("q2_3_output/chr1_13.rds")
chr1_14 <- readRDS("q2_3_output/chr1_14.rds")
chr1_15 <- readRDS("q2_3_output/chr1_15.rds")
chr1_16 <- readRDS("q2_3_output/chr1_16.rds")
chr1_17 <- readRDS("q2_3_output/chr1_17.rds")
chr1_18 <- readRDS("q2_3_output/chr1_18.rds")
chr1_19 <- readRDS("q2_3_output/chr1_19.rds")
chr1_20 <- readRDS("q2_3_output/chr1_20.rds")
chr1_21 <- readRDS("q2_3_output/chr1_21.rds")
chr1_22 <- readRDS("q2_3_output/chr1_22.rds")
chr1_23 <- readRDS("q2_3_output/chr1_23.rds")
chr1_24 <- readRDS("q2_3_output/chr1_24.rds")
chr1_25 <- readRDS("q2_3_output/chr1_25.rds")
chr1_26 <- readRDS("q2_3_output/chr1_26.rds")
chr1_27 <- readRDS("q2_3_output/chr1_27.rds")
chr1_28 <- readRDS("q2_3_output/chr1_28.rds")
chr1_29 <- readRDS("q2_3_output/chr1_29.rds")

# paste("chr1_", 1:29, sep = "", collapse = ", ")
full_chr <- paste(chr1_1, chr1_2, chr1_3, chr1_4, chr1_5, chr1_6, chr1_7, chr1_8, chr1_9, chr1_10, chr1_11, chr1_12, chr1_13, chr1_14, chr1_15, chr1_16, chr1_17, chr1_18, chr1_19, chr1_20, chr1_21, chr1_22, chr1_23, chr1_24, chr1_25, chr1_26, chr1_27, chr1_28, chr1_29, sep = "", collapse = "")
saveRDS(full_chr, "q2_3_output/full_chr.rds")

# ==========================================================
# Load full chr
# ==========================================================

# paste("saveRDS(chr1_", 1:29, ", 'q2_3_output/chr1_", 1:29, ".rds')", sep = "", collapse="  ")
full_chr <- readRDS("q2_3_output/full_chr.rds")

# some random calculations to figure things out
# > nchar(gsub(', |\"|\n', "", paste(chr1[1:x[1]],collapse="")))
# [1] 8594853
# > nchar(full_chr)
# [1] 249250708
# > 249250708/8594853
# [1] 29
# > nchar("TCATGCCCCGAGAGCTGAGTGCAAGGGAGAGGCAGCGCTGTCTGTGCTTC")*171897
# [1] 8594850
# > 171897*50*29
# [1] 249250650
# > 249250708-249250650
# [1] 58
# > substr(full_chr, 1, 10)
# [1] "c(NNNNNNNN"
# > substr(chr1_1, 1, 10)
# [1] "c(NNNNNNNN"

# replace all parentheses - not sure if this will work in previous gsub creations
full_chr_mod <- gsub('c\\(|\\)', "", full_chr)
full_chr_mod <- gsub(')', "", full_chr_mod)
# write.table(full_chr_mod,file="chr1mod.txt")

# str_sub(full_chr_mod, 8594750, 8594900)
# [1] "Cgaaagtttgatgaagatgagatctacacagtcacaaagcaacctcccacaaaatggttattaattaaaaagggaaaaagagacaccgtacaatggcaaaa)cctggcagacatcacgtccagcaactggtcaaagtgaacattaccagca"
# > gsub('//)', "", str_sub(full_chr_mod, 8594750, 8594900))
# [1] "Cgaaagtttgatgaagatgagatctacacagtcacaaagcaacctcccacaaaatggttattaattaaaaagggaaaaagagacaccgtacaatggcaaaa)cctggcagacatcacgtccagcaactggtcaaagtgaacattaccagca"
# > gsub(')', "", str_sub(full_chr_mod, 8594750, 8594900))
# [1] "Cgaaagtttgatgaagatgagatctacacagtcacaaagcaacctcccacaaaatggttattaattaaaaagggaaaaagagacaccgtacaatggcaaaacctggcagacatcacgtccagcaactggtcaaagtgaacattaccagca"
# > full_chr_mod <- gsub('c\\(|)', "", full_chr)
# > 
#   > a_gene_snp$fasta_seq_og <- toupper(str_sub(full_chr_mod, a_gene_snp$SNPWindowStart+11, a_gene_snp$SNPWindowEnd+1+11))
# > # rs10082323 tggaaattgagccttggagAgattaaatgcatggggcatgcc
#   > a_gene_snp$fasta_seq_mod <- a_gene_snp$fasta_seq_og
# > substr(a_gene_snp$fasta_seq_mod, 20, 20) <- a_gene_snp$OtherAllele
# > head(ref_test_group_with_fasta)
# SNPWindow       SNP.x
# 1: chromosome:GRCh37:1:100378631:100378671:1   rs6577147
# 2: chromosome:GRCh37:1:100378631:100378671:1   rs6577147
# 3: chromosome:GRCh37:1:100479240:100479280:1    rs482778
# 4: chromosome:GRCh37:1:100479240:100479280:1    rs482778
# 5: chromosome:GRCh37:1:100479240:100479280:1    rs482778
# 6: chromosome:GRCh37:1:100870398:100870438:1 rs143494014
# FASTA_ID ref_Offset
# 1: chromosome:GRCh37:1:100378631:100378671:1-Dup1         11
# 2: chromosome:GRCh37:1:100378631:100378671:1-Dup1         11
# 3: chromosome:GRCh37:1:100479240:100479280:1-Dup2          6
# 4: chromosome:GRCh37:1:100479240:100479280:1-Dup2          6
# 5: chromosome:GRCh37:1:100479240:100479280:1-Dup2          6
# 6: chromosome:GRCh37:1:100870398:100870438:1-Dup1          0
# Sequence
# 1:         TAGGGTGTGGTT
# 2:         TAGGGTGTGGTT
# 3:      TGTGTATATATATAC
# 4:      TGTGTATATATATAC
# 5:      TGTGTATATATATAC
# 6: AAAAAACAAAACAAAACAAA
# Motif Name ref_Strand
# 1:      EKLF(Zf)/Erythrocyte-Klf1-ChIP-Seq(GSE20478)/Homer          +
#   2:      EKLF(Zf)/Erythrocyte-Klf1-ChIP-Seq(GSE20478)/Homer          +
#   3: OCT:OCT(POU,Homeobox)/NPC-OCT6-ChIP-Seq(GSE43916)/Homer          +
#   4: OCT:OCT(POU,Homeobox)/NPC-OCT6-ChIP-Seq(GSE43916)/Homer          +
#   5: OCT:OCT(POU,Homeobox)/NPC-OCT6-ChIP-Seq(GSE43916)/Homer          +
#   6: FOXA1:AR(Forkhead,NR)/LNCAP-AR-ChIP-Seq(GSE27824)/Homer          +
#   ref_MotifScore            Gene GeneSymbol GeneChr   GenePos       SNP.y
# 1:       11.35154 ENSG00000079335     CDC14A       1 100898208   rs6577147
# 2:       11.35154 ENSG00000156876      SASS6       1 100573815   rs6577147
# 3:       10.43584 ENSG00000079335     CDC14A       1 100898208    rs482778
# 4:       10.43584 ENSG00000122477     LRRC39       1 100629090    rs482778
# 5:       10.43584 ENSG00000156875      HIAT1       1 100526293    rs482778
# 6:       11.90027 ENSG00000079335     CDC14A       1 100898208 rs143494014
# SNPChr    SNPPos AssessedAllele OtherAllele SNPWindowStart SNPWindowEnd
# 1:      1 100378651              T           C      100378631    100378671
# 2:      1 100378651              T           C      100378631    100378671
# 3:      1 100479260              A           C      100479240    100479280
# 4:      1 100479260              A           C      100479240    100479280
# 5:      1 100479260              A           C      100479240    100479280
# 6:      1 100870418              A           C      100870398    100870438
# fasta_seq_og
# 1: ACAAAAATTAGTAGGGTGTGGTTATGTGAGCTTGTAGTCCCA
# 2: ACAAAAATTAGTAGGGTGTGGTTATGTGAGCTTGTAGTCCCA
# 3: TATATGTGTGTATATATATACACACACATACCTATTTTATAT
# 4: TATATGTGTGTATATATATACACACACATACCTATTTTATAT
# 5: TATATGTGTGTATATATATACACACACATACCTATTTTATAT
# 6: AAAAAACAAAACAAAACAAAACAAAACAAAACAAAAAACTAA
# fasta_seq_mod
# 1: ACAAAAATTAGTAGGGTGTCGTTATGTGAGCTTGTAGTCCCA
# 2: ACAAAAATTAGTAGGGTGTCGTTATGTGAGCTTGTAGTCCCA
# 3: TATATGTGTGTATATATATCCACACACATACCTATTTTATAT
# 4: TATATGTGTGTATATATATCCACACACATACCTATTTTATAT
# 5: TATATGTGTGTATATATATCCACACACATACCTATTTTATAT
# 6: AAAAAACAAAACAAAACAACACAAAACAAAACAAAAAACTAA
# > dim(full_chr_mod)
# NULL
# > length(full_chr_mod)
# [1] 1
# > nchar(full_chr_mod)
# [1] 249250621
# > full_chr_mod <- gsub('c\\(|\\)', "", full_chr)
# > 
#   > nchar(full_chr_mod)
# [1] 249250621
# > full_chr_mod <- gsub(')', "", full_chr_mod)
# > 
#   > nchar(full_chr_mod)
# [1] 249250621
# > stri_count(full_chr, regex="\\)")
# Error in stri_count(full_chr, regex = "\\)") : 
#   could not find function "stri_count"
# > library(stringi)
# Warning message:
#   package ‘stringi’ was built under R version 3.5.2 
# > stri_count(full_chr, regex="\\)")
# [1] 29
# > stri_count(full_chr, regex="\\(")
# [1] 29
# > stri_count(full_chr_mod, regex="\\(")
# [1] 0
# > stri_count(full_chr_mod, regex="(")
# Error in stri_count_regex(str, regex, ...) : 
#   Incorrectly nested parentheses in regexp pattern. (U_REGEX_MISMATCHED_PAREN)
# > str_sub(full_chr_mod, 8594750, 8594900)
# [1] "Cgaaagtttgatgaagatgagatctacacagtcacaaagcaacctcccacaaaatggttattaattaaaaagggaaaaagagacaccgtacaatggcaaaacctggcagacatcacgtccagcaactggtcaaagtgaacattaccagcaa"
# > str_sub(full_chr_mod, 8594000, 8594500)
# [1] "GGCTTTGTTTTCCAAAAATACATGTACTCTCAACCATCTCCAATATTTTCCATCTCAGGTAAATGTTTCTGAGAATTCAAATGTGATGGCAGATAAAAGCAAACAGGAAAAAATTTGTATCCAGCTCTAGCTGCTCATCTTCAAAAGCTACACAAAGACCACGAGCTGTGGTGAGTGGAAAAGGGAACTGTTACTCAAATGCAAGTAAATTTTAGAATGTCTCAAATTCAAATTAGTTGAAATACAGGTGTTTTAGATTGTTAACTTTTCCACTTAAGTATAAAGATTTATGTTGAGGAAAAAAATGGTCCCAAGGTCACAGAGCAAAGTGACAGAAAATTaagaaataaagacattaaaaataatgaaataaggtttcttaatgttaaagaagaatgttacacataaggaaaagggaataaactagcagaaattctgtagtactggattgaacttggagaatttcacatgaacttatgaatatacacacataaatacaaa"
# > stri_count(full_chr, regex="c\\(|\\)")
# [1] 58
# > dim(strsplit(full_chr_mod, ","))
# NULL
# > dim(strsplit(full_chr_mod, ")"))
# NULL
# > dim(strsplit(full_chr, ")"))
# NULL
# > dim(strsplit(full_chr, "\\)"))
# NULL
# nchar(full_chr_mod) 249250621

# ==========================================================
# a_gene_snp
# ==========================================================

# a <- franke_cis_data[which(GeneChr==1), c('Gene', 'GeneSymbol', 'GeneChr', 'GenePos')]
# a <- a[!duplicated(a),]

a_gene_snp <- franke_cis_data[which(GeneChr==1), c('Gene', 'GeneSymbol', 'GeneChr', 'GenePos', 'SNP', 'SNPChr', 'SNPPos', 'AssessedAllele', 'OtherAllele')]
# check
# a_gene <- a_gene_snp[!duplicated(a_gene_snp$Gene),]

a_gene_snp$SNPWindowStart <- a_gene_snp$SNPPos-20
a_gene_snp$SNPWindowEnd <- a_gene_snp$SNPPos+20

a_gene_snp$fasta_seq_og <- toupper(str_sub(full_chr_mod, a_gene_snp$SNPWindowStart, a_gene_snp$SNPWindowEnd+1))
a_gene_snp$fasta_seq_mod <- a_gene_snp$fasta_seq_og
# NEW based on previous test: make bp position 21, and make it Assessed Allele
substr(a_gene_snp$fasta_seq_mod, 21, 21) <- a_gene_snp$AssessedAllele

fwrite(a_gene_snp, "q2_3_output/a_gene_snp.txt", sep="\t")

# the 1 at the end should be the + or - strand...fix this
# doesn't really work
# write.fasta(a_gene_snp$fasta_seq_mod, paste("chromosome:GRCh37:",a_gene_snp$SNPChr, ":", a_gene_snp$SNPWindowStart, ":", a_gene_snp$SNPWindowEnd,":1", sep=""), "q2_3_output/fasta_master1.fa", open = "w", nbchar = 60, as.string = F)

write.fasta(a_gene_snp$fasta_seq_mod[1], paste("chromosome:GRCh37:",a_gene_snp$SNPChr[1], ":", a_gene_snp$SNPWindowStart[1], ":", a_gene_snp$SNPWindowEnd[1],":1", sep=""), "q2_3_output/fasta_master_alt.fa", open = "w", nbchar = 60, as.string = F)
for (i in 2:nrow(a_gene_snp)){
  write.fasta(a_gene_snp$fasta_seq_mod[i], paste("chromosome:GRCh37:",a_gene_snp$SNPChr[i], ":", a_gene_snp$SNPWindowStart[i], ":", a_gene_snp$SNPWindowEnd[i],":1", sep=""), "q2_3_output/fasta_master_alt.fa", open = "a", nbchar = 60, as.string = F)
}

write.fasta(a_gene_snp$fasta_seq_og[1], paste("chromosome:GRCh37:",a_gene_snp$SNPChr[1], ":", a_gene_snp$SNPWindowStart[1], ":", a_gene_snp$SNPWindowEnd[1],":1", sep=""), "q2_3_output/fasta_master_og.fa", open = "w", nbchar = 60, as.string = F)
for (i in 2:nrow(a_gene_snp)){
  write.fasta(a_gene_snp$fasta_seq_og[i], paste("chromosome:GRCh37:",a_gene_snp$SNPChr[i], ":", a_gene_snp$SNPWindowStart[i], ":", a_gene_snp$SNPWindowEnd[i],":1", sep=""), "q2_3_output/fasta_master_og.fa", open = "a", nbchar = 60, as.string = F)
}

# new fasta file!

a_gene_snp_new <- a_gene_snp[which(str_sub(a_gene_snp$fasta_seq_og, 21, 21) == a_gene_snp$OtherAllele), ] # 710,183 x 14
a_gene_snp_new2 <- a_gene_snp_new[!duplicated(a_gene_snp_new$SNP),] # 269,450 x 14 non duplicated

write.fasta(a_gene_snp_new2$fasta_seq_mod[1], paste("chromosome:GRCh37:",a_gene_snp_new2$SNPChr[1], ":", a_gene_snp_new2$SNPWindowStart[1], ":", a_gene_snp_new2$SNPWindowEnd[1],":1", sep=""), "q2_3_output/fasta_master_alt.fa", open = "w", nbchar = 60, as.string = F)
for (i in 2:nrow(a_gene_snp_new2)){
  write.fasta(a_gene_snp_new2$fasta_seq_mod[i], paste("chromosome:GRCh37:",a_gene_snp_new2$SNPChr[i], ":", a_gene_snp_new2$SNPWindowStart[i], ":", a_gene_snp_new2$SNPWindowEnd[i],":1", sep=""), "q2_3_output/fasta_master_alt.fa", open = "a", nbchar = 60, as.string = F)
}

write.fasta(a_gene_snp_new2$fasta_seq_og[1], paste("chromosome:GRCh37:",a_gene_snp_new2$SNPChr[1], ":", a_gene_snp_new2$SNPWindowStart[1], ":", a_gene_snp_new2$SNPWindowEnd[1],":1", sep=""), "q2_3_output/fasta_master_og.fa", open = "w", nbchar = 60, as.string = F)
for (i in 2:nrow(a_gene_snp_new2)){
  write.fasta(a_gene_snp_new2$fasta_seq_og[i], paste("chromosome:GRCh37:",a_gene_snp_new2$SNPChr[i], ":", a_gene_snp_new2$SNPWindowStart[i], ":", a_gene_snp_new2$SNPWindowEnd[i],":1", sep=""), "q2_3_output/fasta_master_og.fa", open = "a", nbchar = 60, as.string = F)
}


# test:
z <- cbind(str_sub(a_gene_snp$fasta_seq_og, 21, 21), a_gene_snp$AssessedAllele, a_gene_snp$OtherAllele, a_gene_snp$SNP, a_gene_snp$SNPWindowStart, a_gene_snp$SNPWindowEnd, a_gene_snp$SNPPos, str_sub(full_chr_mod,a_gene_snp$SNPPos, a_gene_snp$SNPPos), a_gene_snp$fasta_seq_og)
colnames(z) <- c("str_sub(a_gene_snp$fasta_seq_og, 21, 21)", "a_gene_snp$AssessedAllele", "a_gene_snp$OtherAllele", "a_gene_snp$SNP", "a_gene_snp$SNPWindowStart", "a_gene_snp$SNPWindowEnd", "a_gene_snp$SNPPos", "str_sub(full_chr_mod,a_gene_snp$SNPPos,a_gene_snp$SNPPos", "a_gene_snp$fasta_sequence_og")
zz <- z[which(z[,2]==toupper(z[,7]) & toupper(z[,7])==toupper(z[,1])), ]
colnames(zz) <- c("str_sub(a_gene_snp$fasta_seq_og, 20, 20)", "a_gene_snp$AssessedAllele", "a_gene_snp$SNP", "a_gene_snp$SNPWindowStart", "a_gene_snp$SNPWindowEnd", "a_gene_snp$SNPPos", "str_sub(full_chr_mod,a_gene_snp$SNPPos,a_gene_snp$SNPPos", "a_gene_snp$fasta_sequence_og")
zzz <- z[which(z[,2]!=toupper(z[,7])), ]
colnames(zzz) <- c("str_sub(a_gene_snp$fasta_seq_og, 20, 20)", "a_gene_snp$AssessedAllele", "a_gene_snp$SNP", "a_gene_snp$SNPWindowStart", "a_gene_snp$SNPWindowEnd", "a_gene_snp$SNPPos", "str_sub(full_chr_mod,a_gene_snp$SNPPos,a_gene_snp$SNPPos", "a_gene_snp$fasta_sequence_og")

# http://bioinfo.cipf.es/apps-beta/bam-viewer/0.0.1/ use this to check allele stuff
# > head(a_gene_snp)
# Gene GeneSymbol GeneChr   GenePos        SNP SNPChr
# 1: ENSG00000000457      SCYL3       1 169842606 rs10082323      1
# 2: ENSG00000000457      SCYL3       1 169842606  rs1011266      1
# 3: ENSG00000000457      SCYL3       1 169842606 rs10127867      1
# 4: ENSG00000000457      SCYL3       1 169842606 rs10157246      1
# 5: ENSG00000000457      SCYL3       1 169842606 rs10157266      1
# 6: ENSG00000000457      SCYL3       1 169842606 rs10157398      1
# SNPPos AssessedAllele OtherAllele SNPWindowStart
# 1: 169643176              A           G      169643156
# 2: 169647046              A           G      169647026
# 3: 170018276              C           T      170018256
# 4: 169756203              G           A      169756183
# 5: 169756389              G           A      169756369
# 6: 169756548              C           G      169756528
# SNPWindowEnd                               fasta_seq_og
# 1:    169643196 AGATTAAATGCATGGGGCATGCCATTTGACTAGAAACTGGAA
# 2:    169647066 GACTATTTTTCCTTCTTGCCGATTTTTATCTGGTTTTTAAAT
# 3:    170018296 TCTGGGATGATCGAGCTTGGTTGGGGGAGGTGGGTCTGCCAT
# 4:    169756223 AGAACTATGCTGCTGCTGCTGCAGTGTAGCCAGGACGCACAG
# 5:    169756409 CTTCTGTAGCCCTAATTTCCGGTTCAAACTCTGCATTCACCT
# 6:    169756568 CAAGAATCCCCACCTCAAAAGTCACTATCTCCCTCCCTGGTA
# fasta_seq_mod
# 1: AGATTAAATGCATGGGGCAGGCCATTTGACTAGAAACTGGAA
# 2: GACTATTTTTCCTTCTTGCGGATTTTTATCTGGTTTTTAAAT
# 3: TCTGGGATGATCGAGCTTGTTTGGGGGAGGTGGGTCTGCCAT
# 4: AGAACTATGCTGCTGCTGCAGCAGTGTAGCCAGGACGCACAG
# 5: CTTCTGTAGCCCTAATTTCAGGTTCAAACTCTGCATTCACCT
# 6: CAAGAATCCCCACCTCAAAGGTCACTATCTCCCTCCCTGGTA
# > head(z)
# [,1] [,2] [,3]         [,4]        [,5]        [,6]       
# [1,] "T"  "A"  "rs10082323" "169643156" "169643196" "169643176"
# [2,] "C"  "A"  "rs1011266"  "169647026" "169647066" "169647046"
# [3,] "G"  "C"  "rs10127867" "170018256" "170018296" "170018276"
# [4,] "T"  "G"  "rs10157246" "169756183" "169756223" "169756203"
# [5,] "C"  "G"  "rs10157266" "169756369" "169756409" "169756389"
# [6,] "A"  "C"  "rs10157398" "169756528" "169756568" "169756548"
# [,7]
# [1,] "g" 
# [2,] "g" 
# [3,] "t" 
# [4,] "G" 
# [5,] "G" 
# [6,] "g" 
# > View(a_gene_snp$fasta_seq_og)
# > head(z <- cbind(str_sub(a_gene_snp$fasta_seq_og, 20, 20), a_gene_snp$AssessedAllele, a_gene_snp$SNP, a_gene_snp$SNPWindowStart, a_gene_snp$SNPWindowEnd, a_gene_snp$SNPPos, str_sub(full_chr_mod,a_gene_snp$SNPPos, a_gene_snp$SNPPos), a_gene_snp$fasta_sequence_og))
# [,1] [,2] [,3]         [,4]        [,5]        [,6]       
# [1,] "T"  "A"  "rs10082323" "169643156" "169643196" "169643176"
# [2,] "C"  "A"  "rs1011266"  "169647026" "169647066" "169647046"
# [3,] "G"  "C"  "rs10127867" "170018256" "170018296" "170018276"
# [4,] "T"  "G"  "rs10157246" "169756183" "169756223" "169756203"
# [5,] "C"  "G"  "rs10157266" "169756369" "169756409" "169756389"
# [6,] "A"  "C"  "rs10157398" "169756528" "169756568" "169756548"
# [,7]
# [1,] "g" 
# [2,] "g" 
# [3,] "t" 
# [4,] "G" 
# [5,] "G" 
# [6,] "g" 
# > 
#   > View(a_gene_snp)
# > "AGACCCCCTGCCAGGGAGACCCCAGGCGCCTGAATGGCCACG"
# [1] "AGACCCCCTGCCAGGGAGACCCCAGGCGCCTGAATGGCCACG"
# > "TGCCAGGGAGACTCCAGGCGCCTGA"
# [1] "TGCCAGGGAGACTCCAGGCGCCTGA"
# > nchar("TGCCAGGGAGACTCCAGGCGCCTGA")
# [1] 25
# > substr("AGACCCCCTGCCAGGGAGACCCCAGGCGCCTGAATGGCCACG", 20, 20)
# [1] "C"
# > substr("AGACCCCCTGCCAGGGAGACCCCAGGCGCCTGAATGGCCACG", 20, 20)
# [1] "C"
# > substr("AGACCCCCTGCCAGGGAGACCCCAGGCGCCTGAATGGCCACG", 18, 22)
# [1] "GACCC"
# > 
#   > substr("AGATTAAATGCATGGGGCATGCCATTTGACTAGAAACTGGAA", 20, 20)
# [1] "T"
# > substr("AGATTAAATGCATGGGGCATGCCATTTGACTAGAAACTGGAA", 19, 21)
# [1] "ATG"
# > substr("AGATTAAATGCATGGGGCATGCCATTTGACTAGAAACTGGAA", 21, 21)
# [1] "G"
# > substr("TTTCTTTACCAGCACTCTTAATAATCATTCTTCTGTGCATCT", 21, 21)
# [1] "A"
# > substr("TTTCTTTACCAGCACTCTTAATAATCATTCTTCTGTGCATCT", 20, 20)
# [1] "A"
# > substr('CAAGAATCCCCACCTCAAAAGTCACTATCTCCCTCCCTGGTA', 21, 21)
# [1] "G"


# ==========================================================
# create homer output
# ==========================================================

# findMotifs.pl /Users/kathryntsai/OneDrive\ -\ Villanova\ University/College/2018-2019/Summer\ 2019/TFs_eQTLs_Research/RProjects/Project2/q2_3_output/fasta_master1.fa fasta /Users/kathryntsai/OneDrive\ -\ Villanova\ University/College/2018-2019/Summer\ 2019/TFs_eQTLs_Research/RProjects/Project2/q2_3_output/ > output.txt
# findMotifs.pl /Users/kathryntsai/OneDrive\ -\ Villanova\ University/College/2018-2019/Summer\ 2019/TFs_eQTLs_Research/RProjects/Project2/q2_3_output/fasta_master_alt.fa fasta /Users/kathryntsai/OneDrive\ -\ Villanova\ University/College/2018-2019/Summer\ 2019/TFs_eQTLs_Research/RProjects/Project2/q2_3_output/redo -find ~/homer/data/knownTFs/vertebrates/known.motifs > homer_master_alt.txt
# findMotifs.pl /Users/kathryntsai/OneDrive\ -\ Villanova\ University/College/2018-2019/Summer\ 2019/TFs_eQTLs_Research/RProjects/Project2/q2_3_output/fasta_master_og.fa fasta /Users/kathryntsai/OneDrive\ -\ Villanova\ University/College/2018-2019/Summer\ 2019/TFs_eQTLs_Research/RProjects/Project2/q2_3_output/redo -find ~/homer/data/knownTFs/vertebrates/known.motifs > homer_master_ref.txt

# ==========================================================
# analyze homer output
# ==========================================================

ref <- fread("q2_3_input/homer_master_ref.txt", sep='\t', header=T)
colnames(ref) <- c("FASTA_ID", "ref_Offset", "Sequence", "Motif Name", "ref_Strand", "ref_MotifScore")
alt <- fread("q2_3_input/homer_master_alt.txt", sep='\t', header=T)
colnames(alt) <- c("FASTA_ID", "alt_Offset", "Sequence", "Motif Name", "alt_Strand", "alt_MotifScore")

# > nrow(ref)
# [1] 6144198
# > nrow(alt)
# [1] 6178496

#ref_alt_common <- merge(ref, alt, by=c("Sequence", "Motif Name"))
#fwrite(ref_alt_common, "q2_3_output/ref_alt_common.txt", sep="\t")

# > nrow(ref_uniqueTFs <-ref[!duplicated(ref$"Motif Name"),])
# [1] 426
# > nrow(alt_uniqueTFs <- alt[!duplicated(alt$"Motif Name"),])
# [1] 426

# not by sequence!
#in_ref_not_alt <- anti_join(ref, alt, by=c("Sequence", "Motif Name"))
#in_ref_not_alt <- ref[!(ref$"Motif Name" %in% alt$"Motif Name" & ref$"Sequence" %in% alt$"Sequence"),]
#in_alt_not_ref <- anti_join(alt, ref, by=c("Sequence", "Motif Name"))
#in_ref_not_alt <- alt[!(alt$"Motif Name" %in% ref$"Motif Name"),]

# ~20,000 each

# in_ref_not_alt <- anti_join(ref, alt, by="Motif Name")
# in_alt_not_ref <- anti_join(alt, ref, by="Motif Name")

# 0
# 0 

ref$FASTA_ID_Abb <- substr(ref$FASTA_ID, 1, regexpr("\\-", ref$FASTA_ID)-1) # didn't work: gsub(" \\-.*", "", ref$FASTA_ID)
alt$FASTA_ID_Abb <- substr(alt$FASTA_ID, 1, regexpr("\\-", alt$FASTA_ID)-1) # didn't work: gsub(" \\-.*", "", alt$FASTA_ID)
colnames(ref) <- c("FASTA_ID", "ref_Offset", "Sequence", "Motif Name", "ref_Strand", "ref_MotifScore", "FASTA_ID_Abb")
colnames(alt) <- c("FASTA_ID", "alt_Offset", "Sequence", "Motif Name", "alt_Strand", "alt_MotifScore", "FASTA_ID_Abb")
in_ref_not_alt <- anti_join(ref, alt, by=c("FASTA_ID_Abb", "Motif Name"))
in_alt_not_ref <- anti_join(alt, ref, by=c("FASTA_ID_Abb", "Motif Name"))
# fwrite(in_ref_not_alt, "q2_3_output/in_ref_not_alt_by_motif_name_fastaid.txt", sep="\t")
# fwrite(in_alt_not_ref, "q2_3_output/in_alt_not_ref_by_motif_name_fastaid.txt", sep="\t")

# nrow(in_alt_not_ref)
# [1] 556934
# > nrow(in_ref_not_alt)
# [1] 535066

in_ref_not_alt_unique <- in_ref_not_alt[!duplicated(in_ref_not_alt$FASTA_ID_Abb),] 
in_alt_not_ref_unique <- in_alt_not_ref[!duplicated(in_alt_not_ref$FASTA_ID_Abb),] 
# fwrite(in_ref_not_alt_unique, "q2_3_output/in_ref_not_alt_unique.txt", sep="\t")
# fwrite(in_alt_not_ref_unique, "q2_3_output/in_alt_not_ref_unique.txt", sep="\t")

# > nrow(in_ref_not_alt_unique)
# [1] 70550
# > nrow(in_alt_not_ref_unique)
# [1] 72078

# not necessary
a_gene_snp <- fread("q2_3_output/a_gene_snp.txt", sep="\t")
a_gene_snp_new <- data.table(a_gene_snp[,"SNP"], paste("chromosome:GRCh37:",a_gene_snp$SNPChr,":", a_gene_snp$SNPWindowStart,":", a_gene_snp$SNPWindowEnd,":1",sep=""))
colnames(a_gene_snp_new)[2] <- "SNPWindow"
in_ref_not_alt_unique_SNP <- merge(a_gene_snp_new, in_ref_not_alt_unique, by.x="SNPWindow", by.y="FASTA_ID_Abb", all.y=T) # inner_join doesn't work?
in_alt_not_ref_unique_SNP <- merge(a_gene_snp_new, in_alt_not_ref_unique, by.x="SNPWindow", by.y="FASTA_ID_Abb", all.y=T) # inner_join doesn't work?
fwrite(in_ref_not_alt_unique_SNP, "q2_3_output/in_ref_not_alt_unique_SNP.txt", sep="\t")
fwrite(in_alt_not_ref_unique_SNP, "q2_3_output/in_alt_not_ref_unique_SNP.txt", sep="\t")

# > nrow(in_ref_not_alt_unique_SNP)
# [1] 264905
# > nrow(in_alt_not_ref_unique_SNP)
# [1] 269853

# dim(in_ref_not_alt_unique_SNP[!duplicated(in_ref_not_alt_unique_SNP),])
# [1] 70550     8
# dim(in_alt_not_ref_unique_SNP[!duplicated(in_alt_not_ref_unique_SNP),])
# [1] 72078     8
# This checks out!

# ==========================================================
# read homer output
# ==========================================================

in_ref_not_alt_unique_SNP <- fread("q2_3_output/in_ref_not_alt_unique_SNP.txt")
in_alt_not_ref_unique_SNP <- fread("q2_3_output/in_ref_not_alt_unique_SNP.txt")

x1 <- in_ref_not_alt_unique_SNP[!duplicated(in_ref_not_alt_unique_SNP),] # 70550
x2 <- in_alt_not_ref_unique_SNP[!duplicated(in_alt_not_ref_unique_SNP),] # 72078

# can you pull out a couple examples of TFs (maybe 1 or 2) from these 73K and
# identify which part of the motif is affected by the SNP? I bet it
# will be the part of a motif where there is a very heavy weight on one
# particular nucleotide and messing that up results in failure to identify the
# motif with the alternative allele. also, could you show one or two examples of
# TFs from the 69K set, where searching the reference genome fails to identify a
# motif but searching with the alternative allele finds it? Similarly, I bet the
# SNPs occur in the super important core part of the motif with high
# probabilities (from the PWM) for that nucleotide. Do you have experience
# making sequence logos? they’re fun plots to show the
# important of each nucleotide in a motif


