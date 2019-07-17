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

# replace all parentheses - not sure if this will work in previous gsub creations
full_chr_mod <- gsub('c\\(|\\)', "", full_chr)
full_chr_mod <- gsub(')', "", full_chr_mod)

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
# write.fasta(a_gene_snp$fasta_seq_mod, paste("chromosome:GRCh37:",a_gene_snp$SNPChr, ":", a_gene_snp$SNPWindowStart, ":", a_gene_snp$SNPWindowEnd,":1", sep=""), "q2_3_output/fasta_master1.fa", open = "w", nbchar = 60, as.string = F)

write.fasta(a_gene_snp$fasta_seq_mod[1], paste("chromosome:GRCh37:",a_gene_snp$SNPChr[1], ":", a_gene_snp$SNPWindowStart[1], ":", a_gene_snp$SNPWindowEnd[1],":1", sep=""), "q2_3_output/fasta_master_alt.fa", open = "w", nbchar = 60, as.string = F)
for (i in 2:nrow(a_gene_snp)){
  write.fasta(a_gene_snp$fasta_seq_mod[i], paste("chromosome:GRCh37:",a_gene_snp$SNPChr[i], ":", a_gene_snp$SNPWindowStart[i], ":", a_gene_snp$SNPWindowEnd[i],":1", sep=""), "q2_3_output/fasta_master_alt.fa", open = "a", nbchar = 60, as.string = F)
}

write.fasta(a_gene_snp$fasta_seq_og[1], paste("chromosome:GRCh37:",a_gene_snp$SNPChr[1], ":", a_gene_snp$SNPWindowStart[1], ":", a_gene_snp$SNPWindowEnd[1],":1", sep=""), "q2_3_output/fasta_master_og.fa", open = "w", nbchar = 60, as.string = F)
for (i in 2:nrow(a_gene_snp)){
  write.fasta(a_gene_snp$fasta_seq_og[i], paste("chromosome:GRCh37:",a_gene_snp$SNPChr[i], ":", a_gene_snp$SNPWindowStart[i], ":", a_gene_snp$SNPWindowEnd[i],":1", sep=""), "q2_3_output/fasta_master_og.fa", open = "a", nbchar = 60, as.string = F)
}

# new fasta file! get rid of duplicates and make sure SNPs fit

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
# z <- cbind(str_sub(a_gene_snp$fasta_seq_og, 21, 21), a_gene_snp$AssessedAllele, a_gene_snp$OtherAllele, a_gene_snp$SNP, a_gene_snp$SNPWindowStart, a_gene_snp$SNPWindowEnd, a_gene_snp$SNPPos, str_sub(full_chr_mod,a_gene_snp$SNPPos, a_gene_snp$SNPPos), a_gene_snp$fasta_seq_og)
# colnames(z) <- c("str_sub(a_gene_snp$fasta_seq_og, 21, 21)", "a_gene_snp$AssessedAllele", "a_gene_snp$OtherAllele", "a_gene_snp$SNP", "a_gene_snp$SNPWindowStart", "a_gene_snp$SNPWindowEnd", "a_gene_snp$SNPPos", "str_sub(full_chr_mod,a_gene_snp$SNPPos,a_gene_snp$SNPPos", "a_gene_snp$fasta_sequence_og")
# zz <- z[which(z[,2]==toupper(z[,7]) & toupper(z[,7])==toupper(z[,1])), ]
# colnames(zz) <- c("str_sub(a_gene_snp$fasta_seq_og, 20, 20)", "a_gene_snp$AssessedAllele", "a_gene_snp$SNP", "a_gene_snp$SNPWindowStart", "a_gene_snp$SNPWindowEnd", "a_gene_snp$SNPPos", "str_sub(full_chr_mod,a_gene_snp$SNPPos,a_gene_snp$SNPPos", "a_gene_snp$fasta_sequence_og")
# zzz <- z[which(z[,2]!=toupper(z[,7])), ]
# colnames(zzz) <- c("str_sub(a_gene_snp$fasta_seq_og, 20, 20)", "a_gene_snp$AssessedAllele", "a_gene_snp$SNP", "a_gene_snp$SNPWindowStart", "a_gene_snp$SNPWindowEnd", "a_gene_snp$SNPPos", "str_sub(full_chr_mod,a_gene_snp$SNPPos,a_gene_snp$SNPPos", "a_gene_snp$fasta_sequence_og")

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

# ref[which(ref$"FASTA_ID_Abb" == "chromosome:GRCh37:1:1429976:1430016:1" & ref$`Motif Name` == "n-Myc(bHLH)/mES-nMyc-ChIP-Seq(GSE11431)/Homer"),]

# nrow(in_alt_not_ref)
# [1] 556934
# > nrow(in_ref_not_alt)
# [1] 535066

#in_ref_not_alt_unique <- in_ref_not_alt[!duplicated(in_ref_not_alt$FASTA_ID_Abb),] 
#in_alt_not_ref_unique <- in_alt_not_ref[!duplicated(in_alt_not_ref$FASTA_ID_Abb),] 
# fwrite(in_ref_not_alt_unique, "q2_3_output/in_ref_not_alt_unique.txt", sep="\t")
# fwrite(in_alt_not_ref_unique, "q2_3_output/in_alt_not_ref_unique.txt", sep="\t")

# > nrow(in_ref_not_alt_unique)
# [1] 70550
# > nrow(in_alt_not_ref_unique)
# [1] 72078

# not necessary
a_gene_snp <- fread("q2_3_output/a_gene_snp.txt", sep="\t")
a_gene_snp_new <- data.table(a_gene_snp[,"SNP"], paste("chromosome:GRCh37:",a_gene_snp$SNPChr,":", a_gene_snp$SNPWindowStart,":", a_gene_snp$SNPWindowEnd,":1",sep=""))
colnames(a_gene_snp_new)[2] <- "FASTA_ID_Abb"
a_gene_snp_new <- a_gene_snp_new[!duplicated(a_gene_snp_new),]
in_ref_not_alt_SNP <- right_join(a_gene_snp_new, in_ref_not_alt) # still 556934
in_alt_not_ref_SNP <- right_join(a_gene_snp_new, in_alt_not_ref) # still 535066
in_ref_not_alt_SNP <- in_ref_not_alt_SNP[!duplicated(in_ref_not_alt_SNP),]
in_alt_not_ref_SNP <- in_alt_not_ref_SNP[!duplicated(in_alt_not_ref_SNP),]
fwrite(in_ref_not_alt_SNP, "q2_3_output/in_ref_not_alt_SNP.txt", sep="\t")
fwrite(in_alt_not_ref_SNP, "q2_3_output/in_alt_not_ref_SNP.txt", sep="\t")

# > nrow(in_ref_not_alt_unique_SNP)
# [1] 264905
# > nrow(in_alt_not_ref_unique_SNP)
# [1] 269853

# dim(in_ref_not_alt_unique_SNP[!duplicated(in_ref_not_alt_unique_SNP),])
# [1] 70550     8
# dim(in_alt_not_ref_unique_SNP[!duplicated(in_alt_not_ref_unique_SNP),])
# [1] 72078     8
# This checks out!

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


