# what's going on with homer??
# ==========================================================
# HYPOTHESIS 6: Quads Analysis
# ==========================================================
franke_cis_data_modified <- franke_cis_data[,c("SNP", "Gene", "Zscore", "SNPChr", "SNPPos")] 
colnames(franke_cis_data_modified) = c("SNP", "cisGene", "cisZscore", "cisSNPChr", "cisSNPPos")
franke_trans_data_modified <- franke_trans_data[,c("SNP", "Gene", "Zscore", "SNPChr", "SNPPos")] 
colnames(franke_trans_data_modified) = c("SNP", "transGene", "transZscore", "transSNPChr", "transSNPPos")

franke_cis_trans_common_snps_gene_zscore <- inner_join(franke_cis_data_modified, franke_trans_data_modified, by="SNP")

#setup: restricted ourselves to SNPs that are both cis and trans eQTLs.
#did not restrict ourselves to cis Genes that are TFs.
#goal: want to identify if any cis eQTL Genes are TFs that can bind in the promoter of the trans eQTL gene
dat <- fread("q2_1_input/homer/trans_output.txt", header = T, sep = "\t", dec = ".")
#336K by 16

#table with SNP, cis Gene, trans Gene, list of Motifs
TFs <- sapply(1:nrow(dat), function(x) strsplit(dat$`Motif Name`[x], split = "[(]")[[1]][1])
length(unique(TFs)) #407

#snps <- fread("q2_1_input/franke_cis_trans_common_snps_gene.csv", header = T, sep = ",")
snps <- franke_cis_trans_common_snps_gene_zscore
# 601288  x    5

# franke_cis_data <- fread("input_data/cis-eQTL_significant_20181017.txt", header = T, sep = "\t", dec = ".")
common_cisgenename <- franke_cis_data$GeneSymbol[match(snps$cisGene, franke_cis_data$Gene)]

#did not restrict ourselves to cis Genes that are TFs.
#the only possibility for explaining this first regime is if the cis eQTL gene is a TF gene
m <- match(toupper(common_cisgenename), toupper(unique(TFs)))
#important to match case! # UPDATED LIST: these numbers don't match
length(which(is.na(m))) #590719# NOW 590847
length(which(!is.na(m))) #10569 (can reduce our search space to this) # NOW 10441

w <- which(!is.na(m))
snps_interesting <- snps[w,]
common_cisgenename_interesting <- common_cisgenename[w]
#sanity check
z <- match(toupper(common_cisgenename_interesting), toupper(unique(TFs)))
length(which(is.na(z))) == 0

#########
#motif_matches <- sapply(1:nrow(snps), function(x) paste(TFs[which(dat$Ensembl == as.character(snps[x,3]))], collapse = ","))
matches <- c()
for (i in 1:nrow(snps_interesting)){
  matches[i] <- match(toupper(common_cisgenename_interesting[i]), toupper(TFs[which(dat$Ensembl == as.character(snps_interesting[i,"transGene"]))]))
}
length(which(!is.na(matches))) #458, NOW 812

#snps_interesting_TFmatch <- cbind(snps_interesting, ifelse(is.na(matches), 0, 1))
#w <- which(snps_interesting_TFmatch$V2 == 1)
snps_interesting_TFmatch <- snps_interesting[which(!is.na(matches)),]
snps_interesting_TFmatch <- cbind(snps_interesting_TFmatch, common_cisgenename_interesting[which(!is.na(matches))])
#snps_interesting_TFmatch <- snps_interesting_TFmatch[,c(1,2,4,3)]
write.table(snps_interesting_TFmatch, "q2_1_written/snps_interesting_TFmatch.txt", row.names = F, col.names = T)

# ==========================================================
# HYPOTHESIS 6: QUADS
# ==========================================================

snps_interesting_TFmatch <- fread("q2_1_written/snps_interesting_TFmatch.txt")
colnames(snps_interesting_TFmatch)[1] <- "SNP1"
colnames(snps_interesting_TFmatch)[10] <- "cisGene_commonname"
# Table of Quads for Bed File
z <- merge(snps_interesting_TFmatch, franke_cis_data, by.x="transGene", by.y="Gene") # 1749034 x 7
colnames(z)[12] <- "SNP2"

# SNP1 is trans to transGene, cis to cisGene which is also a TF Gene
# SNP2 is cis to transGene, so close enough to be in the motif of the TF

quads_interesting <- z[!duplicated(z$transGene),] # 72 x 19

quads_interesting_bed_file <- data.table(paste0("chr",quads_interesting$SNPChr),
                                         quads_interesting$SNPPos - 20,
                                         quads_interesting$SNPPos + 20,
                                         quads_interesting$SNP2,
                                         0,
                                         0)
write.table(quads_interesting_bed_file, "hypothesis6_output/quads_interesting_bed_file.txt", col.names = F, row.names=F, quote=F, sep='\t')

# Extracting 3 sequences from chr1
# Extracting 3 sequences from chr2
# Extracting 4 sequences from chr3
# Extracting 1 sequences from chr4
# Extracting 34 sequences from chr6
# Extracting 1 sequences from chr10
# Extracting 1 sequences from chr11
# Extracting 13 sequences from chr12
# Extracting 1 sequences from chr13
# Extracting 3 sequences from chr16
# Extracting 4 sequences from chr17
# Extracting 3 sequences from chr20
# Extracting 1 sequences from chr21

# works: findMotifsGenome.pl quads_interesting_bed_file.txt hg19 output/ -find ~/homer/data/knownTFs/vertebrates/known.motifs > quads_interesting.csv
# doesn't work: findMotifsGenome.pl ~/Users/kathryntsai/OneDrive\ -\ Villanova\ University/College/2018-2019/Summer\ 2019/TFs_eQTLs_Research/RProjects/eQTLPart1_2/hypothesis6_output/quads_interesting_bed_file.txt hg19 output/ -find ~/homer/data/knownTFs/vertebrates/known.motifs > quads_interesting.txt

quads_interesting_analysis <- read.csv("hypothesis6_input/quads_interesting.csv", sep="\t") # 2913 x 6
#TFs <- sapply(1:nrow(quads_interesting_analysis), function(x) strsplit(quads_interesting_analysis$`Motif Name`[x], split = "[(]")[[1]][1])
unique_TFs <- unique(unique(sub("*\\(.*", "", quads_interesting_analysis$'Motif.Name'))) # 320 total TFs
length(unique_TFs) # 320 unique TFs
unique_SNPs <- unique(sub("*\\-.*", "", quads_interesting_analysis$'PositionID')) # 57 total SNPs
length(unique_SNPs) # 57 unique SNPs

# Rename Homer Output PositionID to SNP2, Abbreviate column
colnames(quads_interesting_analysis)[1] <- "SNP2"
quads_interesting_analysis$SNP2_Abbreviated <- sub("*\\-.*", "", quads_interesting_analysis$'SNP2')

# Select and reorganize columns from original Z file
zz <- z[,c("SNP1", "cisGene", "SNP2", "transGene", "cisZscore", "cisSNPChr", "cisSNPPos", "transZscore", "transSNPChr", "transSNPPos", "Zscore", "SNPChr", "SNPPos")] # 646241 x 13
colnames(zz) <- c("SNP1", "cisGene", "SNP2", "transGene", "cisZscore", "cisSNPChr", "cisSNPPos", "transZscore", "transSNPChr", "transSNPPos", "Zscore", "SNPChr", "SNPPos")
zzz <- merge(quads_interesting_analysis, zz, by.x = "SNP2_Abbreviated", by.y = "SNP2") # Why does 2913 x 6, 646241 x 13...yield 1627 x 18? Because SNP2 isn't right, it should be the abbreviated one!
# Merge information from homer output and original Z file
zzz <- zz[match(quads_interesting_analysis$SNP2_Abbreviated, zz$SNP2)] # 2913x13
zzzz <- zzz[!is.na(zzz$SNP1)] # 609 x 13

# write.table(zzz, "SNP2_Offset_Sequence_Motif.name_Strand_MotifScore_SNP1_cisGene_transGene_Positions_1627_rows.txt", row.names = F, col.names = T)
# write.table(quads_interesting_analysis, "quads_interesting_analysis_analyzed_2913_rows.txt", row.names = F, col.names = T)
# 
# write.table(z, "hypothesis6_output/z.csv", row.names = F, col.names = T, quote=F, sep = "\t")
# write.table(zz, "hypothesis6_output/zz.csv", row.names = F, col.names = T, quote=F, sep = "\t")
# write.table(zzz, "hypothesis6_output/zzz.csv", row.names = F, col.names = T, quote=F, sep = "\t")
# write.table(zzzz, "hypothesis6_output/zzzz.csv", row.names = F, col.names = T, quote=F, sep = "\t")

# unique_SNPs from quads_interesting_analysis # 57 unique SNPs
# [1] "rs1024467"       "rs1000778"       "chr3:46757116"   "rs1001007"       "chr21:44483233"  "rs1005455"      
# [7] "rs1042657"       "chr12:110719585" "rs1009382"       "rs1008438"       "chr3:129238996"  "rs1004439"      
# [13] "rs1005380"       "rs10047617"      "chr12:10588581"  "rs1012753"       "rs1041981"       "chr6:31238984"  
# [19] "chr6:31238135"   "chr6:29913074"   "rs1003376"       "rs1011546"       "rs10082244"      "rs10159433"     
# [25] "chr20:1585445"   "rs1000424"       "rs1045493"       "rs1008723"       "rs1004047"       "rs10083738"     
# [31] "rs1005902"       "rs10082928"      "rs10161486"      "rs11759553"      "chr6:31996524"   "chr6:31237779"  
# [37] "chr6:31237771"   "rs1042148"       "rs10013187"      "chr6:29910719"   "chr6:29910717"   "chr6:29910682"  
# [43] "rs1020848"       "chr6:31783755"   "rs10507484"      "chr12:10763175"  "rs10046213"      "rs10084359"     
# [49] "chr2:136148401"  "rs1014448"       "chr6:31783863"   "rs1051788"       "chr6:31325028"   "rs1007976"      
# [55] "rs1042127"       "rs10049468"      "chr6:31797272"  

# zzz[!duplicated(zzz$SNP2),"SNP2"] # 12 unique SNPs from joining
# [1] "rs1024467"      "rs1000778"      "chr3:46757116"  "rs1001007"      "rs1005455"      "rs1008438"      "rs1004439"      "rs10159433"     "rs10161486"    
# [10] "chr12:10763175" "chr6:31783863"  "chr6:31797272" 

# zz["chr21:44483233",]
# SNP1 cisGene SNP2      transGene cisZscore cisSNPChr cisSNPPos transZscore transSNPChr transSNPPos Zscore SNPChr SNPPos
# 1: <NA>    <NA> <NA> chr21:44483233        NA        NA        NA          NA          NA          NA     NA     NA     NA

# z["chr21:44483233",]
# transGene SNP1 cisGene cisZscore cisSNPChr cisSNPPos transZscore transSNPChr transSNPPos cisGene_commonname Pvalue SNP2 SNPChr SNPPos AssessedAllele OtherAllele
# 1: chr21:44483233 <NA>    <NA>        NA        NA        NA          NA          NA          NA               <NA>     NA <NA>     NA     NA           <NA>        <NA>
#   Zscore GeneSymbol GeneChr GenePos NrCohorts NrSamples FDR
# 1:     NA       <NA>      NA      NA        NA        NA  NA

# Match transgene to cisgene
common_transgenename <- franke_trans_data$GeneSymbol[match(zzzz$transGene, franke_trans_data$Gene)] # 1627
TFs <- sub("*\\(.*", "", zzzz$'Motif.Name') 
m <- match(toupper(common_transgenename), toupper(unique(TFs)))
length(which(is.na(m))) #1627
length(which(!is.na(m))) #0 

w <- which(!is.na(m))
snps_interesting <- snps[w,]
common_cisgenename_interesting <- common_cisgenename[w]
#sanity check
z <- match(toupper(common_cisgenename_interesting), toupper(unique(TFs)))
length(which(is.na(z))) == 0

#########
#motif_matches <- sapply(1:nrow(snps), function(x) paste(TFs[which(dat$Ensembl == as.character(snps[x,3]))], collapse = ","))
matches <- c()
for (i in 1:nrow(snps_interesting)){
  matches[i] <- match(toupper(common_cisgenename_interesting[i]), toupper(TFs[which(dat$Ensembl == as.character(snps_interesting[i,"transGene"]))]))
}
length(which(!is.na(matches))) #458, NOW 812

#snps_interesting_TFmatch <- cbind(snps_interesting, ifelse(is.na(matches), 0, 1))
#w <- which(snps_interesting_TFmatch$V2 == 1)
snps_interesting_TFmatch <- snps_interesting[which(!is.na(matches)),]
snps_interesting_TFmatch <- cbind(snps_interesting_TFmatch, common_cisgenename_interesting[which(!is.na(matches))])
#snps_interesting_TFmatch <- snps_interesting_TFmatch[,c(1,2,4,3)]
write.table(snps_interesting_TFmatch, "q2_1_written/snps_interesting_TFmatch.txt", row.names = F, col.names = T)

read.csv("hypothesis6/allgenes.txt")