## homer output for Kathryn

#setup: restricted ourselves to SNPs that are both cis and trans eQTLs. 
#did not restrict ourselves to cis Genes that are TFs. 
#goal: want to identify if any cis eQTL Genes are TFs that can bind in the promoter of the trans eQTL gene 
dat <- fread("~/Downloads/output_trans.txt", header = T, sep = "\t", dec = ".")
#336K by 16 

#table with SNP, cis Gene, trans Gene, list of Motifs 
TFs <- sapply(1:nrow(dat), function(x) strsplit(dat$`Motif Name`[x], split = "[(]")[[1]][1])
length(unique(TFs)) #407 

snps <- fread("~/Downloads/franke_cis_trans_common_snps_gene.txt", header = T, sep = ",")
snps <- snps[,-1]
colnames(snps) <- c("SNP","cisGene","transGene")
# 601288  x    3

franke_data <- fread("~/Documents/SRLAB/Kathryn/cis-eQTL_significant_20181017.txt.gz", header = T, sep = "\t", dec = ".")
common_cisgenename <- franke_data$GeneSymbol[match(snps$cisGene, franke_data$Gene)]

#did not restrict ourselves to cis Genes that are TFs.
#the only possibility for explaining this first regime is if the cis eQTL gene is a TF gene 
m <- match(toupper(common_cisgenename), toupper(unique(TFs)))
#important to match case! 
length(which(is.na(m))) #590719
length(which(!is.na(m))) #10569 (can reduce our search space to this)

w <- which(!is.na(m))
snps_interesting <- snps[w,]
common_cisgenename_interesting <- common_cisgenename[w]
#sanity check 
length(which(is.na(match(toupper(common_cisgenename_interesting), toupper(unique(TFs)))))) == 0

#########
#motif_matches <- sapply(1:nrow(snps), function(x) paste(TFs[which(dat$Ensembl == as.character(snps[x,3]))], collapse = ","))
matches <- c()
for (i in 1:nrow(snps_interesting)){
  if (i %% 1000 == 0){print(i)}
  matches[i] <- match(toupper(common_cisgenename_interesting[i]), toupper(TFs[which(dat$Ensembl == as.character(snps_interesting[i,3]))]))
}
length(which(!is.na(matches))) #458 

#snps_interesting_TFmatch <- cbind(snps_interesting, ifelse(is.na(matches), 0, 1))
#w <- which(snps_interesting_TFmatch$V2 == 1)
snps_interesting_TFmatch <- snps_interesting[which(!is.na(matches)),]
head(snps_interesting_TFmatch)
snps_interesting_TFmatch <- cbind(snps_interesting_TFmatch, common_cisgenename_interesting[which(!is.na(matches))])
colnames(snps_interesting_TFmatch)[4] <- "cisGene_commonname"
snps_interesting_TFmatch <- snps_interesting_TFmatch[,c(1,2,4,3)]
write.table(snps_interesting_TFmatch, "~/Documents/SRLAB/Kathryn/snps_interesting_TFmatch.txt", row.names = F, col.names = T)

