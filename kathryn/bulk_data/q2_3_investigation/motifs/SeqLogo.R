# SeqLogo - Tiffany's code

## IMPORTING INFORMATION, SETTING WORKING DIRECTORY
# setwd("/Users/amariuta/documents/SRLAB/Emma_Homer")
# 
source("https://bioconductor.org/biocLite.R")
biocLite("seqLogo")

library(seqLogo)
# 
# IRF1_jaspar <- read.table("IRF1.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE) 
# IRF2_jaspar <- read.table("IRF2.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE) 
# ISRE_jaspar <- read.table("ISRE.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE) 
# stat4_jaspar <- read.table("stat4.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# tbet_jaspar <- read.table("Tbet_Th1.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# 
# 
# IRF2_jaspar <- IRF2_jaspar[2:nrow(IRF2_jaspar),]
# stat4_jaspar <- stat4_jaspar[1:(nrow(stat4_jaspar)-3),]
# 
# #/data/srlab/amariuta/jobs/Emma_Homer/HomerFiles/IFN_damp_nlen0/knownResults/known9.motif
# IRF1_enriched <- read.table("IRF1_enriched.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE) 
# #data/srlab/amariuta/jobs/Emma_Homer/HomerFiles/IFN_damp_nlen0/knownResults/known1.motif
# IRF2_enriched <- read.table("IRF2_enriched.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE) 
# #/data/srlab/amariuta/jobs/Emma_Homer/HomerFiles/IFN_damp_nlen0/knownResults/known8.motif
# ISRE_enriched <- read.table("ISRE_enriched.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE) 
# #/data/srlab/amariuta/jobs/Emma_Homer/HomerFiles/IFN_mag_dedup_nlen0/knownResults/known11.motif
# stat4_enriched <- read.table("stat4_enriched.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE) 
# #/data/srlab/amariuta/jobs/Emma_Homer/HomerFiles/Tcell_damp_nlen0_uniquersid/knownResults/known17.motif
# stat3_enriched <- read.table("stat3_enriched.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE) 
# 
# IRF4_motif <- read.table("known1.motif", sep = "\t", header = F, stringsAsFactors = FALSE) 
# 
# 
# setwd("/Users/amariuta/documents/SRLAB/SNPs_and_TFs")
# PU1IRF <- read.table("PU_1_IRF.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE) 
# 
# Tbet_homer <- read.table("Tbet.txt", sep = "\t", header = F, stringsAsFactors = FALSE)
# Gata3_homer <- read.table("Gata3.txt", sep = "\t", header = F, stringsAsFactors = FALSE)
# Bcl6_homer <- read.table("Bcl6.txt", sep = "\t", header = F, stringsAsFactors = FALSE)
# FOXP3_jaspar <- read.table("Foxp3.txt", sep = "\t", header = F, stringsAsFactors = FALSE)
# Runx1_homer <- read.table("Runx1.txt", sep = "\t", header = F, stringsAsFactors = FALSE)
# Runx1_homer <- Runx1_homer[,1:4]

## CREATING THRESHOLD

## not sure what this is
# logoddsthresh <- c(x1$ref_MotifScore[1], x2$alt_MotifScore[1])
# AUC_status1 <- c(0.657, 0.645, 0.727, 0.49, 0.852)
# cor.test(logoddsthresh, AUC_status1) #0.9244963, p-val = 0.0755

pe <- makePWM(t(Runx1_homer))
pdf("Runx1_bitsscore.pdf")
seqLogo(pe, ic.scale = T) #T : bits score, F : probability
dev.off()

pe <- makePWM(t(ISRE_jaspar))
pdf("ISRE_jaspar.pdf")
seqLogo(pe, ic.scale = T) #T : bits score, F : probability
dev.off()

pe <- makePWM(t(tbet_jaspar))
pdf("TBET_jaspar.pdf")
seqLogo(pe, ic.scale = T) #T : bits score, F : probability
dev.off()

stat4_enriched <- stat4_enriched[2:nrow(stat4_enriched),]

pe <- makePWM(t(IRF1_enriched))
pj <- makePWM(t(IRF1_jaspar))

pdf("IRF1e_bitsscore.pdf")
seqLogo(pe, ic.scale = T) #T : bits score, F : probability
dev.off()
pdf("IRF1j_bitsscore.pdf")
seqLogo(pj, ic.scale = T)
dev.off()

pdf("IRF1e_prob.pdf")
seqLogo(pe, ic.scale = F) #T : bits score, F : probability
dev.off()
pdf("IRF1j_prob.pdf")
seqLogo(pj, ic.scale = F)
dev.off()



pe <- makePWM(t(IRF2_enriched))
pj <- makePWM(t(IRF2_jaspar))

pdf("IRF2e_bitsscore.pdf")
seqLogo(pe, ic.scale = T) #T : bits score, F : probability
dev.off()
pdf("IRF2j_bitsscore.pdf")
seqLogo(pj, ic.scale = T)
dev.off()

pdf("IRF2e_prob.pdf")
seqLogo(pe, ic.scale = F) #T : bits score, F : probability
dev.off()
pdf("IRF2j_prob.pdf")
seqLogo(pj, ic.scale = F)
dev.off()



pe <- makePWM(t(ISRE_enriched))
pj <- makePWM(t(ISRE_jaspar))

pdf("ISREe_bitsscore.pdf")
seqLogo(pe, ic.scale = T) #T : bits score, F : probability
dev.off()
pdf("ISREj_bitsscore.pdf")
seqLogo(pj, ic.scale = T)
dev.off()

pdf("ISREe_prob.pdf")
seqLogo(pe, ic.scale = F) #T : bits score, F : probability
dev.off()
pdf("ISREj_prob.pdf")
seqLogo(pj, ic.scale = F)
dev.off()


pe <- makePWM(t(stat4_enriched))
pj <- makePWM(t(stat4_jaspar))

pdf("stat4e_bitsscore.pdf")
seqLogo(pe, ic.scale = T) #T : bits score, F : probability
dev.off()
pdf("stat4j_bitsscore.pdf")
seqLogo(pj, ic.scale = T)
dev.off()

pdf("stat4e_prob.pdf")
seqLogo(pe, ic.scale = F) #T : bits score, F : probability
dev.off()
pdf("stat4j_prob.pdf")
seqLogo(pj, ic.scale = F)
dev.off()