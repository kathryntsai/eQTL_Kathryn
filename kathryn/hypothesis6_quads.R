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

quads_interesting <- z[!duplicated(z$transGene),] # 257 x 23

quads_interesting_bed_file <- data.table(paste0("chr",quads_interesting$SNPChr),
                                         quads_interesting$SNPPos - 20,
                                         quads_interesting$SNPPos + 20,
                                         quads_interesting$SNP2,
                                         0,
                                         0)
write.table(quads_interesting_bed_file, "hypothesis6_output/quads_interesting_bed_file.txt", col.names = F, row.names=F, quote=F, sep='\t')

# Version 1
# Position file = quads_interesting_bed_file.txt
# Genome = hg19
# Output Directory = output/
#   Will find motif(s) in /Users/kathryntsai/homer/data/knownTFs/vertebrates/known.motifs
# Found mset for "human", will check against vertebrates motifs
# Peak/BED file conversion summary:
#   BED/Header formatted lines: 72
# peakfile formatted lines: 0
# 
# Peak File Statistics:
#   Total Peaks: 72
# Redundant Peak IDs: 15
# Peaks lacking information: 0 (need at least 5 columns per peak)
# Peaks with misformatted coordinates: 0 (should be integer)
# Peaks with misformatted strand: 0 (should be either +/- or 0/1)
# 
# Redunant Peaks found: Remove or rename these or some programs may have trouble...
# 
# 15 duplicate peak IDs out of 72 total peaks
# Background files for 200 bp fragments found.
# Custom genome sequence directory: /Users/kathryntsai/homer/.//data/genomes/hg19//
#   
# Extracting sequences from directory: /Users/kathryntsai/homer/.//data/genomes/hg19//
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
# 
# 
# Reading input files...
# 72 total sequences read
# 428 motifs loaded
# Finding instances of 428 motif(s)
# |0%                                    50%                                  100%|
#   =================================================================================
#   Cleaning up tmp files...


# Version 2
# Position file = quads_interesting_bed_file.txt
# Genome = hg19
# Output Directory = output/
#   Will find motif(s) in /Users/kathryntsai/homer/data/knownTFs/vertebrates/known.motifs
# Found mset for "human", will check against vertebrates motifs
# Peak/BED file conversion summary:
#   BED/Header formatted lines: 257
# peakfile formatted lines: 0
# 
# Peak File Statistics:
#   Total Peaks: 257
# Redundant Peak IDs: 3
# Peaks lacking information: 0 (need at least 5 columns per peak)
# Peaks with misformatted coordinates: 0 (should be integer)
# Peaks with misformatted strand: 0 (should be either +/- or 0/1)
# 
# Redunant Peaks found: Remove or rename these or some programs may have trouble...
# 
# 3 duplicate peak IDs out of 257 total peaks
# Background files for 200 bp fragments found.
# Custom genome sequence directory: /Users/kathryntsai/homer/.//data/genomes/hg19//
#   
#   Extracting sequences from directory: /Users/kathryntsai/homer/.//data/genomes/hg19//
#   Extracting 32 sequences from chr1
# Extracting 17 sequences from chr2
# Extracting 14 sequences from chr3
# Extracting 9 sequences from chr4
# Extracting 14 sequences from chr5
# Extracting 12 sequences from chr6
# Extracting 13 sequences from chr7
# Extracting 8 sequences from chr8
# Extracting 11 sequences from chr9
# Extracting 11 sequences from chr10
# Extracting 14 sequences from chr11
# Extracting 21 sequences from chr12
# Extracting 2 sequences from chr13
# Extracting 7 sequences from chr14
# Extracting 5 sequences from chr15
# Extracting 11 sequences from chr16
# Extracting 16 sequences from chr17
# Extracting 5 sequences from chr18
# Extracting 15 sequences from chr19
# Extracting 10 sequences from chr20
# Extracting 1 sequences from chr21
# Extracting 9 sequences from chr22
# 
# 
# Reading input files...
# 257 total sequences read
# 428 motifs loaded
# Finding instances of 428 motif(s)
# |0%                                    50%                                  100%|
#   =================================================================================
#   Cleaning up tmp files...
# 

# works: findMotifsGenome.pl quads_interesting_bed_file.txt hg19 output/ -find ~/homer/data/knownTFs/vertebrates/known.motifs > quads_interesting.csv
# doesn't work: findMotifsGenome.pl ~/Users/kathryntsai/OneDrive\ -\ Villanova\ University/College/2018-2019/Summer\ 2019/TFs_eQTLs_Research/RProjects/eQTLPart1_2/hypothesis6_output/quads_interesting_bed_file.txt hg19 output/ -find ~/homer/data/knownTFs/vertebrates/known.motifs > quads_interesting.txt

quads_interesting_analysis <- read.csv("hypothesis6_input/quads_interesting.csv", sep="\t") # 10562 x 6
#TFs <- sapply(1:nrow(quads_interesting_analysis), function(x) strsplit(quads_interesting_analysis$`Motif Name`[x], split = "[(]")[[1]][1])
unique_TFs <- unique(sub("*\\(.*", "", quads_interesting_analysis$'Motif.Name')) # 386 total TFs

# glue_collapse(unique_TFs, sep=", ")
# AP-1, AP-2gamma, AP-2alpha, Ap4, FOXA1:AR, ARE, AR-halfsite, Arnt:Ahr, Ascl1, Ascl2, Atf1, Atf2, Atf3, Atf4, Atf7, Atoh1, Bach1, Bach2, Bapx1, Barx1, IRF:BATF, BATF, Bcl11a, Bcl6, BHLHA15, bHLHE40, bHLHE41, BMAL1, BMYB, BORIS, Brachyury, Brn1, Brn2, bZIP:IRF, Cdx2, CDX4, CEBP:AP1, CEBP:CEBP, CEBP, Chop, CHR, CLOCK, c-Myc, COUP-TFII, CRE, CRX, CTCF, CUX1, Cux2, Dlx3, Mouse_Recombination_Hotspot, DMRT1, DMRT6, RAR:RXR, DUX4, Duxbl, E2A, E2F1, E2F3, E2F4, E2F6, E2F7, EBF2, EBF, EBF1, EBNA1, E-box, Egr1, Egr2, EHF, EKLF, ELF1, ELF3, Elf4, ELF5, Elk1, Elk4, Eomes, ERE, ERG, Erra, ERRg, Unknown-ESC-element, Esrrb, ETS1, Ets1-distal, ETS:E-box, ETS, ETS:RUNX, ETV1, Etv2, ETV4, EWS:ERG-fusion, EWS:FLI1-fusion, Fli1, Fosl2, FOXA1, Foxa2, Foxa3, FoxD3, Fox:Ebox, Foxf1, Foxh1, FOXK1, FOXK2, FoxL2, FOXM1, Foxo1, Foxo3, FOXP1, Fra1, Fra2, FXR, GABPA, Gata2, GATA3, GATA, Gata4, Gata6, Gata1, Gfi1b, GFY, Gli2, GLI3, GLIS3, GRE, GRHL2, GSC, Hand2, HEB, HIC1, HIF-1a, HIF-1b, HIF2a, HINFP, HLF, HNF1b, Hnf1, HNF4a, Hnf6b, HNF6, Hoxa10, Hoxa11, Hoxa13, HOXA1, HOXA2, Hoxa9, HOXB13, Hoxb4, Hoxc9, Hoxd10, Hoxd11, Hoxd12, Hoxd13, HRE, IRF1, IRF2, IRF3, IRF4, IRF8, Isl1, ISRE, Jun-AP1, JunB, c-Jun-CRE, JunD, KLF10, KLF14, KLF3, Klf4, KLF5, KLF6, Klf9, LEF1, Lhx1, Lhx2, Lhx3, LXH9, Unknown, LRF, Nr5a2, LXRE, MafA, MafB, MafF, MafK, Max, Maz, Mef2a, Mef2b, Mef2c, Mef2d, Meis1, MITF, MNT, AMYB, MYB, Myf5, MYNN, MyoD, MyoG, Nanog, NeuroD1, NeuroG2, NF1:FOXA1, NF1-halfsite, NF1, NFAT:AP1, NFAT, NFE2L2, NF-E2, NFIL3, NFY, Nkx2.1, Nkx2.2, Nkx2.5, Nkx3.1, Nkx6.1, n-Myc, NPAS2, NPAS, EAR2, Nrf2, NRF, Nur77, Oct11, Oct2, Oct4, Oct4:Sox17, OCT4-SOX2-TCF-NANOG, Oct6, OCT:OCT, OCT:OCT-short, Olig2, Otx2, NFkB-p50,p52, p53, p63, NFkB-p65, PAX3:FKHR-fusion, PAX5, PAX6, Pax7, Pax8, PBX2, Pbx3, Pdx1, PGR, Phox2a, Phox2b, Pit1+1bp, Pit1, Pitx1:Ebox, Pitx1, Pknox1, PPARa, PPARE, PRDM10, PRDM14, PRDM15, PRDM1, PRDM9, PR, Prop1, PSE, Ptf1a, PU.1:IRF8, PU.1-IRF, PU.1, RARa, RBPJ:Ebox, Rbpj1, Reverb, Rfx1, Rfx2, RFX, Rfx5, Rfx6, RORa, RORgt, RORg, RUNX, RUNX1, RUNX2, RUNX-AML, RXR, SCL, SCRT1, SF1, Six1, Six2, Six4, Smad2, Smad3, Smad4, Slug, Sox10, Sox15, Sox17, Sox2, Sox3, Sox4, Sox6, Sox9, Sp1, Sp2, Sp5, SPDEF, SpiB, Srebp1a, Srebp2, CArG, ZNF143|STAF, STAT1, Stat3+il21, Stat3, STAT4, STAT5, STAT6, TATA-Box, Tbet, Tbox:Smad, Tbr1, Tbx20, Tbx21, Tbx5, Tbx6, Tcf12, Tcf21, Tcf3, TCF4, TCFL2, Tcf7, Tcfcp2l1, TEAD1, TEAD2, TEAD3, TEAD4, TEAD, Tgif1, Tgif2, THRb, Tlx?, TR4, THRa, TRPS1, Twist, Twist2, USF1, Usf2, VDR, WT1, X-box, YY1, Zac1, ZBTB12, ZBTB18, ZBTB33, ZEB1, ZEB2, Zfp281, Zfp57, Zfp809, ZFX, Zic3, Zic, ZKSCAN1, ZNF136, ZNF165, ZNF189, Znf263, ZNF264, ZNF322, ZNF341, ZNF415, ZNF416, ZNF467, ZNF519, ZNF652/HepG2-ZNF652.Flag-ChIP-Seq, ZNF669, ZNF675, ZNF692, ZNF711, ZNF768, ZNF7, ZSCAN22

unique_SNPs <- unique(sub("*\\-.*", "", quads_interesting_analysis$'PositionID')) # 254 unique SNPs
length(unique_SNPs) # 254 unique SNPs

# glue_collapse(unique_SNPs, sep=", ")
# chr22:32003939, rs1000042, rs10402909, rs10163743, rs1003904, rs1011009, rs1024467, rs1000778, rs1008727, rs10086941, rs10085974, rs10086660, rs1000215, rs10238746, rs10046370, rs1000287, chr3:46757116, rs10048736, rs10084408, chr1:161483656, chr1:25617251, chr22:42322087, rs1005455, chr19:47104678, chr19:7755472, rs10412640, rs1009570, rs117163045, chr16:74455003, rs1053293, chr16:347928, chr15:79066932, rs1007381, rs1001852, rs1045546, rs10160257, chr11:60610065, rs1039839, rs10219104, rs10117795, rs1008278, rs1010790, chr6:26368938, chr3:52552916, rs10510667, rs10167565, rs10012, rs10157091, rs10157139, rs1040010, rs1001511, rs114090926, rs1005536, rs111346628, rs1004439, chr16:632225, rs10444445, rs1018696, rs10114049, rs10114028, chr8:145541424, rs1003696, rs1003000, rs10022956, rs10008770, rs1004588, rs10049226, rs10159433, rs1005770, rs10493007, rs1011476, chr7:66660247, rs10000905, rs1032465, rs1002408, rs10084634, rs1003680, rs1000728, rs1004017, rs1016754, rs1006938, rs1001609, rs10153482, chr19:15918647, rs1007352, chr19:6935113, chr19:4511278, rs10406751, rs10163741, rs1015760, rs1018242, rs1003394, rs1005207, rs1004770, rs1018190, chr16:1551629, rs10128872, rs10160844, rs10459059, rs10161486, rs10082568, rs10212, rs10159713, chr10:71168653, rs1000875, chr10:5005761, rs1002095, rs1007394, rs10086255, rs10086568, rs1005335, rs1019227, rs1012371, rs1007626, rs1016822, rs10080924, rs1003973, rs1075846, rs10036474, rs10000173, rs10000509, rs1002410, rs10165492, rs10153633, rs10166707, rs1002297, chr1:212142055, rs1006310, rs10158415, rs1008449, rs1001207, rs10442705, rs59319728, rs10047805, rs10113984, rs10084748, rs10164896, rs1001513, rs1005413, rs1003580, rs108816, rs10403773, rs10221279, rs10747696, rs10171775, chr1:149813866, rs10048276, rs10130401, rs10086976, rs10080355, chr12:10763175, rs10511920, rs1002935, rs10035211, rs10172772, rs10158888, rs1016316, rs1001425, chr16:71956505, rs112551196, rs10068, rs10047970, rs1002579, chr12:123186880, rs10450723, rs10492289, rs1000933, rs1005559, rs10047657, chr12:9083164, chr12:6909388, rs10047467, chr11:60165353, rs10501381, rs1001956, rs116324892, rs10039597, chr5:102518916, rs10454291, rs10041573, rs10043273, rs10000266, rs10000034, rs10000433, rs10045, rs10049413, rs10048628, rs11240666, rs10157502, rs10081982, rs1029083, rs10153412, rs1011339, rs1005949, rs10224583, rs1014444, rs10049045, rs1026916, rs10114702, rs10095753, rs10046529, chr6:31783863, rs1003649, chr2:71047649, rs111742659, rs10160019, rs10036532, rs1003855, rs1007992, rs111863112, rs1012621, rs10160732, rs10114108, rs111228463, rs10736397, chr12:8693418, rs1020364, rs12440789, rs1004028, rs10000553, rs1000107, rs1003348, rs1018051, rs1034995, rs1030826, rs1001362, rs10456, rs1008562, rs10158440, rs10449322, rs10132100, chr16:66968283, rs10035526, rs10036258, rs1007923, rs1003431, rs1029300, rs10115374, rs10432735, rs1016140, chr14:101004844, rs1008137, rs1004467, rs10159070, rs1005752, chr12:52842736, rs1000768, rs10047326, rs1009227, rs10081672, chr11:63719761, rs1029474, rs1012912, rs10035677

# Rename Homer Output PositionID to SNP2, Abbreviate column
colnames(quads_interesting_analysis)[1] <- "SNP2"
quads_interesting_analysis$SNP2_Abbreviated <- sub("*\\-.*", "", quads_interesting_analysis$'SNP2')

# Select and reorganize columns from original Z file
zz <- z[,c("SNP1", "cisGene", "SNP2", "transGene", "cisZscore", "cisSNPChr", "cisSNPPos", "transZscore", "transSNPChr", "transSNPPos", "Zscore", "SNPChr", "SNPPos", "cisGene_commonname")] # 646241 x 14
colnames(zz) <- c("SNP1", "cisGene", "SNP2", "transGene", "cisZscore", "cisSNPChr", "cisSNPPos", "transZscore", "transSNPChr", "transSNPPos", "Zscore", "SNPChr", "SNPPos", "cisGene_commonname")

zzz <- merge(quads_interesting_analysis, zz, by.x = "SNP2_Abbreviated", by.y = "SNP2") # Why does 10562 x 6, 646241 x 13...yield 34933 x 20?
# Verify: length(zzz[!duplicated(zzz$SNP2),"SNP2"]) # 257 unique SNPs from joining, not the same as 254 unique SNPs from unique_SNPs

zzz$SNP2 <- sub("*\\-.*", "", zzz$'SNP2')
# Verify: length(zzz[!duplicated(zzz$SNP2),"SNP2"]) # 254 unique SNPs from joining, now same as 254 unique SNPs from unique_SNPs
zzzz <- zzz[!is.na(zzz$SNP2),] # 34933 x 20

# don't worry about this
# # zzz <- merge(quads_interesting_analysis, zz, by.x = "SNP2_Abbreviated", by.y = "SNP2", all.x = T, all.y=F) # Why does 10562 x 6, 646241 x 13...yield 34933 x 19?   
# # Verify: zzz[!duplicated(zzz$SNP2),"SNP2"] # 257 unique SNPs from joining
# # #zzz chr16:347928-Dup1, rs1001425-Dup1, rs1001425-Dup2
# 
# # Match information from homer output and original Z file into original zz file
# # zzz <- zz[match(quads_interesting_analysis$SNP2_Abbreviated, zz$SNP2)] # 10562 x 13
# # zzzz <- zzz[!is.na(zzz$SNP1)] # 10562 x 13
# 
# # zzz <- match(quads_interesting_analysis$SNP2_Abbreviated, zz$SNP2) # 10562 x 13 # can't make it into quads_interesting_analysis
# # zzzz <- zzz[!is.na(zzz$SNP1)] # 10562 x 13

# Match transgene to cisgene
zzzz$MotifName_Abbreviated <- toupper(sub("*\\(.*", "", zzzz$'Motif.Name'))

doesItMatch <- function(cisGene, TFMotif){ 
  if (cisGene==TFMotif)
    return(1)
  else
    return(0)
}

zzzz$DoesItMatch <- 0
for(i in 1:nrow(zzzz)){
  zzzz$DoesItMatch[i] <- doesItMatch(toupper(zzzz$cisGene_commonname[i]), toupper(zzzz$MotifName_Abbreviated[i]))
}

# View(data.table(zzzz$cisGene_commonname,zzzz$MotifName_Abbreviated,zzzz$DoesItMatch))

sum(zzzz$DoesItMatch) # 121 yes
# length(zzzz$DoesItMatch) - sum(zzzz$DoesItMatch) = 34812
# zzzz[which(zzzz$DoesItMatch == 1),] or zzzz[(zzzz$DoesItMatch == 1),]

# Not sure why this method doesn't work
# # m <- match(toupper(zzzz$cisGene_commonname), toupper(TFs))
# # length(which(is.na(m))) #0
# # length(which(!is.na(m))) #349333
# View(data.table(toupper(zzzz$cisGene_commonname), toupper(TFs), zzzz$DoesItMatch, match(toupper(zzzz$cisGene_commonname), toupper(TFs))))
# # 
# # w <- which(!is.na(m))
# # snps_interesting <- snps[w,]
# # common_cisgenename_interesting <- common_cisgenename[w]
# # sanity check
# # z <- match(toupper(common_cisgenename_interesting), toupper(unique(TFs)))
# # length(which(is.na(z))) == 0
 
# #########
# #motif_matches <- sapply(1:nrow(snps), function(x) paste(TFs[which(dat$Ensembl == as.character(snps[x,3]))], collapse = ","))
# # matches <- c()
# # for (i in 1:nrow(snps_interesting)){
# #   matches[i] <- match(toupper(common_cisgenename_interesting[i]), toupper(TFs[which(dat$Ensembl == as.character(snps_interesting[i,"transGene"]))]))
# # }
# # length(which(!is.na(matches))) #458, NOW 812
# 
# #snps_interesting_TFmatch <- cbind(snps_interesting, ifelse(is.na(matches), 0, 1))
# #w <- which(snps_interesting_TFmatch$V2 == 1)
# snps_interesting_TFmatch <- snps_interesting[which(!is.na(matches)),]
# snps_interesting_TFmatch <- cbind(snps_interesting_TFmatch, common_cisgenename_interesting[which(!is.na(matches))])
# snps_interesting_TFmatch <- snps_interesting_TFmatch[,c(1,2,4,3)]

# 121 Matches
write.table(zzzz[(zzzz$DoesItMatch == 1),], "hypothesis6_output/interesting_information.txt", row.names = F, col.names = T, sep="\t", quote=F)

read.csv("hypothesis6/allgenes.txt")