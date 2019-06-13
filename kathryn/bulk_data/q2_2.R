# ==========================================================
# QUESTION 2_2: ARE THE SNPS NEAR A BINDING MOTIF?
# ==========================================================

# ==========================================================
# 2_2-A CREATION OF BED FILE WITH EQTLGEN AND PFIZER DATA FOR ANALYSIS IN HOMER
# ==========================================================

franke_cis_bed_file <- data.table(paste0("chr",franke_cis_data$SNPChr),
                                  franke_cis_data$SNPPos - 20,
                                  franke_cis_data$SNPPos + 20,
                                  franke_cis_data$SNP, # is this right?
                                  0,
                                  0)
colnames(franke_cis_bed_file) <- c("c1","c2","c3","c4","c5","c6")

franke_cis_bed_file_small <- franke_cis_bed_file[!duplicated(franke_cis_bed_file$c4)]
# write.table(franke_cis_bed_file_small, file="q2_2_bedfiles/2/franke_cis_bed_file_small.txt", col.names = F, row.names=F, quote=F, sep='\t')

# 10 M sequences - runs for too long
#findMotifsGenome.pl ~/Documents/franke_cis_bed_file.txt hg19 output/ -find ~/homer/data/knownTFs/vertebrates/known.motifs  > franke_output1.txt

# 4 M sequences - still runs for a while
#findMotifsGenome.pl <input> hg19 output/ -find ~/homer/data/knownTFs/vertebrates/known.motifs  > franke_output1.txt

# https://www.google.com/search?q=IMPORT+large+11+gb+file+to+r&oq=IMPORT+large+11+gb+file+to+r&aqs=chrome..69i57j33l4.6046j0j1&sourceid=chrome&ie=UTF-8
# https://stackoverflow.com/questions/15030910/randomly-sample-a-percentage-of-rows-within-a-data-frame
# https://www.google.com/search?ei=uer_XKyQPOjn_Qa3_bXYBA&q=IMPORT+large+11+gb+file+to+r++caret&oq=IMPORT+large+11+gb+file+to+r++caret&gs_l=psy-ab.3..33i22i29i30.477.4130..4170...0.0..0.804.1097.1j1j6-1......0....1..gws-wiz.......0i71j35i302i39.nT1C0kZk4dU

franke_cis_bed_file_smaller1 <- franke_cis_bed_file_small %>% sample_frac(dim(davenport_data)[1] / dim(franke_cis_bed_file_small)[1])
franke_cis_bed_file_smaller2 <- franke_cis_bed_file_small %>% sample_frac(dim(davenport_data)[1] / dim(franke_cis_bed_file_small)[1])
franke_cis_bed_file_smaller3 <- franke_cis_bed_file_small %>% sample_frac(dim(davenport_data)[1] / dim(franke_cis_bed_file_small)[1])
franke_cis_bed_file_smaller4 <- franke_cis_bed_file_small %>% sample_frac(dim(davenport_data)[1] / dim(franke_cis_bed_file_small)[1])
franke_cis_bed_file_smaller5 <- franke_cis_bed_file_small %>% sample_frac(dim(davenport_data)[1] / dim(franke_cis_bed_file_small)[1])

write.table(franke_cis_bed_file_smaller1, file="q2_2_bedfiles/3/franke_cis_bed_file_smaller1.txt", col.names = F, row.names=F, quote=F, sep='\t')
write.table(franke_cis_bed_file_smaller2, file="q2_2_bedfiles/3/franke_cis_bed_file_smaller2.txt", col.names = F, row.names=F, quote=F, sep='\t')
write.table(franke_cis_bed_file_smaller3, file="q2_2_bedfiles/3/franke_cis_bed_file_smaller3.txt", col.names = F, row.names=F, quote=F, sep='\t')
write.table(franke_cis_bed_file_smaller4, file="q2_2_bedfiles/3/franke_cis_bed_file_smaller4.txt", col.names = F, row.names=F, quote=F, sep='\t')
write.table(franke_cis_bed_file_smaller5, file="q2_2_bedfiles/3/franke_cis_bed_file_smaller5.txt", col.names = F, row.names=F, quote=F, sep='\t')

# findMotifsGenome.pl /Users/kathryntsai/OneDrive\ -\ Villanova\ University/College/2018-2019/Summer\ 2019/TFs_eQTLs_Research/RProjects/eQTLPart1_2/q2_2_bedfiles/3/franke_cis_bed_file_smaller1.txt  hg19 output/ -find ~/homer/data/knownTFs/vertebrates/known.motifs > franke_downsampled_output1.txt
# findMotifsGenome.pl /Users/kathryntsai/OneDrive\ -\ Villanova\ University/College/2018-2019/Summer\ 2019/TFs_eQTLs_Research/RProjects/eQTLPart1_2/q2_2_bedfiles/3/franke_cis_bed_file_smaller2.txt  hg19 output/ -find ~/homer/data/knownTFs/vertebrates/known.motifs > franke_downsampled_output2.txt
# findMotifsGenome.pl /Users/kathryntsai/OneDrive\ -\ Villanova\ University/College/2018-2019/Summer\ 2019/TFs_eQTLs_Research/RProjects/eQTLPart1_2/q2_2_bedfiles/3/franke_cis_bed_file_smaller3.txt  hg19 output/ -find ~/homer/data/knownTFs/vertebrates/known.motifs > franke_downsampled_output3.txt
# findMotifsGenome.pl /Users/kathryntsai/OneDrive\ -\ Villanova\ University/College/2018-2019/Summer\ 2019/TFs_eQTLs_Research/RProjects/eQTLPart1_2/q2_2_bedfiles/3/franke_cis_bed_file_smaller4.txt  hg19 output/ -find ~/homer/data/knownTFs/vertebrates/known.motifs > franke_downsampled_output4.txt
# findMotifsGenome.pl /Users/kathryntsai/OneDrive\ -\ Villanova\ University/College/2018-2019/Summer\ 2019/TFs_eQTLs_Research/RProjects/eQTLPart1_2/q2_2_bedfiles/3/franke_cis_bed_file_smaller5.txt  hg19 output/ -find ~/homer/data/knownTFs/vertebrates/known.motifs > franke_downsampled_output5.txt

# ----

davenport_bed_file <- data.table(paste0("chr", davenport_data$SNP_chr), 
                                 davenport_data$SNP_pos - 20, 
                                 davenport_data$SNP_pos + 20,
                                 davenport_data$SNP, # is this right?
                                 0,
                                 0)
colnames(davenport_bed_file) <- c("c1","c2","c3","c4","c5","c6")
write.table(davenport_bed_file, file="q2_2_bedfiles/davenport_bed_file.txt", col.names = F, row.names=F, quote=F, sep='\t')

# findMotifsGenome.pl /Users/kathryntsai/OneDrive\ -\ Villanova\ University/College/2018-2019/Summer\ 2019/TFs_eQTLs_Research/RProjects/eQTLPart1_2/q2_2_bedfiles/davenport_bed_file.txt hg19 output/ -find ~/homer/data/knownTFs/vertebrates/known.motifs -size 40  > davenport_output.txt

# ==========================================================
# 2_2-B CALCULATE PROPORTIONS
# ==========================================================

# New code - does TF bind exactly to SNP?  Comment this out if there are not exact answers
decide <- function(x){
  #offset >= 20 - len(motif) & offset < 20 & strand = "+"
  if(x['Offset'] >= 20 - nchar(x['Sequence']) & x['Offset'] < 20 & x['Strand'] == "+")
    return (1)
  #offset >= 19 - len(motif) & offset < 21 & strand = "-"
  else if (x['Offset'] >= 19 - nchar(x['Sequence']) & x['Offset'] < 21 & x['Strand'] == "-")
    return (1)
  else
    return (0)
}

davenport_output <- fread("q2_2_input/homer/davenport_output.txt")

# New code
davenport_output$DoesItBindExactly <- apply(davenport_output, 1, decide)
davenport_output <- davenport_output %>% filter(DoesItBindExactly == 1)

# analyze unique because it's faster
davenport_output_unique <- data.frame(unique(davenport_output$PositionID))
davenport_output_unique <- as.data.frame(lapply(davenport_output_unique, FUN = function(x) (gsub("\\-.*$", "", x))))
davenport_output_unique <- unique(davenport_output_unique)
davenport_bed_file_total <- data.table(unique(davenport_bed_file$c4))
# FOR PFIZER DATA:
dim(davenport_output_unique)[1] / dim(davenport_bed_file_total)[1] # 4181 / 4818 * 100 = 92.39779% !! // Exact Binding: 3 / 4525 = 0.06629%

# ----

franke_downsampled_output1 <- fread("q2_2_input/homer/3/franke_downsampled_output1.txt") # 31938
franke_downsampled_output2 <- fread("q2_2_input/homer/3/franke_downsampled_output2.txt") # 30998
franke_downsampled_output3 <- fread("q2_2_input/homer/3/franke_downsampled_output3.txt") # 30838
franke_downsampled_output4 <- fread("q2_2_input/homer/3/franke_downsampled_output4.txt") # 31117
franke_downsampled_output5 <- fread("q2_2_input/homer/3/franke_downsampled_output5.txt") # 31058

# ERROR: Homer sampled from -100 to 100, not -20 to 20: need to rerun analysis

# findMotifsGenome.pl /Users/kathryntsai/OneDrive\ -\ Villanova\ University/College/2018-2019/Summer\ 2019/TFs_eQTLs_Research/RProjects/eQTLPart1_2/q2_2_bedfiles/3/franke_cis_bed_file_smaller1.txt  hg19 ~/Users/kathryntsai/OneDrive\ -\ Villanova\ University/College/2018-2019/Summer\ 2019/TFs_eQTLs_Research/RProjects/eQTLPart1_2/input_data_homer/q2_2/3/ -find ~/homer/data/knownTFs/vertebrates/known.motifs -size 40 > franke_downsampled_output1.txt
# findMotifsGenome.pl /Users/kathryntsai/OneDrive\ -\ Villanova\ University/College/2018-2019/Summer\ 2019/TFs_eQTLs_Research/RProjects/eQTLPart1_2/q2_2_bedfiles/3/franke_cis_bed_file_smaller2.txt  hg19 ~/Users/kathryntsai/OneDrive\ -\ Villanova\ University/College/2018-2019/Summer\ 2019/TFs_eQTLs_Research/RProjects/eQTLPart1_2/input_data_homer/q2_2/3/ -find ~/homer/data/knownTFs/vertebrates/known.motifs -size 40 > franke_downsampled_output2.txt
# findMotifsGenome.pl /Users/kathryntsai/OneDrive\ -\ Villanova\ University/College/2018-2019/Summer\ 2019/TFs_eQTLs_Research/RProjects/eQTLPart1_2/q2_2_bedfiles/3/franke_cis_bed_file_smaller3.txt  hg19 ~/Users/kathryntsai/OneDrive\ -\ Villanova\ University/College/2018-2019/Summer\ 2019/TFs_eQTLs_Research/RProjects/eQTLPart1_2/input_data_homer/q2_2/3/ -find ~/homer/data/knownTFs/vertebrates/known.motifs -size 40 > franke_downsampled_output3.txt
# findMotifsGenome.pl /Users/kathryntsai/OneDrive\ -\ Villanova\ University/College/2018-2019/Summer\ 2019/TFs_eQTLs_Research/RProjects/eQTLPart1_2/q2_2_bedfiles/3/franke_cis_bed_file_smaller4.txt  hg19 ~/Users/kathryntsai/OneDrive\ -\ Villanova\ University/College/2018-2019/Summer\ 2019/TFs_eQTLs_Research/RProjects/eQTLPart1_2/input_data_homer/q2_2/3/ -find ~/homer/data/knownTFs/vertebrates/known.motifs -size 40 > franke_downsampled_output4.txt
# findMotifsGenome.pl /Users/kathryntsai/OneDrive\ -\ Villanova\ University/College/2018-2019/Summer\ 2019/TFs_eQTLs_Research/RProjects/eQTLPart1_2/q2_2_bedfiles/3/franke_cis_bed_file_smaller5.txt  hg19 ~/Users/kathryntsai/OneDrive\ -\ Villanova\ University/College/2018-2019/Summer\ 2019/TFs_eQTLs_Research/RProjects/eQTLPart1_2/input_data_homer/q2_2/3/ -find ~/homer/data/knownTFs/vertebrates/known.motifs -size 40 > franke_downsampled_output5.txt

# franke_downsampled_output1 %>% filter(franke_downsampled_output1$Offset>=30 & franke_downsampled_output1$Offset < 70)
# franke_downsampled_output2 %>% filter(franke_downsampled_output1$Offset>=30 & franke_downsampled_output1$Offset < 70)
# franke_downsampled_output3 %>% filter(franke_downsampled_output1$Offset>=30 & franke_downsampled_output1$Offset < 70)
# franke_downsampled_output4 %>% filter(franke_downsampled_output1$Offset>=30 & franke_downsampled_output1$Offset < 70)
# franke_downsampled_output5 %>% filter(franke_downsampled_output1$Offset>=30 & franke_downsampled_output1$Offset < 70)

franke_downsampled_output1$DoesItBindExactly <- apply(franke_downsampled_output1, 1, decide)
franke_downsampled_output2$DoesItBindExactly <- apply(franke_downsampled_output2, 1, decide)
franke_downsampled_output3$DoesItBindExactly <- apply(franke_downsampled_output3, 1, decide)
franke_downsampled_output4$DoesItBindExactly <- apply(franke_downsampled_output4, 1, decide)
franke_downsampled_output5$DoesItBindExactly <- apply(franke_downsampled_output5, 1, decide)

# New code - filter homer output for whether the TF binds exactly.  Comment this out if there are not exact answers
franke_downsampled_output1 <- franke_downsampled_output1 %>% filter(DoesItBindExactly == 1)
franke_downsampled_output2 <- franke_downsampled_output2 %>% filter(DoesItBindExactly == 1)
franke_downsampled_output3 <- franke_downsampled_output3 %>% filter(DoesItBindExactly == 1)
franke_downsampled_output4 <- franke_downsampled_output4 %>% filter(DoesItBindExactly == 1)
franke_downsampled_output5 <- franke_downsampled_output5 %>% filter(DoesItBindExactly == 1)

franke_downsampled_output1_unique <- data.frame(unique(franke_downsampled_output1$PositionID))
franke_downsampled_output1_unique <- as.data.frame(lapply(franke_downsampled_output1_unique, FUN = function(x) (gsub("\\-.*$", "", x))))
franke_downsampled_output1_unique <- unique(franke_downsampled_output1_unique)
franke_cis_bed_file_smaller1_total <- data.table(unique(franke_cis_bed_file_smaller1[,4])) #same for all 5 samples
dim(franke_downsampled_output1_unique)[1] / dim(franke_cis_bed_file_smaller1_total)[1] # 4449/4818 = 92.34122% // Exact Binding: 9 .1867995%

franke_downsampled_output2_unique <- data.frame(unique(franke_downsampled_output2$PositionID))
franke_downsampled_output2_unique <- as.data.frame(lapply(franke_downsampled_output2_unique, FUN = function(x) (gsub("\\-.*$", "", x))))
franke_downsampled_output2_unique <- unique(franke_downsampled_output2_unique)
dim(franke_downsampled_output2_unique)[1] / dim(franke_cis_bed_file_smaller1_total)[1] # 4480 92.98464% // 5 .1037775%

franke_downsampled_output3_unique <- data.frame(unique(franke_downsampled_output3$PositionID))
franke_downsampled_output3_unique <- as.data.frame(lapply(franke_downsampled_output3_unique, FUN = function(x) (gsub("\\-.*$", "", x))))
franke_downsampled_output3_unique <- unique(franke_downsampled_output3_unique)
dim(franke_downsampled_output3_unique)[1] / dim(franke_cis_bed_file_smaller1_total)[1] # 4497 93.33748% // 3 .0622665%

franke_downsampled_output4_unique <- data.frame(unique(franke_downsampled_output4$PositionID))
franke_downsampled_output4_unique <- as.data.frame(lapply(franke_downsampled_output4_unique, FUN = function(x) (gsub("\\-.*$", "", x))))
franke_downsampled_output4_unique <- unique(franke_downsampled_output4_unique)
dim(franke_downsampled_output4_unique)[1] / dim(franke_cis_bed_file_smaller1_total)[1] # 4476 92.90162% // 4 .083022%

franke_downsampled_output5_unique <- data.frame(unique(franke_downsampled_output5$PositionID))
franke_downsampled_output5_unique <- as.data.frame(lapply(franke_downsampled_output5_unique, FUN = function(x) (gsub("\\-.*$", "", x))))
franke_downsampled_output5_unique <- unique(franke_downsampled_output5_unique)
dim(franke_downsampled_output5_unique)[1] / dim(franke_cis_bed_file_smaller1_total)[1] # 4490 93.1922% // 5 .1037775%


# ==========================================================
# 2_2-C CALCULATE PROPORTION FOR TRANS EQTLS
# ==========================================================

franke_trans_bed_file <- data.table(paste0("chr",franke_trans_data$SNPChr), #59786 x 6
                                    franke_trans_data$SNPPos - 20,
                                    franke_trans_data$SNPPos + 20,
                                    franke_trans_data$SNP,
                                    0,
                                    0)
colnames(franke_trans_bed_file) <- c("c1","c2","c3","c4","c5","c6")

franke_trans_bed_file_small <- franke_trans_bed_file[!duplicated(franke_trans_bed_file$c4)] #3853 x 6
franke_trans_bed_file_smaller1 <- franke_trans_bed_file_small %>% sample_frac(.2) #771 x 6
franke_trans_bed_file_smaller2 <- franke_trans_bed_file_small %>% sample_frac(.2)
franke_trans_bed_file_smaller3 <- franke_trans_bed_file_small %>% sample_frac(.2)
franke_trans_bed_file_smaller4 <- franke_trans_bed_file_small %>% sample_frac(.2)
franke_trans_bed_file_smaller5 <- franke_trans_bed_file_small %>% sample_frac(.2)

write.table(franke_trans_bed_file_smaller1, file="q2_2_bedfiles/3/franke_trans_bed_file_smaller1.txt", col.names = F, row.names=F, quote=F, sep='\t')
write.table(franke_trans_bed_file_smaller2, file="q2_2_bedfiles/3/franke_trans_bed_file_smaller2.txt", col.names = F, row.names=F, quote=F, sep='\t')
write.table(franke_trans_bed_file_smaller3, file="q2_2_bedfiles/3/franke_trans_bed_file_smaller3.txt", col.names = F, row.names=F, quote=F, sep='\t')
write.table(franke_trans_bed_file_smaller4, file="q2_2_bedfiles/3/franke_trans_bed_file_smaller4.txt", col.names = F, row.names=F, quote=F, sep='\t')
write.table(franke_trans_bed_file_smaller5, file="q2_2_bedfiles/3/franke_trans_bed_file_smaller5.txt", col.names = F, row.names=F, quote=F, sep='\t')

# findMotifsGenome.pl /Users/kathryntsai/OneDrive\ -\ Villanova\ University/College/2018-2019/Summer\ 2019/TFs_eQTLs_Research/RProjects/eQTLPart1_2/q2_2_bedfiles/3/franke_trans_bed_file_smaller1.txt   hg19 output/ -find ~/homer/data/knownTFs/vertebrates/known.motifs > franke_trans_downsampled_output1.txt
# findMotifsGenome.pl /Users/kathryntsai/OneDrive\ -\ Villanova\ University/College/2018-2019/Summer\ 2019/TFs_eQTLs_Research/RProjects/eQTLPart1_2/q2_2_bedfiles/3/franke_trans_bed_file_smaller2.txt   hg19 output/ -find ~/homer/data/knownTFs/vertebrates/known.motifs > franke_trans_downsampled_output2.txt
# findMotifsGenome.pl /Users/kathryntsai/OneDrive\ -\ Villanova\ University/College/2018-2019/Summer\ 2019/TFs_eQTLs_Research/RProjects/eQTLPart1_2/q2_2_bedfiles/3/franke_trans_bed_file_smaller3.txt   hg19 output/ -find ~/homer/data/knownTFs/vertebrates/known.motifs > franke_trans_downsampled_output3.txt
# findMotifsGenome.pl /Users/kathryntsai/OneDrive\ -\ Villanova\ University/College/2018-2019/Summer\ 2019/TFs_eQTLs_Research/RProjects/eQTLPart1_2/q2_2_bedfiles/3/franke_trans_bed_file_smaller4.txt   hg19 output/ -find ~/homer/data/knownTFs/vertebrates/known.motifs > franke_trans_downsampled_output4.txt
# findMotifsGenome.pl /Users/kathryntsai/OneDrive\ -\ Villanova\ University/College/2018-2019/Summer\ 2019/TFs_eQTLs_Research/RProjects/eQTLPart1_2/q2_2_bedfiles/3/franke_trans_bed_file_smaller5.txt   hg19 output/ -find ~/homer/data/knownTFs/vertebrates/known.motifs > franke_trans_downsampled_output5.txt

# Note: these variables are reused from the cis analysis
franke_downsampled_output1 <- fread("q2_2_input/homer/3/franke_trans_downsampled_output1.txt") # 31347
franke_downsampled_output2 <- fread("q2_2_input/homer/3/franke_trans_downsampled_output2.txt") # 32005
franke_downsampled_output3 <- fread("q2_2_input/homer/3/franke_trans_downsampled_output3.txt") # 31387
franke_downsampled_output4 <- fread("q2_2_input/homer/3/franke_trans_downsampled_output4.txt") # 31717
franke_downsampled_output5 <- fread("q2_2_input/homer/3/franke_trans_downsampled_output5.txt") # 31727

# New code - does TF bind exactly to SNP?  Comment this out if there are not exact answers
franke_downsampled_output1$DoesItBindExactly <- apply(franke_downsampled_output1, 1, decide)
franke_downsampled_output2$DoesItBindExactly <- apply(franke_downsampled_output2, 1, decide)
franke_downsampled_output3$DoesItBindExactly <- apply(franke_downsampled_output3, 1, decide)
franke_downsampled_output4$DoesItBindExactly <- apply(franke_downsampled_output4, 1, decide)
franke_downsampled_output5$DoesItBindExactly <- apply(franke_downsampled_output5, 1, decide)

# New code - filter homer output for whether the TF binds exactly.  Comment this out if there are not exact answers
franke_downsampled_output1 <- franke_downsampled_output1 %>% filter(DoesItBindExactly == 1)
franke_downsampled_output2 <- franke_downsampled_output2 %>% filter(DoesItBindExactly == 1)
franke_downsampled_output3 <- franke_downsampled_output3 %>% filter(DoesItBindExactly == 1)
franke_downsampled_output4 <- franke_downsampled_output4 %>% filter(DoesItBindExactly == 1)
franke_downsampled_output5 <- franke_downsampled_output5 %>% filter(DoesItBindExactly == 1)

franke_downsampled_output1_unique <- data.frame(unique(franke_downsampled_output1$PositionID))
franke_downsampled_output1_unique <- as.data.frame(lapply(franke_downsampled_output1_unique, FUN = function(x) (gsub("\\-.*$", "", x))))
franke_downsampled_output1_unique <- unique(franke_downsampled_output1_unique)
franke_trans_bed_file_smaller1_total <- data.table(franke_trans_bed_file_smaller1[,4]) #same for all 5 samples
dim(franke_downsampled_output1_unique)[1] / dim(franke_trans_bed_file_smaller1_total)[1] # 721/771 93.51492% // 1 0.1297017%

franke_downsampled_output2_unique <- data.frame(unique(franke_downsampled_output2$PositionID))
franke_downsampled_output2_unique <- as.data.frame(lapply(franke_downsampled_output2_unique, FUN = function(x) (gsub("\\-.*$", "", x))))
franke_downsampled_output2_unique <- unique(franke_downsampled_output2_unique)
dim(franke_downsampled_output2_unique)[1] / dim(franke_trans_bed_file_smaller1)[1] # 734/771 94.94163% // 3 0.3891051%

franke_downsampled_output3_unique <- data.frame(unique(franke_downsampled_output3$PositionID))
franke_downsampled_output3_unique <- as.data.frame(lapply(franke_downsampled_output3_unique, FUN = function(x) (gsub("\\-.*$", "", x))))
franke_downsampled_output3_unique <- unique(franke_downsampled_output3_unique)
dim(franke_downsampled_output3_unique)[1] / dim(franke_trans_bed_file_smaller1)[1] # 732/771 94.42283% // 1 0.1297017%

franke_downsampled_output4_unique <- data.frame(unique(franke_downsampled_output4$PositionID))
franke_downsampled_output4_unique <- as.data.frame(lapply(franke_downsampled_output4_unique, FUN = function(x) (gsub("\\-.*$", "", x))))
franke_downsampled_output4_unique <- unique(franke_downsampled_output4_unique)
dim(franke_downsampled_output4_unique)[1] / dim(franke_trans_bed_file_smaller1)[1] #728/771 0.9338521 // 0 0%

franke_downsampled_output5_unique <- data.frame(unique(franke_downsampled_output5$PositionID))
franke_downsampled_output5_unique <- as.data.frame(lapply(franke_downsampled_output5_unique, FUN = function(x) (gsub("\\-.*$", "", x))))
franke_downsampled_output5_unique <- unique(franke_downsampled_output5_unique)
dim(franke_downsampled_output5_unique)[1] / dim(franke_trans_bed_file_smaller1)[1] # 720/771 100% // 0 0%

