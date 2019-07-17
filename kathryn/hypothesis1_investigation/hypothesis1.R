HiC_data_2 <- fread(paste0(my_path, "hypothesis1_input/chicp_javierre/PCHiC_peak_matrix_cutoff5.tsv"), sep="\t", header=T) #728,838 x 30
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525

# have you tried using awk to select only the super significant interactions? (I.e. Rows with p<X)

allgenes <- fread(paste0(my_path, "hypothesis1_input/allgenes.txt"), sep="\t", header=F)
colnames(allgenes) <- c("GeneChr", "GeneStart", "GeneEnd", "GeneSymbol", "0", "Strand")

allgenes <- allgenes[,1:4]
allgenes <- allgenes[!duplicated(allgenes),]
allgenes$Size <- allgenes$GeneEnd - allgenes$GeneStart
allgenes <- allgenes %>% group_by(GeneSymbol) %>% filter(Size ==max(Size))

franke_trans_data_copy <- franke_trans_data
franke_trans_data_copy$rank <- 1:nrow(franke_trans_data_copy)
franke_trans_data <- merge(franke_trans_data_copy, allgenes[,c("GeneStart", "GeneEnd", "GeneSymbol")], by="GeneSymbol") # lose 59786 - 54178 = 5608 genes lost
franke_trans_data$rank <- 1:nrow(franke_trans_data)

y <- HiC_data_2 %>%
  filter(franke_trans_data$SNPChr[1] == baitChr) %>%
  filter(franke_trans_data$SNPPos[1] >= baitStart) %>%
  filter(franke_trans_data$SNPPos[1] <= baitEnd)

# if (dim(y)[1] != 0){
#   z <- cbind(franke_trans_data[["rank"]], y[, c("baitChr", "baitStart", "baitEnd", "baitID", "baitName", "oeChr", "oeStart", "oeEnd", "oeID", "oeName")], franke_trans_data[["SNP"]], franke_trans_data[["SNPChr"]], franke_trans_data[["SNPPos"]], franke_trans_data[["Gene"]], franke_trans_data[["GeneSymbol"]], franke_trans_data[["GeneChr"]], franke_trans_data[["GeneStart"]], franke_trans_data[["GeneEnd"]])
#   colnames(z) <- c("row_eQTLGen", "baitChr", "baitStart", "baitEnd", "baitID", "baitName", "oeChr_poss", "oeStart_poss", "oeEnd_poss", "oeID_poss", "oeName_poss", "SNP", "SNPChr", "SNPPos", "Gene", "GeneSymbol", "GeneChr", "GeneStart", "GeneEnd")
#   zzz <- z %>% filter(GeneChr == oeChr_poss) %>% filter(GeneStart >= oeStart_poss | GeneEnd <= oeEnd_poss)
#   
#   if (dim(z)[1] != 0) {
#     matches <- zzz
#     } else 
#     {
#       matches <- -1
#     }
#   } else {
#     matches <- -2
#   }

  if (dim(y)[1] == 0) {
    matches <- data.table(matrix(-2, 1, 19))
  } else {
    z <- cbind(franke_trans_data$rank[1], 
               y[, c("baitChr", "baitStart", "baitEnd", "baitID", "baitName", "oeChr", "oeStart", "oeEnd", "oeID", "oeName")], 
               franke_trans_data[1,c("SNP", "SNPChr", "SNPPos", "Gene", "GeneSymbol", "GeneChr", "GeneStart", "GeneEnd")]
    )
    #[["SNP"]], franke_trans_data[["SNPChr"]], franke_trans_data[["SNPPos"]], franke_trans_data[["Gene"]], franke_trans_data[["GeneSymbol"]], franke_trans_data[["GeneChr"]], franke_trans_data[["GeneStart"]], franke_trans_data[["GeneEnd"]])
    colnames(z) <- c("row_eQTLGen", "baitChr", "baitStart", "baitEnd", "baitID", "baitName", "oeChr_poss", "oeStart_poss", "oeEnd_poss", "oeID_poss", "oeName_poss", "SNP", "SNPChr", "SNPPos", "Gene", "GeneSymbol", "GeneChr", "GeneStart", "GeneEnd")
    zzz <- z %>% filter(GeneChr == oeChr_poss) %>% filter(GeneStart >= oeStart_poss | GeneEnd <= oeEnd_poss)
    if (dim(z)[1] == 0) {
      matches <- data.table(matrix(-1, 1, 19))
    } else {
      matches <- zzz
    }
  }

for(i in 2:10){#nrow(franke_trans_data)){
  y <- HiC_data_2 %>%
    filter(franke_trans_data$SNPChr[i] == baitChr) %>%
    filter(franke_trans_data$SNPPos[i] >= baitStart) %>%
    filter(franke_trans_data$SNPPos[i] <= baitEnd)
  
  if (dim(y)[1] != 0){
    z <- cbind(franke_trans_data$rank[i], 
               y[, c("baitChr", "baitStart", "baitEnd", "baitID", "baitName", "oeChr", "oeStart", "oeEnd", "oeID", "oeName")], 
               franke_trans_data[i,c("SNP", "SNPChr", "SNPPos", "Gene", "GeneSymbol", "GeneChr", "GeneStart", "GeneEnd")]
               )
               #[["SNP"]], franke_trans_data[["SNPChr"]], franke_trans_data[["SNPPos"]], franke_trans_data[["Gene"]], franke_trans_data[["GeneSymbol"]], franke_trans_data[["GeneChr"]], franke_trans_data[["GeneStart"]], franke_trans_data[["GeneEnd"]])
    colnames(z) <- c("row_eQTLGen", "baitChr", "baitStart", "baitEnd", "baitID", "baitName", "oeChr_poss", "oeStart_poss", "oeEnd_poss", "oeID_poss", "oeName_poss", "SNP", "SNPChr", "SNPPos", "Gene", "GeneSymbol", "GeneChr", "GeneStart", "GeneEnd")
    zzz <- z %>% filter(GeneChr == oeChr_poss) %>% filter(GeneStart >= oeStart_poss | GeneEnd <= oeEnd_poss)
    if (dim(z)[1] != 0)
      matches <- rbind(matches, zzz, fill=T)
    else
      matches <- rbind(matches, data.table(matrix(-1, 1, 19)), fill=T)
  } else matches <- rbind(matches, data.table(matrix(-2, 1, 19)), fill=T)
}