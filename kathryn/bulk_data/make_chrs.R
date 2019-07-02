# ==========================================================
# Load full chr1
# ==========================================================

full_chr1 <- readRDS("q4_input/full_chr1.rds")
full_chr1_mod <- gsub('c\\(|\\)', "", full_chr1)
full_chr1_mod <- gsub(')', "", full_chr1_mod)
loadRDS(full_chr1_mod, "hypothesis1_output/full_chr1.rds")

# x <- 4985013
# div <- seq_len(x)
# paste(div[x %% div == 0], collapse=", ")
# 1, 3, 11, 29, 33, 87, 319, 957, 5209, 15627, 57299, 151061, 171897, 453183, 1661671, 4985013

# ==========================================================
# Load chr2 - chr22
# ==========================================================

# read in chromosomes from hg19 genome, taken from homer/data/genomes/hg19 from q2_3_input folder
chr2 <- fread("q4_input/chr2.fa")
div <- seq_len(nchar(chr2))
chr2_factors <- paste(div[x %% div == 0], collapse=", ")

chr3 <- fread("q4_input/chr3.fa")
div <- seq_len(nchar(chr3))
chr3_factors <- paste(div[x %% div == 0], collapse=", ")

chr4 <- fread("q4_input/chr4.fa")
div <- seq_len(nchar(chr4))
chr4_factors <- paste(div[x %% div == 0], collapse=", ")

chr5 <- fread("q4_input/chr5.fa")
div <- seq_len(nchar(chr5))
chr5_factors <- paste(div[x %% div == 0], collapse=", ")

chr6 <- fread("q4_input/chr6.fa")
div <- seq_len(nchar(chr6))
chr6_factors <- paste(div[x %% div == 0], collapse=", ")

chr7 <- fread("q4_input/chr7.fa")
div <- seq_len(nchar(chr7))
chr7_factors <- paste(div[x %% div == 0], collapse=", ")

chr8 <- fread("q4_input/chr8.fa")
div <- seq_len(nchar(chr8))
chr8_factors <- paste(div[x %% div == 0], collapse=", ")

chr9 <- fread("q4_input/chr9.fa")
div <- seq_len(nchar(chr9))
chr9_factors <- paste(div[x %% div == 0], collapse=", ")

chr10 <- fread("q4_input/chr10.fa")
div <- seq_len(nchar(chr10))
chr10_factors <- paste(div[x %% div == 0], collapse=", ")

chr11 <- fread("q4_input/chr11.fa")
div <- seq_len(nchar(chr11))
chr11_factors <- paste(div[x %% div == 0], collapse=", ")

chr12 <- fread("q4_input/chr12.fa")
div <- seq_len(nchar(chr12))
chr12_factors <- paste(div[x %% div == 0], collapse=", ")

chr13 <- fread("q4_input/chr13.fa")
div <- seq_len(nchar(chr13))
chr13_factors <- paste(div[x %% div == 0], collapse=", ")

chr14 <- fread("q4_input/chr14.fa")
div <- seq_len(nchar(chr14))
chr14_factors <- paste(div[x %% div == 0], collapse=", ")

chr15 <- fread("q4_input/chr15.fa")
div <- seq_len(nchar(chr15))
chr15_factors <- paste(div[x %% div == 0], collapse=", ")

chr16 <- fread("q4_input/chr16.fa")
div <- seq_len(nchar(chr16))
chr16_factors <- paste(div[x %% div == 0], collapse=", ")

chr17 <- fread("q4_input/chr17.fa")
div <- seq_len(nchar(chr17))
chr17_factors <- paste(div[x %% div == 0], collapse=", ")

chr18 <- fread("q4_input/chr18.fa")
div <- seq_len(nchar(chr18))
chr18_factors <- paste(div[x %% div == 0], collapse=", ")

chr19 <- fread("q4_input/chr19.fa")
div <- seq_len(nchar(chr19))
chr19_factors <- paste(div[x %% div == 0], collapse=", ")

chr20 <- fread("q4_input/chr20.fa")
div <- seq_len(nchar(chr20))
chr20_factors <- paste(div[x %% div == 0], collapse=", ")

chr21 <- fread("q4_input/chr21.fa")
div <- seq_len(nchar(chr21))
chr21_factors <- paste(div[x %% div == 0], collapse=", ")

chr22 <- fread("q4_input/chr22.fa")
div <- seq_len(nchar(chr22))
chr22_factors <- paste(div[x %% div == 0], collapse=", ")

# # ==========================================================
# # Create Chr1 Information
# # ==========================================================
# 
# # subset chr1
# chosen_factor=171897 # need 29 arrays of 171897 rows of 50 chars each
# counter=1
# total_rows <- nrow(chr1) # 4985013
# x <- matrix(0, nrow=total_rows/chosen_factor*2-2, ncol=1)
# for (i in 1:nrow(x)){
#   #print(i)
#   if (i %% 2 == 1){
#     x[i,1] <- chosen_factor * counter 
#   }
#   else{
#     x[i,1] <- chosen_factor * counter + 1
#     counter = counter + 1
#   }
# }
# 
# each
chr1 <- gsub(', |\"|\n', "", paste(chr1[1:x[1]],collapse="")) #should be 171897*50=8594850, then nchar(chr1_1) = 8594853
chr2 <- gsub(', |\"|\n', "", paste(chr1[x[2]:x[3]],collapse=""))
chr3 <- gsub(', |\"|\n', "", paste(chr1[x[4]:x[5]],collapse=""))
chr4 <- gsub(', |\"|\n', "", paste(chr1[x[6]:x[7]],collapse=""))
chr5 <- gsub(', |\"|\n', "", paste(chr1[x[8]:x[9]],collapse=""))
chr6 <- gsub(', |\"|\n', "", paste(chr1[x[10]:x[11]],collapse=""))
chr7 <- gsub(', |\"|\n', "", paste(chr1[x[12]:x[13]],collapse=""))
chr8 <- gsub(', |\"|\n', "", paste(chr1[x[14]:x[15]],collapse=""))
chr9 <- gsub(', |\"|\n', "", paste(chr1[x[16]:x[17]],collapse=""))
chr10 <- gsub(', |\"|\n', "", paste(chr1[x[18]:x[19]],collapse=""))
chr11 <- gsub(', |\"|\n', "", paste(chr1[x[20]:x[21]],collapse=""))
chr12 <- gsub(', |\"|\n', "", paste(chr1[x[22]:x[23]],collapse=""))
chr13 <- gsub(', |\"|\n', "", paste(chr1[x[24]:x[25]],collapse=""))
chr14 <- gsub(', |\"|\n', "", paste(chr1[x[26]:x[27]],collapse=""))
chr15 <- gsub(', |\"|\n', "", paste(chr1[x[28]:x[29]],collapse=""))
chr16 <- gsub(', |\"|\n', "", paste(chr1[x[30]:x[31]],collapse=""))
chr17 <- gsub(', |\"|\n', "", paste(chr1[x[32]:x[33]],collapse=""))
chr18 <- gsub(', |\"|\n', "", paste(chr1[x[34]:x[35]],collapse=""))
chr19 <- gsub(', |\"|\n', "", paste(chr1[x[36]:x[37]],collapse=""))
chr20 <- gsub(', |\"|\n', "", paste(chr1[x[38]:x[39]],collapse=""))
chr21 <- gsub(', |\"|\n', "", paste(chr1[x[40]:x[41]],collapse=""))
chr22 <- gsub(', |\"|\n', "", paste(chr1[x[42]:x[43]],collapse=""))

saveRDS(chr1, 'q2_3_output/chr1.rds')  
saveRDS(chr2, 'q2_3_output/chr2.rds')  
saveRDS(chr3, 'q2_3_output/chr3.rds')  
saveRDS(chr4, 'q2_3_output/chr4.rds')  
saveRDS(chr5, 'q2_3_output/chr5.rds')  
saveRDS(chr6, 'q2_3_output/chr6.rds')  
saveRDS(chr7, 'q2_3_output/chr7.rds')  
saveRDS(chr8, 'q2_3_output/chr8.rds')  
saveRDS(chr9, 'q2_3_output/chr9.rds')  
saveRDS(chr10, 'q2_3_output/chr10.rds')  
saveRDS(chr11, 'q2_3_output/chr11.rds')  
saveRDS(chr12, 'q2_3_output/chr12.rds')  
saveRDS(chr13, 'q2_3_output/chr13.rds')  
saveRDS(chr14, 'q2_3_output/chr14.rds')  
saveRDS(chr15, 'q2_3_output/chr15.rds')  
saveRDS(chr16, 'q2_3_output/chr16.rds')  
saveRDS(chr17, 'q2_3_output/chr17.rds')  
saveRDS(chr18, 'q2_3_output/chr18.rds')  
saveRDS(chr19, 'q2_3_output/chr19.rds')  
saveRDS(chr20, 'q2_3_output/chr20.rds')  
saveRDS(chr21, 'q2_3_output/chr21.rds')  
saveRDS(chr22, 'q2_3_output/chr22.rds')  
