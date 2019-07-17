# ==========================================================
# Load full chr1
# ==========================================================

# full_chr1 <- readRDS("q4_input/full_chr1.rds")
# full_chr1_mod <- gsub('c\\(|\\)', "", full_chr1)
# full_chr1_mod <- gsub(')', "", full_chr1_mod)
# saveRDS(full_chr1_mod, "hypothesis1_output/full_chr1.rds")

# x <- 4985013
# div <- seq_len(x)
# paste(div[x %% div == 0], collapse=", ")
# 1, 3, 11, 29, 33, 87, 319, 957, 5209, 15627, 57299, 151061, 171897, 453183, 1661671, 4985013

# ==========================================================
# Load chr2 - chr22
# ==========================================================

# read in chromosomes from hg19 genome, taken from homer/data/genomes/hg19 from q2_3_input folder
chr2 <- fread("q4_input/chr2.fa")
# div <- seq_len(nchar(chr2))
# chr2_factors <- paste(div[x %% div == 0], collapse=", ")

chr3 <- fread("q4_input/chr3.fa")
# div <- seq_len(nchar(chr3))
# chr3_factors <- paste(div[x %% div == 0], collapse=", ")

chr4 <- fread("q4_input/chr4.fa")
# div <- seq_len(nchar(chr4))
# chr4_factors <- paste(div[x %% div == 0], collapse=", ")

chr5 <- fread("q4_input/chr5.fa")
# div <- seq_len(nchar(chr5))
# chr5_factors <- paste(div[x %% div == 0], collapse=", ")

chr6 <- fread("q4_input/chr6.fa")
# div <- seq_len(nchar(chr6))
# chr6_factors <- paste(div[x %% div == 0], collapse=", ")

chr7 <- fread("q4_input/chr7.fa")
# div <- seq_len(nchar(chr7))
# chr7_factors <- paste(div[x %% div == 0], collapse=", ")

chr8 <- fread("q4_input/chr8.fa")
# div <- seq_len(nchar(chr8))
# chr8_factors <- paste(div[x %% div == 0], collapse=", ")

chr9 <- fread("q4_input/chr9.fa")
# div <- seq_len(nchar(chr9))
# chr9_factors <- paste(div[x %% div == 0], collapse=", ")

chr10 <- fread("q4_input/chr10.fa")
# div <- seq_len(nchar(chr10))
# chr10_factors <- paste(div[x %% div == 0], collapse=", ")

chr11 <- fread("q4_input/chr11.fa")
# div <- seq_len(nchar(chr11))
# chr11_factors <- paste(div[x %% div == 0], collapse=", ")

chr12 <- fread("q4_input/chr12.fa")
# div <- seq_len(nchar(chr12))
# chr12_factors <- paste(div[x %% div == 0], collapse=", ")

chr13 <- fread("q4_input/chr13.fa")
# div <- seq_len(nchar(chr13))
# chr13_factors <- paste(div[x %% div == 0], collapse=", ")

chr14 <- fread("q4_input/chr14.fa")
# div <- seq_len(nchar(chr14))
# chr14_factors <- paste(div[x %% div == 0], collapse=", ")

chr15 <- fread("q4_input/chr15.fa")
# div <- seq_len(nchar(chr15))
# chr15_factors <- paste(div[x %% div == 0], collapse=", ")

chr16 <- fread("q4_input/chr16.fa")
# div <- seq_len(nchar(chr16))
# chr16_factors <- paste(div[x %% div == 0], collapse=", ")

chr17 <- fread("q4_input/chr17.fa")
# div <- seq_len(nchar(chr17))
# chr17_factors <- paste(div[x %% div == 0], collapse=", ")

chr18 <- fread("q4_input/chr18.fa")
# div <- seq_len(nchar(chr18))
# chr18_factors <- paste(div[x %% div == 0], collapse=", ")

chr19 <- fread("q4_input/chr19.fa")
# div <- seq_len(nchar(chr19))
# chr19_factors <- paste(div[x %% div == 0], collapse=", ")

chr20 <- fread("q4_input/chr20.fa")
# div <- seq_len(nchar(chr20))
# chr20_factors <- paste(div[x %% div == 0], collapse=", ")

chr21 <- fread("q4_input/chr21.fa")
# div <- seq_len(nchar(chr21))
# chr21_factors <- paste(div[x %% div == 0], collapse=", ")

chr22 <- fread("q4_input/chr22.fa")
# div <- seq_len(nchar(chr22))
# chr22_factors <- paste(div[x %% div == 0], collapse=", ")

# # ==========================================================
# # Create Information Part 1
# # ==========================================================
# 
chr2 <- gsub(', |\"|\n', "", chr2)
full_chr2 <- gsub('c\\(|\\)', "", full_chr2)
full_chr2 <- gsub(')', "", full_chr2)

chr3 <- gsub(', |\"|\n', "", chr3)
full_chr3 <- gsub('c\\(|\\)', "", full_chr3)
full_chr3 <- gsub(')', "", full_chr3)

chr4 <- gsub(', |\"|\n', "", chr4)
full_chr4 <- gsub('c\\(|\\)', "", full_chr4)
full_chr4 <- gsub(')', "", full_chr4)

chr5 <- gsub(', |\"|\n', "", chr5)
full_chr5 <- gsub('c\\(|\\)', "", full_chr5)
full_chr5 <- gsub(')', "", full_chr5)

chr6 <- gsub(', |\"|\n', "", chr6)
full_chr6 <- gsub('c\\(|\\)', "", full_chr6)
full_chr6 <- gsub(')', "", full_chr6)

chr7 <- gsub(', |\"|\n', "", chr7)
full_chr7 <- gsub('c\\(|\\)', "", full_chr7)
full_chr7 <- gsub(')', "", full_chr7)

chr8 <- gsub(', |\"|\n', "", chr8)
full_chr8 <- gsub('c\\(|\\)', "", full_chr8)
full_chr8 <- gsub(')', "", full_chr8)

chr9 <- gsub(', |\"|\n', "", chr9)
full_chr9 <- gsub('c\\(|\\)', "", full_chr9)
full_chr9 <- gsub(')', "", full_chr9)

chr10 <- gsub(', |\"|\n', "", chr10)
full_chr10 <- gsub('c\\(|\\)', "", full_chr10)
full_chr10 <- gsub(')', "", full_chr10)

chr11 <- gsub(', |\"|\n', "", chr11)
full_chr11 <- gsub('c\\(|\\)', "", full_chr11)
full_chr11 <- gsub(')', "", full_chr11)

chr12 <- gsub(', |\"|\n', "", chr12)
full_chr12 <- gsub('c\\(|\\)', "", full_chr12)
full_chr12 <- gsub(')', "", full_chr12)

chr13 <- gsub(', |\"|\n', "", chr13)
full_chr13 <- gsub('c\\(|\\)', "", full_chr13)
full_chr13 <- gsub(')', "", full_chr13)

chr14 <- gsub(', |\"|\n', "", chr14)
full_chr14 <- gsub('c\\(|\\)', "", full_chr14)
full_chr14 <- gsub(')', "", full_chr14)

chr15 <- gsub(', |\"|\n', "", chr15)
full_chr15 <- gsub('c\\(|\\)', "", full_chr15)
full_chr15 <- gsub(')', "", full_chr15)

chr16 <- gsub(', |\"|\n', "", chr16)
full_chr16 <- gsub('c\\(|\\)', "", full_chr16)
full_chr16 <- gsub(')', "", full_chr16)

chr17 <- gsub(', |\"|\n', "", chr17)
full_chr17 <- gsub('c\\(|\\)', "", full_chr17)
full_chr17 <- gsub(')', "", full_chr17)

chr18 <- gsub(', |\"|\n', "", chr18)
full_chr18 <- gsub('c\\(|\\)', "", full_chr18)
full_chr18 <- gsub(')', "", full_chr18)

chr19 <- gsub(', |\"|\n', "", chr19)
full_chr19 <- gsub('c\\(|\\)', "", full_chr19)
full_chr19 <- gsub(')', "", full_chr19)

chr20 <- gsub(', |\"|\n', "", chr20)
full_chr20 <- gsub('c\\(|\\)', "", full_chr20)
full_chr20 <- gsub(')', "", full_chr20)

chr21 <- gsub(', |\"|\n', "", chr21)
full_chr21 <- gsub('c\\(|\\)', "", full_chr21)
full_chr21 <- gsub(')', "", full_chr21)

chr22 <- gsub(', |\"|\n', "", chr22)
full_chr22 <- gsub('c\\(|\\)', "", full_chr22)
full_chr22 <- gsub(')', "", full_chr22)


saveRDS(full_chr1, 'hypothesis1_output/full_chr1.rds')  
saveRDS(full_chr2, 'hypothesis1_output/full_chr2.rds')  
saveRDS(full_chr3, 'hypothesis1_output/full_chr3.rds')  
saveRDS(full_chr4, 'hypothesis1_output/full_chr4.rds')  
saveRDS(full_chr5, 'hypothesis1_output/full_chr5.rds')  
saveRDS(full_chr6, 'hypothesis1_output/full_chr6.rds')  
saveRDS(full_chr7, 'hypothesis1_output/full_chr7.rds')  
saveRDS(full_chr8, 'hypothesis1_output/full_chr8.rds')  
saveRDS(full_chr9, 'hypothesis1_output/full_chr9.rds')  
saveRDS(full_chr10, 'hypothesis1_output/full_chr10.rds')  
saveRDS(full_chr11, 'hypothesis1_output/full_chr11.rds')  
saveRDS(full_chr12, 'hypothesis1_output/full_chr12.rds')  
saveRDS(full_chr13, 'hypothesis1_output/full_chr13.rds')  
saveRDS(full_chr14, 'hypothesis1_output/full_chr14.rds')  
saveRDS(full_chr15, 'hypothesis1_output/full_chr15.rds')  
saveRDS(full_chr16, 'hypothesis1_output/full_chr16.rds')  
saveRDS(full_chr17, 'hypothesis1_output/full_chr17.rds')  
saveRDS(full_chr18, 'hypothesis1_output/full_chr18.rds')  
saveRDS(full_chr19, 'hypothesis1_output/full_chr19.rds')  
saveRDS(full_chr20, 'hypothesis1_output/full_chr20.rds')  
saveRDS(full_chr21, 'hypothesis1_output/full_chr21.rds')  
saveRDS(full_chr22, 'hypothesis1_output/full_chr22.rds')  
