#entire STRING database
setwd("C:/Users/Tiffany Amariuta/Documents/SRLAB/STRINGDATA") 
ensp <- read.table("ENSG_ENSP_conversion.txt", sep="\t", header = T, stringsAsFactors = F)
#protein_links <- read.table("9606.protein.links.v10.txt", sep="\t", header = T, stringsAsFactors = F)
#protein_actions <- read.table("9606.protein.actions.v10.txt", sep="\t", header = T, stringsAsFactors = F)
ultmatrix <- read.table("ultmatrix_2kB_nodecimal.txt", sep="\t", header = T, stringsAsFactors = F)
#ultmatrix <- read.table("TPM_QC_fulltimecourse_bulk_nornalabel.txt", sep="\t", header = T, stringsAsFactors = F)
proteinlinks <- read.table("proteinlinks.txt", sep="\t", header = T, stringsAsFactors = F)

#6029/6152 ENSG IDs converted to ENSP IDs 

#### done one time
# proteinlinksscorematrix <- matrix(0, nrow=nrow(protein_links), ncol = 3)
# for (j in 1:nrow(protein_links)){
#   score_vec <- unlist(strsplit(protein_links$protein1.protein2.combined_score[j], split = " "))
#   proteinlinksscorematrix[j,1] <- strsplit(score_vec[1], "[.]")[[1]][2]
#   proteinlinksscorematrix[j,2] <- strsplit(score_vec[2], "[.]")[[1]][2]
#   proteinlinksscorematrix[j,3] <- as.numeric(score_vec[3])
# }
# write.table(proteinlinksscorematrix, file = "proteinlinks.txt", sep="\t")
####

#SET THIS#
genes <- 1000 
##########
residual_sums_weighting <- matrix(0, nrow = 1, ncol = 6)
residual_info_weighting <- matrix(0, nrow = 3, ncol = 6)

for (l in 1:6){
  assign(paste0("feature_matrix_",l), matrix(0, genes, 2)) 
  assign(paste0("sum_",l), matrix(0, genes, 1))
  assign(paste0("golden_matrix_",l), matrix(0, genes, 1))
}
#Ax = b; A is padded feature, x is what we solve for, b is golden
feature <- list(feature_matrix_1,feature_matrix_2,feature_matrix_3,feature_matrix_4,feature_matrix_5,feature_matrix_6)
sum <- list(sum_1, sum_2, sum_3, sum_4, sum_5, sum_6)
golden <- list(golden_matrix_1, golden_matrix_2, golden_matrix_3, golden_matrix_4, golden_matrix_5, golden_matrix_6)

transcript_1 <- ultmatrix$X0hr
transcript_2 <- ultmatrix$X4hr
transcript_3 <- ultmatrix$X8hr
transcript_4 <- ultmatrix$X12hr
transcript_5 <- ultmatrix$X24hr
transcript_6 <- ultmatrix$X72hr

for (w in 1:6){
  feature[[w]][,1] <- 1 
}

#note: the STRING pairing is not comprehensive, not every single protein interaction is modeled  
for (i in 1:genes){
  for (z in 1:6){golden[[z]][i] <- get(paste0("transcript_", z))[i]}
  
  gene_name <- ultmatrix$Gene_Name[i]  
  gene_ID <- ultmatrix$Gene_ID[i]   
  protein_matches <- ensp$Ensembl.Protein.ID[which(gene_ID == ensp$Ensembl.Gene.ID)]
  for (n in 1:length(protein_matches)){
    protein_match <- protein_matches[n]
    interacting_proteins_indicies <- which(protein_match == proteinlinks[,1])
    if (length(interacting_proteins_indicies) == 0){
    }
    else{
      interacting_proteins <- proteinlinks$V2[interacting_proteins_indicies]
      interacting_proteins_scores <- proteinlinks[interacting_proteins_indicies,3]
      for (p in 1:length(interacting_proteins)){
        m <- match(interacting_proteins[p], ensp$Ensembl.Protein.ID, incomparables = NULL)
        if (is.na(m)){
        }
        else{
          matching_interacting_gene <- ensp$Ensembl.Gene.ID[m]
          for (t in 1:6){
          matching_gene_transcript <- get(paste0("transcript_", t))[match(matching_interacting_gene, ultmatrix$Gene_ID)]
          feature[[t]][i,2] <- feature[[t]][i,2] + matching_gene_transcript*interacting_proteins_scores[p]
          sum[[t]][i] <- sum[[t]][i] + interacting_proteins_scores[p]
          }
        }
      }
    }
  }
  for (h in 1:6){
    if (sum[[h]][i] == 0){
      feature[[h]][i,2] <- 0 
    } 
    else{
      feature[[h]][i,2] <- feature[[h]][i,2]/sum[[h]][i]
    }
  }
}

#at this point all the real transcript values are filled in -> golden 
#and all the weighted averages based on interactions have been calculated -> feature

h <- 1 
b <- matrix(unlist(golden[h]))   #real transcript from bulk 
A <- matrix(unlist(feature[h]), ncol = 2, byrow = FALSE)   #weighted transcript from interactions
lsfit <- lsfit(A, b, intercept = FALSE) #lsfit$residuals has all the errors 
residual_sums_weighting[1,h] <- sum(abs(lsfit$residuals))/genes
residual_info_weighting[1,h] <- var(lsfit$residuals)
residual_info_weighting[2,h] <- max(lsfit$residuals)
residual_info_weighting[3,h] <- min(abs(lsfit$residuals))

#    }
#  }
#}


#### RESIDUAL PLOT ####

#this is a higher error than for using the tsv file in STRING_imputation_weights_padding
#residual_total/10   #1.82 average error for 10 genes at 0 hr (w/o lsfit)

x <- c(1,2,3,4,5,6)
y <- c(1.156335, 1.359033, 1.297591, 1.423481, 1.29539, 1.162241)  #10 genes
z <- c(1.487769, 1.822165, 1.876369, 1.912533, 1.884362, 1.481904) #50 genes
w <- c(1.461369, 1.703524, 1.705558, 1.722738, 1.786906, 1.550239) #100 genes
v <- c(1.541727, 1.85663, 1.790134, 1.734982, 1.671338, 1.583771)  #500 genes
u <- c(1.580951, 1.834349, 1.784913, 1.733783, 1.68408, 1.591917)  #1000 genes

pdf("residual_plots_wholeDB.pdf")
plot(x,y, col = 5, bg = 5, xlab = "time points", ylab = "average residual error", main = "RNA Seq imputation using STRING interactions", ylim = c(min(v,w,y,z), max(v,w,y,z) + 0.3))
lines(x,y, col = 5, bg = 5)
points(x,z, col = 4, bg = 4)
lines(x,z, col = 4, bg = 4)
points(x,w, col = 3, bg = 3)
lines(x,w, col = 3, bg = 3)
points(x,v, col = 2, bg = 2)
lines(x,v, col = 2, bg = 2)
points(x,u, col = 6, bg = 6)
lines(x,u, col = 6, bg = 6)
legend("topleft", c("10 genes", "50 genes", "100 genes", "500 genes", "1000 genes"), fill = c(5,4,3,2,6))
dev.off()


#### VARIANCE PLOT ####
x <- c(1,2,3,4,5,6)
y <- c(2.385655e+00, 2.812627e+00, 2.247302e+00, 2.835941e+00, 2.857566e+00, 2.254396e+00)  #10 genes
z <- c(3.485763e+00, 4.567034e+00, 4.693235e+00, 4.924042e+00, 5.014068e+00, 3.240128e+00) #50 genes
w <- c(3.388739e+00, 4.218965e+00, 4.140600e+00, 4.268371e+00, 4.565712e+00, 3.653413e+00) #100 genes
v <- c(3.847550e+00, 5.118678e+00, 4.910070e+00, 4.700169e+00, 4.469773e+00, 4.089306e+00)  #500 genes
#u <- c()  #1000 genes

pdf("varofresidual_plots_wholeDB.pdf")
plot(x,y, col = 5, bg = 5, xlab = "time points", ylab = "average residual error", main = "RNA Seq imputation using STRING interactions", ylim = c(min(v,w,y,z), max(v,w,y,z) + 0.3))
lines(x,y, col = 5, bg = 5)
points(x,z, col = 4, bg = 4)
lines(x,z, col = 4, bg = 4)
points(x,w, col = 3, bg = 3)
lines(x,w, col = 3, bg = 3)
points(x,v, col = 2, bg = 2)
lines(x,v, col = 2, bg = 2)
#points(x,u, col = 6, bg = 6)
#lines(x,u, col = 6, bg = 6)
legend("topleft", c("10 genes", "50 genes", "100 genes", "500 genes"), fill = c(5,4,3,2))
dev.off()
