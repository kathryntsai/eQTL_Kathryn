
features <- read.table("q4_output/features_forKathryn.txt", header = F, stringsAsFactors = F, sep = " ")
features <- features$V1

mycoef <- coef(ENet_fit)
mycoef1 <- mycoef[,1]
pdf(paste0("q4_output/beta_plot_cisvtrans_eqtl.pdf"), height = 5)
plotBetas_pretty(mycoef1, "cis v trans plot")
dev.off()

plotBetas_pretty <- function(mean_b, thisname){
  betas <- mean_b[-1] #remove intercept 
  betas[abs(betas) > 6] <- 6
  #betas <- betas/min(betas)
  #betas <- rowSums(mean_b)/ncol(mean_b)
  mean_beta <- mean(betas) 
  sd_beta <- sd(betas) 
  abs_betas <- abs(betas)
  w <- which(abs_betas != 0) #(mean_beta+1*sd_beta)) #control stringency for plotting #mean_beta+1*sd_beta
  abs_betas_sorted <- sort(abs_betas, decreasing = T)
  #at most take 10 w 
  cutoff <- min(length(w), 15)
  abs_betas_sorted <- abs_betas_sorted[1:cutoff]
  w <- match(abs_betas_sorted, abs_betas) #w update
  
  allfeatures_thresh <- features[w]
  betas_thresh <- betas[w]
  
  #significance 
  #distance from mean beta value
  num_standarddevs <- floor(abs(mean_beta - betas_thresh)/sd_beta)
  enrichment_stars <- num_standarddevs #dummy vector
  #star system -> "" = 1 st dev away, * = 2, ** = 3, *** >= 4
  for (i in 1:length(num_standarddevs)){
    if (num_standarddevs[i] <= 1){enrichment_stars[i] <- ""} #btwn 1 and 2
    if (num_standarddevs[i] == 2){enrichment_stars[i] <- "*"} #btwn 2 and 3
    if (num_standarddevs[i] == 3){enrichment_stars[i] <- "**"} #btwn 3 and 4 
    if (num_standarddevs[i] >= 4){enrichment_stars[i] <- "***"} #btwn 4 and above
  }
  beta_matrix_df <- data.frame(allfeatures_thresh, betas_thresh, enrichment_stars)
  beta_matrix_df$betas_thresh <- as.numeric(as.character(beta_matrix_df$betas_thresh))
  beta_matrix_df <- beta_matrix_df[order(beta_matrix_df$betas_thresh),]
  
  
  
  w_posbeta <- which(betas_thresh > 0)
  w_negbeta <- which(betas_thresh <= 0)
  plotbetas <- as.numeric(beta_matrix_df$betas_thresh)
  
  if (length(w_negbeta) > 0 & length(w_posbeta) > 0){ #mix of both 
    mneg <- max(which(beta_matrix_df$betas_thresh < 0))
    plotbetas <- c(plotbetas[1:mneg],0,plotbetas[(mneg+1):length(plotbetas)])
  }
  #add hash here
  allfeat <- sub("_", " ", beta_matrix_df$allfeatures_thresh)
  allfeat <- sub("_", " ", allfeat)
  allfeat <- sub("_", " ", allfeat)
  allfeat <- sub("_", " ", allfeat)
  labelsize <- 0.8
  
  if (length(w_negbeta) > 0 & length(w_posbeta) == 0){ #only negative betas 
    y <- barplot(plotbetas,horiz = T,cex.axis = 1, xlim= c(min(plotbetas)-1, max(plotbetas)+2),
                 cex.lab = 1, cex.main = 1, yaxt = "n", xlab = "Feature Coefficient", main = paste0("IMPACT fit ",thisname), 
                 col = c(rep("#F8766D",length(w_negbeta)+1), rep("#00BFC4", length(w_posbeta))))
    text(y = y[1:length(w_negbeta),1], x = 0.05, col = "black", allfeat[1:length(w_negbeta)], srt = 0, cex = labelsize,xpd = T, adj = 0)
    #text(y = y[(length(w_negbeta)+2):nrow(y),1], x = plotbetas[(length(w_negbeta)+2):nrow(y)] + 0.05, col = "black", allfeat[(1+length(w_negbeta)):length(betas_thresh)], 
    #     srt = 0, cex = labelsize, xpd = T, adj = 0)
    #segments(x0 = -0.05, y0 = y[which(plotbetas==0)]-0.1, x1 = 0.05, y1 = y[which(plotbetas==0)]+0.3, lwd = 2) #hash
    #segments(x0 = -0.05, y0 = y[which(plotbetas==0)]-0.3, x1 = 0.05, y1 = y[which(plotbetas==0)]+0.1, lwd = 2) #hash
    
    
  }
  
  if (length(w_negbeta) > 0 & length(w_posbeta) > 0){ #both negative and positive betas 
    y <- barplot(plotbetas,horiz = T,cex.axis = 1, xlim= c(min(plotbetas)-1, max(plotbetas)+2),
                 cex.lab = 1, cex.main = 1, yaxt = "n", xlab = "Feature Coefficient", main = paste0("IMPACT fit ",thisname), 
                 col = c(rep("#F8766D",length(w_negbeta)+1), rep("#00BFC4", length(w_posbeta))))
    text(y = y[1:length(w_negbeta),1], x = 0.05, col = "black", allfeat[1:length(w_negbeta)], srt = 0, cex = labelsize,xpd = T, adj = 0)
    text(y = y[(length(w_negbeta)+2):nrow(y),1], x = plotbetas[(length(w_negbeta)+2):nrow(y)] + 0.05, col = "black", allfeat[(1+length(w_negbeta)):length(betas_thresh)], 
         srt = 0, cex = labelsize, xpd = T, adj = 0)
    segments(x0 = -0.05, y0 = y[which(plotbetas==0)]-0.1, x1 = 0.05, y1 = y[which(plotbetas==0)]+0.3, lwd = 2) #hash
    segments(x0 = -0.05, y0 = y[which(plotbetas==0)]-0.3, x1 = 0.05, y1 = y[which(plotbetas==0)]+0.1, lwd = 2) #hash
    
    
  }
  
  if (length(w_posbeta) > 0 & length(w_negbeta) == 0){ #only positive betas, no negatives
    y <- barplot(plotbetas,horiz = T,cex.axis = 1, xlim= c(0, max(plotbetas)+1),
                 cex.lab = 1, cex.main = 1, yaxt = "n", xlab = "Feature Coefficient", main = paste0("IMPACT fit ",thisname), 
                 col = rep("#00BFC4", length(w_posbeta)))
    text(y = y[,1], x = plotbetas + 0.05, col = "black", allfeat[1:length(betas_thresh)], 
         srt = 0, cex = labelsize, xpd = T, adj = 0)
    
  }
  #segments(x0 = -0.05, y0 = y[which(plotbetas==0)]-0.1, x1 = 0.05, y1 = y[which(plotbetas==0)]+0.3, lwd = 2) #hash
  #segments(x0 = -0.05, y0 = y[which(plotbetas==0)]-0.3, x1 = 0.05, y1 = y[which(plotbetas==0)]+0.1, lwd = 2) #hash
}