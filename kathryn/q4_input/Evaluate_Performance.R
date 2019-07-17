args = commandArgs(trailingOnly=TRUE)
#List of Functions
library(glmnet)
library(ROCR)
library(GenomicRanges)

#### Load_and_Train_Model_AUC_CVglmnet_LabeledData ####
Load_and_Train_Model_AUC_CVglmnet_LabeledData <- function(classifier_all, TP_training_regions, TN_training_regions, specificity_cutoff){
  
  train <- classifier_all[1:(TP_training_regions + TN_training_regions),2:ncol(classifier_all)] 
  train_labels <- classifier_all[1:(TP_training_regions + TN_training_regions),1]
  test <- classifier_all[(TP_training_regions + TN_training_regions+1):nrow(classifier_all),2:ncol(classifier_all)] 
  test_labels <- classifier_all[(TP_training_regions + TN_training_regions+1):nrow(classifier_all),1]
  
  ENet_fit <- cv.glmnet(x=as.matrix(train[complete.cases(train),]), y= train_labels[complete.cases(train)], family = "binomial", type.measure = "auc", alpha = 0.5)
  ENet_pred_lambdamin <- predict(ENet_fit,as.matrix(test[complete.cases(test),]),s="lambda.min", type = "response") #type = response ensures that the scale is from 0 to 1 
  pred <- prediction(ENet_pred_lambdamin, test_labels[complete.cases(test)])

  #----------------------------------------------------------------------------
  sa_matrix <- matrix(0,4,1)
  #recall-precision curve from https://stackoverflow.com/questions/8499361/easy-way-of-counting-precision-recall-and-f1-score-in-r
  PR.perf <- performance(pred, "prec", "rec") 
  w <- which(PR.perf@y.values[[1]] == "NaN")
  if (length(w)>0){PR.perf@y.values[[1]][1] <- 1}
  f <- approxfun(data.frame(PR.perf@x.values , PR.perf@y.values)) #https://stackoverflow.com/questions/39268159/auc-of-a-precision-recall-curve-by-using-package-rocr  
  auprc <- integrate(f, 0, 1)$value
  #PRC <- cbind(PR.perf@x.values[[1]], PR.perf@y.values[[1]])

  perf <- performance(pred, 'sens', 'spec') #x is spec, y is sens
  ix <- which.min(abs(perf@x.values[[1]] - specificity_cutoff))
  sensitivity <- perf@y.values[[1]][ix]
  auc <- round(slot(performance(pred, measure = "auc"), "y.values")[[1]],4)

  #matthews correlation coefficient. 
  #what is cutoff for tp and tn?? 
  thresh <- mean(ENet_pred_lambdamin)
  true_values <- test_labels[complete.cases(test)]
  predicted_values <- ENet_pred_lambdamin
 
  tp <- length(which(true_values == 1 & predicted_values > thresh))
  tn <- length(which(true_values == 0 & predicted_values <= thresh))
  fp <- length(which(true_values == 0 & predicted_values > thresh))
  fn <- length(which(true_values == 1 & predicted_values <= thresh))
  
  #mcc <- (tp*tn)-(fp*fn) / ((tp+fp) * (tp+fn) * (tn + fp) * (tn + fn))^0.5 #bug in wiki
  n <- tp + fn + tn + fp
  s = (tp + fn)/n
  p = (tp + fp)/n
  mcc = (tp/n - s * p ) / (p*s*(1-s)*(1-p))^0.5

  sa_matrix[1,1] <- sensitivity
  sa_matrix[2,1] <- auc # sensitivity vs specificoty - too easy
  sa_matrix[3,1] <- auprc # preferred, not both based on negative set proportion
  sa_matrix[4,1] <- mcc 

  #newList <- list("betas" = coef(ENet_fit, s = "lambda.min"), "prediction" = ENet_pred_lambdamin, "four_metrics" = sa_matrix) #, "PRC" = PRC)
  newList <- list("four_metrics" = sa_matrix, "perf_spec_sens" = perf, "perf_prec_rec" = PR.perf)

  return(newList)
}



	#CV
	numiters <- 25
	Sens_AUC_collection <- matrix(0,4,numiters) 
	myclassifier <- read.table("/Users/kathryntsai/OneDrive\ -\ Villanova\ University/College/2018-2019/Summer\ 2019/TFs_eQTLs_Research/GitHub/eQTL_Kathryn/kathryn/bulk_data/q4_output/Classifier_binary_ENCODE_cis_v_trans_eQTLs.txt.gz",sep = "\t", header = F, stringsAsFactors = FALSE) #might need \t for classifiers that were made from data table (new version), needed " " for those made from Granges
	w <- which(is.na(myclassifier), arr.ind = T) #for tf = 254, NA is for conservation. can we go back and fix these? or are these signs that the genome is wrong? no, genome is definitely right. Perhaps conservation score couldn't be retrieved. 
	if(length(w) > 0){myclassifier <- myclassifier[-w[,1],]}
	count <- 0
        while (count < numiters){
	#replace nrow(positive_bed) with length(which(myclassifier[,1] == 1)) and nrow(negative_bed) with length(which(myclassifier[,1] == 0))
	#mix all rows 
	num_pos <- length(which(myclassifier[,1] == 1))
	num_neg <- length(which(myclassifier[,1] == 0))
        s_p <- sample(1:num_pos, size = num_pos, replace = F) #1000 at most 
        s_n <- sample(1:num_neg, size = num_neg, replace = F) #10000

        classifier_TP_train <- myclassifier[1:num_pos,]
        classifier_TN_train <- myclassifier[(num_pos + 1):nrow(myclassifier),]
        classifier <- rbind(classifier_TP_train[s_p[1:(0.8*length(s_p))],], classifier_TN_train[s_n[1:(0.8*length(s_n))],], classifier_TP_train[s_p[(0.8*length(s_p)+1):(1*length(s_p))],],classifier_TN_train[s_n[(0.8*length(s_n)+1):(1*length(s_n))],])

        tryCatch({LATM <- Load_and_Train_Model_AUC_CVglmnet_LabeledData(classifier, 0.8*length(s_p), 0.8*length(s_n), 0.99)
        count <- count + 1
        print(count)
	Sens_AUC_collection[,count]<- LATM$four_metrics
	print(Sens_AUC_collection[,count])
#        Betas <- LATM$betas #want to show that betas don't change much between samplings 
#        write.table(Betas[,1],paste0("PredictionOutput/Betas_Jan4Revisions2_",TFs[tf],"_trial",count,"_UNBALANCED_AUPRC.txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)
	#write.table(LATM$PRC, paste0("PredictionOutput/PRcurve_coords_Jan4Revisions2_",TFs[tf],"_trial",count,"_UNBALANCED_AUPRC.txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)
	
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        }#while loop
        write.table(Sens_AUC_collection, paste0("q4_output/perf_ENCODE_cisvtrans.txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)


	#entire model: up to 1,000 positive, 10,000 negative
	ENet_fit <- cv.glmnet(x=as.matrix(myclassifier[,-1]), y=as.matrix(myclassifier[,1]), family = "binomial", type.measure = "auc", alpha = 0.5)
        assign(paste0("IMPACT_fit_",tf), ENet_fit)

	w1 <- which(ls(1) == paste0("IMPACT_fit_",tf))
	save(list = ls(1)[w1], file = paste0("q4_output/IMPACT_model_ENCODE_cisvtrans.RData"), envir = .GlobalEnv)















