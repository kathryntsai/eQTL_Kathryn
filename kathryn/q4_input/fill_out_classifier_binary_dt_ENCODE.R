#List of Functions
args = commandArgs(trailingOnly=TRUE)
tf <- as.numeric(args)
library(GenomicRanges)
len_datasets <- 500
chrs <- c(1:22,"X","Y")

get_overlap <- function(index,chr,mydat){ #for each gene in list of ENSG_IDs 
  if (index %% 100 == 0){print(index)}
  #my_chrom <- my_chroms[index] #e.g. 1, not chr1
  my_chrom <- chr
  my_start <- mydat$start[index]
  my_end <- mydat$end[index]
  len_gene <- my_end-my_start
  s <- sapply(1:len_datasets, function(x) compute_overlap(x,my_chrom,my_start, my_end,len_gene))
  return(s)
}

compute_overlap <- function(i,my_chrom,my_gene_start,my_gene_end,len_gene){

  prop_full <- 0
  prop_part <- 0
  prop_over <- 0
  dat <- get(paste0("dataset_",i,"_chr",my_chrom)) #this is a data table
  colnames(dat) <- c("start","end")
  ## bug fix ## 
  w_swap <- which(dat$end < dat$start)
  if(length(w_swap) > 0){
    dat_new_end <- dat$start[w_swap]
    dat_new_start <- dat$end[w_swap]
    dat$start[w_swap] <- dat_new_start
    dat$end[w_swap] <- dat_new_end
  }

  w <- which(dat$start >= my_gene_start & dat$end <= my_gene_end)
  prop_full <- sum(dat$end[w]-dat$start[w])/len_gene

  #partial 
  w_l <- which(dat$start < my_gene_start & dat$end > my_gene_start & dat$end < my_gene_end)
  w_r <- which(dat$start < my_gene_end & dat$end > my_gene_end & dat$start > my_gene_start)
  if (length(c(w_l,w_r)) > 0){
    prop1 <- sum(dat$end[w_l]-my_gene_start)/len_gene
    prop2 <- sum(my_gene_end-dat$start[w_r])/len_gene
    prop_part <- prop1+prop2
  }else{prop_part <- 0}

  #over
  w_over <- which(dat$start < my_gene_start & dat$end > my_gene_end)
  if (length(w_over) > 0){prop_over <- 1}

  prop <- prop_full+prop_part+prop_over
  return(prop)
}

        #these text files should have 3 columns and 10K rows each and look like 
        #chr1 15672 15673 (if the SNP is at 72)
        positive_bed <- read.table("/Users/kathryntsai/OneDrive - Villanova University/College/2018-2019/Summer 2019/TFs_eQTLs_Research/RProjects/Project2/q4_output/franke_cis_bed_file_smaller.txt",sep = "\t", header = F, stringsAsFactors = FALSE)
        negative_bed <- read.table("/Users/kathryntsai/OneDrive - Villanova University/College/2018-2019/Summer 2019/TFs_eQTLs_Research/RProjects/Project2/q4_output/franke_trans_bed_file_smaller.txt",sep = "\t", header = F, stringsAsFactors = FALSE)

        head(positive_bed)
        head(negative_bed)
        dim(positive_bed)
        dim(negative_bed)
        
        training_data <- rbind(positive_bed, negative_bed)
        colnames(training_data) <- c('chr','start','end')
        my_chroms <- sapply(1:nrow(training_data), function(x) strsplit(training_data[x,1], split = "chr")[[1]][2])

	      numextrafeatures_pluslabel <- 1 #label (1 if cis, 0 if trans)

        #build feature matrix for positive and negative sets, then sample from each. 
        classifier_TP_train <- matrix(0,nrow(positive_bed),(len_datasets+numextrafeatures_pluslabel))
        classifier_TP_train[,1] <- 1
        classifier_TN_train <- matrix(0,nrow(negative_bed),(len_datasets+numextrafeatures_pluslabel))
        classifier_TN_train[,1] <- 0
        classifier_train_empty_besides2col <- rbind(classifier_TP_train,classifier_TN_train)

	      chromosomes <- c(1:22,"X","Y")
        for (chrom in 1:length(chromosomes)){
                #load chr RData one at a time 
                load(paste0("q4_input/fullIMPACT_ENCODE_dataset_partitioned_for_chr",chromosomes[chrom],"_datatable.RData"))
                print(paste0("chr",chromosomes[chrom]))
                #run
                w <- which(my_chroms == chromosomes[chrom]) #no chr in var name 
                print(paste0("total segments to be analyzed: ", length(w)))
                if (length(w)>0){
                part_classifier_impact <- matrix(0,length(w),len_datasets)
                for (g in 1:length(w)){
                  #print(g)
                  part_classifier_impact[g,] <- get_overlap(g,chromosomes[chrom], training_data[w,])
                }

                #fill it in! 
                classifier_train_empty_besides2col[w,(numextrafeatures_pluslabel+1):ncol(classifier_train_empty_besides2col)] <- part_classifier_impact
                }#end if len w > 0

                remove(list = ls(pattern = "dataset_"))
        }
        write.table(classifier_train_empty_besides2col, "q4_output/Classifier_binary_ENCODE_cis_v_trans_eQTLs.txt", sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)
        system("gzip q4_output/Classifier_binary_ENCODE_cis_v_trans_eQTLs.txt")

