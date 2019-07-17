## LOAD

library(devtools)
devtools::load_all("/Users/kathryntsai/Documents/GitHub/finemapr")

library(magrittr)
library(dplyr)
library(ggplot2)
library(knitr)
library(ggthemes)
theme_set(theme_stata())

for(i in 1:3){#length(chr1_cis_split)){
  ## WHICH GENE

  chr1_cis <- franke_cis_data[which(GeneChr==1),]
  chr1_cis_split <- split(chr1_cis, chr1_cis$GeneSymbol)
  this_gene <- chr1_cis_split[[i]] 
  
  if (nrow(this_gene) >= 1000) {
    print(paste0(this_gene$GeneSymbol[1]))
    next
  }
  
  my_path="/Users/kathryntsai/OneDrive\ -\ Villanova\ University/College/2018-2019/Summer\ 2019/TFs_eQTLs_Research/RProjects/Project2/finemapping/"
  
  ## ZSCORES
  
  new <- this_gene[, c("SNP", "Zscore")]
  colnames(new) <- c("snp", "zscore")
  fwrite(new, paste0(my_path, this_gene$GeneSymbol[1], "_z.z"), sep="\t", col.names=F)
  z1 <- read_zscore(paste0(my_path, this_gene$GeneSymbol[1], "_z.z"))
  
  ## LD MATRIX
  ldin <- glue_collapse(this_gene$SNP, sep="\\n") 
  #sysin <- noquote(gsub("\\n", "\n", gsub('\\\'', "\'", paste0('curl -k -H "Content-Type: application/json" -X POST -d \'{', '"snps"', ": \"", ldin, '", "pop": "CEU","r2_d": "r2"}\'', " 'https://ldlink.nci.nih.gov/LDlinkRest/ldmatrix?&token=bc2543d7f1e5' > ", "'", my_path, this_gene$GeneSymbol[1], "_ld.ld'"), str, fixed=T), str, fixed=T))
  #z <- system(sysin) 
  system(paste0('curl -k -H "Content-Type: application/json" -X POST -d \'{', '"snps"', ": \"", ldin, '", "pop": "CEU","r2_d": "r2"}\'', " 'https://ldlink.nci.nih.gov/LDlinkRest/ldmatrix?&token=bc2543d7f1e5' > ", "'", my_path, this_gene$GeneSymbol[1], "_ld.txt'"))
  
  ld1 <- read_ld(paste0(my_path, this_gene$GeneSymbol[1], "_ld.txt"))
  
  ## CLEAN THIS
  
  z1 <- z1[match(row.names(ld1), z1$snp),]
  
  ld1 <- ld1[z1$snp %in% rownames(ld1), z1$snp %in% colnames(ld1)] # ensure matrix is square
  z1 <- z1[match(row.names(ld1), z1$snp),]
  
  ld1 <- ld1[!is.na(round(diag(ld1), 4)),!is.na(round(diag(ld1), 4))] # ensure diagonal exists for caviar
  z1 <- z1[match(row.names(ld1), z1$snp),]
  
  n1 <- nrow(z1)
  
  ## EXPLORE
  
  z1 %>% arrange(-abs(as.numeric(zscore))) %>% head(5) %>% kable(digits = 1)
  ggplot(z1, aes(zscore)) + geom_histogram(fill='red', color='blue')
  mutate(z1, pval = pchisq(zscore^2, df = 1, lower.tail = FALSE)) %>%
    ggplot(aes(pval)) + geom_histogram(fill='red', color='blue')
  
  ## RUN FINEMAP
  
  options(finemapr_finemap = "/Users/kathryntsai/finemap_v1.1_MacOSX/finemap_v1.1_MacOSX")
  
  out_finemap <- run_finemap(z1, ld1, n1, args = "--n-causal-max 3") # if not using test, make it z1
  out_finemap$num_loci <- length(out_finemap$snp[["snp"]])
  out_finemap$dir_run <- paste0(getwd(), '/run_finemap')
  print(out_finemap)
  plot(out_finemap, label_size = 3, grid_ncol = 1)
  
  ## RUN CAVIAR
  
  options(finemapr_caviar = "/Users/kathryntsai/Documents/GitHub/caviar/CAVIAR-C++/CAVIAR")
  out_caviar <- run_caviar(z1, ld1, args = "-c 3")
  print(out_caviar)
  plot(out_caviar, label_size = 3)
  
  ## RUN PAINTOR
  
  options(finemapr_paintor = "/Users/kathryntsai/Documents/GitHub/PAINTOR_V2.1/PAINTOR")
  out_paintor <- run_paintor(z1, ld1, n1, args = "-enumerate 3")
  out_paintor$num_loci <- length(out_paintor$snp[["snp"]])
  out_paintor <- collect_results(out_paintor)
  print(out_paintor)
  plot(out_paintor, label_size = 3)

}