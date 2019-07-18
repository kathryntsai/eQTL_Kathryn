library(data.table)
library(glue)
library(devtools)
franke_cis_data <- fread("input_data/cis-eQTL_significant_20181017.txt", header = TRUE, sep = "\t", dec = ".")
devtools::load_all("/data/srlab/ktsai/finemapr")
my_path <- "/data/srlab/ktsai/"

for (j in 1:1){
  
  print(paste0("j: ", j))
  chr1_cis <- franke_cis_data[which(GeneChr==j),]
  chr1_cis_split <- split(chr1_cis, chr1_cis$GeneSymbol)
  
  for(i in 1:length(chr1_cis_split)){
    
    this_gene <- chr1_cis_split[[i]]
    if (nrow(this_gene) >= 1000) {
      print(paste0(this_gene$GeneSymbol[1]))
      next
    }
  
    z1 <- read_zscore(paste0(my_path, "chr", j, "/", this_gene$GeneSymbol[1], "_z.z"))
    ld1 <- read_ld(paste0(my_path, "chr", j, "/", this_gene$GeneSymbol[1], "_ld.txt"))
    
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
    
    options(finemapr_finemap = "/data/srlab/ktsai/finemap_v1.1_MacOSX/finemap_v1.1_MacOSX")
    
    out_finemap <- run_finemap(z1, ld1, n1, args = "--n-causal-max 3") # if not using test, make it z1
    out_finemap$num_loci <- length(out_finemap$snp[["snp"]])
    out_finemap$dir_run <- paste0(my_path, 'run_finemap')
    print(out_finemap)
    saveRDS(out_finemap, paste0(my_path, "chr", j, "/",  this_gene$GeneSymbol[1], "_finemap.rds"))
    plot(out_finemap, label_size = 3, grid_ncol = 1)
    
    ## RUN CAVIAR
    
    options(finemapr_caviar = "/data/srlab/ktsai/caviar/CAVIAR-C++/CAVIAR")
    out_caviar <- run_caviar(z1, ld1, args = "-c 3")
    print(out_caviar)
    saveRDS(out_caviar, paste0(my_path, "chr", j, "/", this_gene$GeneSymbol[1], "_caviar.rds"))
    plot(out_caviar, label_size = 3)
    
    ## RUN PAINTOR
    
    options(finemapr_paintor = "/data/srlab/ktsai/PAINTOR_V2.1/PAINTOR")
    out_paintor <- run_paintor(z1, ld1, n1, args = "-enumerate 3")
    out_paintor$num_loci <- length(out_paintor$snp[["snp"]])
    out_paintor <- collect_results(out_paintor)
    print(out_paintor)
    saveRDS(out_paintor, paste0(my_path, "chr", j, "/", this_gene$GeneSymbol[1], "_paintor.rds"))
    plot(out_paintor, label_size = 3)
  
  }
}