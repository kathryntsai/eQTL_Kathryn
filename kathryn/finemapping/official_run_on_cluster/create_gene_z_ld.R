# n=1;
# max=22;
# while [ "$n" -le "$max" ]; do
# rm "my_analysis_$n.Rout"
# n=`expr "$n" + 1`;
# done

library(data.table)
library(glue)

# my_path <- paste0(getwd(), "/") 
franke_cis_data <- fread("input_data/cis-eQTL_significant_20181017.txt", header = TRUE, sep = "\t", dec = ".")

my_path <- "/data/srlab/ktsai/"
print("this works")

for(j in 2:22){
  
  print(paste0("j: ", j))
  chr1_cis <- franke_cis_data[which(GeneChr==j),]
  chr1_cis_split <- split(chr1_cis, chr1_cis$GeneSymbol)
  
  for(i in 1:length(chr1_cis_split)){
    print(paste0("i: ", i))
    ## WHICH GENE
    
    this_gene <- chr1_cis_split[[i]]
    print(paste0(colnames(this_gene), ": ", this_gene[1,]))
    
    if (nrow(this_gene) >= 1000) {
      print(paste0(this_gene$GeneSymbol[1]))
      next
    }
    
    ## ZSCORES
    
    new <- this_gene[, c("SNP", "Zscore")]
    colnames(new) <- c("snp", "zscore")
    fwrite(new, paste0(my_path, "chr", j, "/", this_gene$GeneSymbol[1], "_z.z"), sep="\t", col.names=F)
    
    ## LD MATRIX
    ldin <- glue_collapse(this_gene$SNP, sep="\\n")
    system(paste0('curl -k -H "Content-Type: application/json" -X POST -d \'{', '"snps"', ": \"", ldin, '", "pop": "CEU","r2_d": "r2"}\'', " 'https://ldlink.nci.nih.gov/LDlinkRest/ldmatrix?&token=bc2543d7f1e5' > ", "'", my_path, "chr", j, "/",  this_gene$GeneSymbol[1], "_ld.txt'"))
  }
  
}  

