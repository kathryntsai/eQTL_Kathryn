# ==========================================================
# Trans GeneCards
# ==========================================================

trans_gene_cards <- data.table(snps_interesting_TFmatch[,unique(snps_interesting_TFmatch$transGene)], paste("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", trim(snps_interesting_TFmatch[,unique(snps_interesting_TFmatch$transGene)]), sep=""))

# Install biomaRt
# source("https://bioconductor.org/biocLite.R")
# biocLite("biomaRt")

# https://www.biostars.org/p/178726/
library("biomaRt")

mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)
ens <- snps_interesting_TFmatch[,unique(snps_interesting_TFmatch$transGene)]
ensLookup <- gsub("\\.[0-9]*$", "", ens) # not necessary

annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_transcript_id", "ensembl_gene_id", "gene_biotype", "external_gene_name"),
  filter="ensembl_gene_id",
  values=ensLookup,
  uniqueRows=TRUE)

annotLookup2 <- data.frame(
  ens[match(annotLookup$ensembl_gene_id, ensLookup)],
  annotLookup)

colnames(annotLookup2) <- c(
  "original_id",
  c("ensembl_transcript_id", "ensembl_gene_id", "gene_biotype", "external_gene_name"))
annotLookup_abbrev <- data.table(unique(annotLookup2[, "original_id"]), unique(annotLookup2[, "external_gene_name"]))
fwrite(annotLookup_abbrev, "q2_1_written/trans_gene_commonname.csv")