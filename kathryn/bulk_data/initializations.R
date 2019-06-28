library(data.table)
library(dplyr)
library(glue)
# Read Emma Davenport's eQTL R-variables - not used
# load("input_data/files_for_pfizer_eqtl.rda")

# Read Emma Davenport's Pfizer Data
davenport_data <- read.csv("input_data/pfizer_eqtl_table_1.csv")
pfizer_data <-davenport_data
# Read Lude Franke's cis-eQTL data
franke_cis_data <- fread("input_data/cis-eQTL_significant_20181017.txt", header = TRUE, sep = "\t", dec = ".")

# Read Lude Franke's trans-eQTL data
franke_trans_data <- fread("input_data/trans-eQTL_significant_20181017.txt", header = TRUE, sep = "\t", dec = ".")

library(fastmatch)
library(seqinr)
library(tidyverse)
library(stringr)