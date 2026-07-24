#This script preprocesses the output from GWAMA (meta-analysis) into a format suitable for PRS-CS with QC steps included.
#Author: Jonathan Chan
#Date: 2025-12-22

library(tidyverse)

args <- commandArgs(trailingOnly = T)

# Input and output file paths
input <- args[1]  # Input GWAMA summary statistics file
output <- args[2] # Output PRS-CS formatted summary statistics file

# Read in the GWAS_SSF sumstats
sumstats <- read_tsv(input)

#QC so that if there is EAF data, remove rare variants
if(is.numeric(sumstats$effect_allele_frequency)){
  sumstats_qc <- sumstats %>%
    filter(effect_allele_frequency> 0.01, effect_allele_frequency < 0.99)
} else {
  sumstats_qc <- sumstats
}

#Rename columns that look like rsid e.g. rs_id to rsid
if('rs_id' %in% colnames(sumstats_qc) & !'rsid' %in% colnames(sumstats_qc)){
  sumstats_qc <- sumstats_qc %>%
    dplyr::rename('rsid'='rs_id')
} 

# Select PRS-CS columns
if('rsid' %in% colnames(sumstats_qc) & 'variant_id' %in% colnames(sumstats_qc)){
  sumstats_qc <- sumstats_qc %>%
    select(-variant_id) %>%
    dplyr::rename(variant_id = rsid)
} else if ('rsid' %in% colnames(sumstats_qc) & !'variant_id' %in% colnames(sumstats_qc)){
  sumstats_qc <- sumstats_qc %>%
    dplyr::rename(variant_id = rsid)
}

sumstats_out <- sumstats_qc %>%
  select(SNP = variant_id,
         A1 = effect_allele,
         A2 = other_allele,
         BETA = beta,
         SE = standard_error,
         P = p_value)

# Write out the processed summary statistics for PRS-CS
write_tsv(sumstats_out, output)