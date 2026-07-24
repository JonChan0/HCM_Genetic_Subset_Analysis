#Script to collate the per-chromosome .REGENIE files outputted from Step2 and to convert it to _manhattan_rsid.tsv format (~SNPTEST). 
#It then writes out this .tsv for each phenotype within the .REGENIE files.
#Author: Jonathan Chan
#Date: 2024-06-04

library(tidyverse)

args <- commandArgs(trailingOnly=T)

if (length(args)==0){
  print('No arguments supplied to this script')
  stop()
}

regenie_path <- args[1]
output_prefix <- args[2]


#Import--------------------------------------------------------------------------

regenie_files <- list.files(regenie_path, pattern=str_c(output_prefix,'_\\d{1,2}\\.regenie$'))

if(length(regenie_files)!=22){
  print('Script expected 22 .regenie files with each corresponding to a chromosome')
  stop()
}

regenie_paths <- str_c(regenie_path,'/',regenie_files)
import <- map(regenie_paths, ~data.table::fread(.) %>% as_tibble())#Import in all the .regenie files across chromosomes

#Import in the Ydict files to identify what each of the Y variables corresponds to in terms of phenotype
regenie_ydicts <- list.files(regenie_path, pattern=str_c(output_prefix,'_\\d{1,2}\\.regenie.Ydict'))
ydict_paths <- str_c(regenie_path,'/',regenie_ydicts)
ydicts <- map(ydict_paths, ~read.table(., header=F, col.names=c('Y','pheno')) %>% as_tibble())

#Check if all the ydict list elements are identical
are_tibbles_identical <- function(tibble_list) {
  reference_tibble <- tibble_list[[1]]
  all(lapply(tibble_list[-1], identical, reference_tibble)) #Compares all of the other tibbles to the first tibble and sees if all are TRUE
}

identical_ydicts <- are_tibbles_identical(ydicts)

if(identical_ydicts==F){
  print('Ydict files do not correspond to the same phenotypes across the chromosomes')
  stop()
}

rm(identical_ydicts)

#Reformat-------------------------------------------------------------------------
#For each phenotype, reformat it and write it out, to save on memory (do it sequentially instead of jointly)
#The _manhattan_rsid.tsv format is 
# rsid    snid    chromosome      position        allele_A        allele_B        info    all_total       maf     pval    beta    beta_se full_snid
# rs2631652       19:56683002     19      56683002        G       A       1       2349    0.0940826       0.768802        -0.0139954      0.0476126       19:56683002_G_A
#Allele_B is the effect allele and maf refers to the minor allele generally allele_B
#In REGENIE, Allele1 is the effect allele

pheno_extractor_writer <- function(input_allpheno_list, ydict, pheno, output_path, output_prefix){
  
  y_of_interest <- ydict$Y[ydict$pheno==pheno]
  #print(y_of_interest)
  
  #Grab out the common columns to all phenotypes e.g CHROM; GENPOS; ID etc.
  common_cols <- map(input_allpheno_list, ~select(.,1:9)) %>% bind_rows()
  #Grab the phenotype_specific columns
  pheno_specific_cols <- map(input_allpheno_list, ~select(.,ends_with(y_of_interest))) %>% bind_rows()
  #print(pheno_specific_cols)
  colnames(pheno_specific_cols) <- c('BETA','SE','CHISQ','LOG10P')
  #Bind them togethers
  overall <- bind_cols(common_cols, pheno_specific_cols)
  rm(common_cols, pheno_specific_cols)
  
  #Reformat to manhattan_rsid.tsv
  overall <- overall %>%
    select(rsid=ID, chromosome=CHROM, position=GENPOS, allele_A=ALLELE0, allele_B=ALLELE1, info=INFO, all_total=N, eaf=A1FREQ, pval=LOG10P, beta=BETA, beta_se=SE)%>%
    mutate(snid = str_c(as.character(chromosome),':',as.character(position)),
           maf = ifelse(eaf >= 0.5, 1-eaf,eaf)) %>%
    dplyr::relocate(snid, .after='rsid') %>%
    dplyr::relocate(maf, .after='eaf')
  
  print(str_c('Writing out output files for ',pheno))
  
  write_tsv(overall,str_c(output_path,output_prefix,'_',pheno,'_manhattan_rsid_logpval.tsv'))                                                                                                                                              

  overall <-overall %>%
    mutate(pval=10^-pval) #Convert the -log10P version to the raw version

  write_tsv(overall,str_c(output_path,output_prefix,'_',pheno,'_manhattan_rsid.tsv'))
  
}

#Run over each phenotype
walk(ydicts[[1]]$pheno, ~pheno_extractor_writer(import, ydicts[[1]],., str_c(regenie_path,'formatted/'),output_prefix))
