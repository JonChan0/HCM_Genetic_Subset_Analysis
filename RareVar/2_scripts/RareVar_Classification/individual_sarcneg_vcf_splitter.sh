#This file iterates over chunks to split the overall VCF into individiual vcfs for sarcomere-negative individuals.

input_IDs="/gpfs3/well/PROCARDIS/jchan/hcmr_ukbb/rarevar_subset_analysis/data/hcmr_sarc_neg_rarevar_vcf_hcrids_nogenofilter.tsv"
input_vcf="hcmr_pheno_nogenofilter_recalibrated.vcf.gz"

#split -l 1000 --additional-suffix=.tsv $input_IDs sarcneg_IDs_split


#/gpfs3/well/PROCARDIS/jchan/bin/bcftools/bcftools +split $input_vcf -S sarcneg_IDs_splitaa.tsv -o individual_sarcneg_vcfs


#/gpfs3/well/PROCARDIS/jchan/bin/bcftools/bcftools +split $input_vcf -S sarcneg_IDs_splitab.tsv -o individual_sarcneg_vcfs

#It also splits the VEP file of the overall HCMR into individual_vep/sarcneg

input_vep="/gpfs3/well/PROCARDIS/jchan/hcmr_ukbb/rarevar_subset_analysis/output/hcmr_pheno_recalibrated_vep.vcf.gz"
output_folder="/gpfs3/well/PROCARDIS/jchan/hcmr_ukbb/rarevar_subset_analysis/output/individual_vep/sarcneg"

/gpfs3/well/PROCARDIS/jchan/bin/bcftools/bcftools +split $input_vep -S sarcneg_IDs_splitaa.tsv -o $output_folder

/gpfs3/well/PROCARDIS/jchan/bin/bcftools/bcftools +split $input_vep -S sarcneg_IDs_splitab.tsv -o $output_folder
