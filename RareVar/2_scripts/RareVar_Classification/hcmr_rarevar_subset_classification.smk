'''
Snakemake pipeline to prepare individual-level VCF files from HCMR to classify for each sarcomere positive individual what class of rare variant + the number of doses.
Author: Jonathan Chan
Date: 2023-09-13
'''

''' 
Functions
1) Run VEP on each individual-level VCF file
2) Feed in the VCF output from VEP into R script which uses Kate's OMGL classifications to compute the class of rare variants each carrier (i.e gene type & pLOF or missense) + dosage.
'''
from os import listdir
from os.path import isfile, join
from re import findall

import numpy as np  
import pandas as pd

##Define your input variables here!!!!!!!!!!!!!!!!

#For sarcomere-positive only
mypath = '/well/PROCARDIS/jchan/hcmr_ukbb/hcmr_vcf/individual_sarcpos_vcfs/' #If focusing on only sarcomere-positive individuals
pheno_geno_vcf_filenames_path= '/well/PROCARDIS/jchan/hcmr_ukbb/rarevar_subset_analysis/data/hcmr_sarc_pos_rarevar_vcf_hcrids_nogenofilter.tsv'
output_subfolder =''
sarcneg_indicator='FALSE'

#For sarcomere-negative only
# mypath = '/well/PROCARDIS/jchan/hcmr_ukbb/hcmr_vcf/individual_sarcneg_vcfs/' #If focusing on only sarcomere-negative individuals
# pheno_geno_vcf_filenames_path= '/well/PROCARDIS/jchan/hcmr_ukbb/rarevar_subset_analysis/data/hcmr_sarc_neg_rarevar_vcf_hcrids_nogenofilter.tsv'
# output_subfolder='sarcneg/'
# sarcneg_indicator='TRUE'


ids = [f for f in listdir(mypath) if isfile(join(mypath, f))] #Grabs out all files in mypath
# Uses regexp to grab only the sample ID
id_list_list = [findall('(.+).vcf', p) for p in ids] 
id_nonfiltered = [item for sublist in id_list_list for item in sublist] 

#Only use the VCF files which correspond to individuals who are both phenotyped sarcomere positive = helps define the IDs
pheno_genos_list_list=pd.read_table(pheno_geno_vcf_filenames_path, delimiter='\t', header=None).values.tolist()
pheno_genos = [item for sublist in pheno_genos_list_list for item in sublist] 

id = [x for x in id_nonfiltered if x in pheno_genos] #Only defines ID if it is present in the phenos_genos list

rule all:
    input:
        '../output/hcmr_pheno_recalibrated_vep_b38.vcf.gz',
        '../output/hcmr_pheno_recalibrated_vep_b38.vcf.gz.tbi'
        ,
        expand('../output/individual_vep/'+output_subfolder+'{id}.vcf',id=id),
        expand('../output/individual_class/'+output_subfolder+'{id}_class.tsv', id=id),
        expand('../output/individual_class/'+output_subfolder+'{id}_class2.tsv', id=id),
        expand('../output/individual_class/'+output_subfolder+'overlap_variants_{id}_class.tsv', id=id)

rule liftover:
    input: 
        vcf_input='/well/PROCARDIS/jchan/hcmr_ukbb/hcmr_vcf/hcmr_pheno_nogenofilter_recalibrated.vcf.gz'
    output:
        vcf_output='../output/hcmr_pheno_recalibrated_b38.vcf.gz'
    params:
        chain_filepath="/gpfs3/well/PROCARDIS/jchan/bin/liftover/hg19ToHg38.over.chain.gz",
        b38_fasta="/gpfs3/well/PROCARDIS/jchan/bin/liftover/hg38.fa"
    conda:
        'gms'
    shell:'''
        CrossMap.py vcf {params.chain_filepath} {input.vcf_input} {params.b38_fasta} {output.vcf_output}
    '''

rule joint_vep:
    input: rules.liftover.output.vcf_output
    output: '../output/hcmr_pheno_recalibrated_vep_b38.vcf'
    conda:
        "vep110"
    resources:
        mem_mb=8000
    shell: '''
        vep -i {input} -o {output} -offline --cache --canonical --mane --hgvs --hgvsg --vcf --dir_cache /well/PROCARDIS/jchan/bin/ensembl-vep/vep_data/ \
        --assembly GRCh38 --af_gnomadg --max_af \
        --plugin NMD \
        --plugin LoF,loftee_path:/well/PROCARDIS/cgrace/bin/loftee_b38_v2/loftee,gerp_bigwig:/well/PROCARDIS/cgrace/bin/loftee_b38_v2/files/gerp_conservation_scores.homo_sapiens.GRCh38.bw.1,human_ancestor_fa:/well/PROCARDIS/cgrace/bin/loftee_b38_v2/files/human_ancestor.fa.gz,conservation_file:/well/PROCARDIS/cgrace/bin/loftee_b38_v2/files/loftee.sql \
    '''

rule bgzip_tabix_vep:
    input:rules.joint_vep.output
    output:
        bgzip_output='../output/hcmr_pheno_recalibrated_vep_b38.vcf.gz',
        tabix_output='../output/hcmr_pheno_recalibrated_vep_b38.vcf.gz.tbi'
    conda:
        "gms"
    resources:
        mem_mb=8000
    shell:'''
        bgzip {input}
        /gpfs3/well/PROCARDIS/jchan/bin/bcftools/bcftools index -t {output.bgzip_output} -o {output.tabix_output}
    '''

rule split_vep: #Only splits out those which are in the pheno_geno_vcf_filenames_path i.e either sarc-pos or sarc-neg
    input: 
        vep_bgz = rules.bgzip_tabix_vep.output.bgzip_output
    output: 
        vep_output_folder=directory('../output/individual_vep/'+output_subfolder)
    params:
        hcrids=pheno_geno_vcf_filenames_path
    conda:
        "gms"
    resources:
        mem_mb=8000
    shell:'''
        /gpfs3/well/PROCARDIS/jchan/bin/bcftools/bcftools +split {input.vep_bgz} -S {params.hcrids} -o {output.vep_output_folder}
    '''

#NEED TO RUN THE FIRST 3 RULES FIRST BECAUSE OF ISSUE WITH INPUT NOT MATCHING WELL!
rule rscript_analysis:
    input: 
        accessory_input = rules.split_vep.output.vep_output_folder,
        main_input='../output/individual_vep/'+output_subfolder+'{id}.vcf'
    output: 
        class_file='../output/individual_class/'+output_subfolder+'{id}_class.tsv',
        class2_file='../output/individual_class/'+output_subfolder+'{id}_class2.tsv',
        overlap_variants_file='../output/individual_class/'+output_subfolder+'overlap_variants_{id}_class.tsv'
    conda:
        "gms"
    resources:
        mem_mb=16000
    params:
        kt_acmg_classifications='../data/HCMR_final_AH_150319_completevarlist.csv',
        sarcneg=sarcneg_indicator
    shell:'''
        Rscript hcmr_rarevar_classifier.R {input.main_input} {output.class_file} {output.class2_file} {output.overlap_variants_file} {params.kt_acmg_classifications} {params.sarcneg}
    '''