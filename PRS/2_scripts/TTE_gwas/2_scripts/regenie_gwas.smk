'''
Snakemake script to run REGENIE for TTE GWAS in HCMR
'''

configfile: "config.yaml"

rule all:
    input:
        config['step1_output_path']+config['output_prefix']+"_pred.list"
        # , #The output from the first step of REGENIE
        # expand("{path}{output_prefix}_{chr}.regenie",path=config['step2_output_path'],output_prefix=config['output_prefix'],chr=config['chromosomes'])

rule regenie_step1:
    input:
        ldpruned_hardcalled_snps=config['hardcalled_genotypes_allchr_ldpruned2_plink_prefix']+'.bed', 
        phenoFile=config['pheno_file'],
        covarFile=config['covar_file']
    output:
        map_file=config['step1_output_path']+config['output_prefix']+"_pred.list"
    resources:
        mem_mb=32000
    threads: 8
    conda:
        'regenie4.1.2_env'
    params:
        hardcalled_snps_allchr_ldpruned2_path=config['hardcalled_genotypes_allchr_ldpruned2_plink_prefix'],
        nThreads=8,
        output_prefix=config['step1_output_path']+config['output_prefix'],
        optional_params=config['step1_optional_params'],
        phenoColList=config['event_times'],
        eventColList=config['event_indicators'],
        cont_covariates=config['cont_covariates'],
        cat_covariates=config['cat_covariates']
    shell:'''
        regenie \
        --step 1 \
        --bed {params.hardcalled_snps_allchr_ldpruned2_path} \
        --bsize 1000 \
        --phenoFile {input.phenoFile} \
        --covarFile {input.covarFile} \
        --cv 2 \
        --t2e \
        --phenoColList {params.phenoColList} \
        --eventColList {params.eventColList} \
        --covarColList {params.cont_covariates} \
        --catCovarList {params.cat_covariates} \
        --threads {params.nThreads} \
        --out {params.output_prefix} \
        {params.optional_params}
    '''

# rule regenie_step2_singleAssoc:
#     input:
#         phenoFile=config['pheno_file'],
#         covarFile=config['covar_file'],
#         step1_map_file=rules.regenie_step1.output.map_file,
#         imputed_snps=config['imputed_genotypes_plink_prefix']+'.bed'
#     output:
#         config['step2_output_path']+config['output_prefix']+"_{chr}.regenie"
#     resources:
#         mem_mb=32000
#     threads: 8
#     conda:
#         'regenie4.1.2_env'
#     params:
#         nThreads=8,
#         output_prefix=config['step2_output_path']+config['output_prefix']+"_{chr}",
#         plink_prefix=config['imputed_genotypes_plink_prefix'],
#         phenoColList=config['event_times'],
#         eventColList=config['event_indicators'],
#         cont_covariates=config['cont_covariates'],
#         cat_covariates=config['cat_covariates']
#     shell:'''
#         regenie \
#         --step 2 \
#         --bed {params.plink_prefix}
#         --covarFile {input.covarFile} \
#         --phenoFile {input.phenoFile} \
#         --minINFO 0.7 \
#         --minMAC 50 \
#         --bsize 1000 \
#         --t2e \
#         --phenoColList {params.phenoColList} \
#         --eventColList {params.eventColList} \
#         --covarColList {params.cont_covariates} \
#         --catCovarList {params.cat_covariates} \
#         --pred {input.step1_map_file} \
#         --firth --approx \
#         --out {params.output_prefix} \
#         --no-split \
#         --threads {params.nThreads}
#     '''

# rule chr_merger_pheno_extractor:
#     input:
#         expand("{path}{output_prefix}_{chr}.regenie",path=config['step2_output_path'],output_prefix=config['output_prefix'],chr=config['chromosomes']) #The output from the second step of REGENIE
#     output:
#         expand("{path}formatted/{output_prefix}_{pheno}_manhattan_rsid.tsv",path=config['step2_output_path'],pheno=config['phenotypes'],output_prefix=config['output_prefix']),
#         expand("{path}formatted/{output_prefix}_{pheno}_manhattan_rsid_logpval.tsv",path=config['step2_output_path'],pheno=config['phenotypes'],output_prefix=config['output_prefix'])
#     params:
#         regenie_path = config['step2_output_path'],
#         output_prefix=config['output_prefix']
#     resources:
#         mem_mb=64000
#     shell:'''
#         module purge
#         module load R/4.3.2-gfbf-2023a

#         Rscript REGENIE_collater_formatter.R {params.regenie_path} {params.output_prefix}
#     '''