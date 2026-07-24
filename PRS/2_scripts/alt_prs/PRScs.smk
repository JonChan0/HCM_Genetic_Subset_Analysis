################################################################################
# PRS-CS training pipeline (Snakemake)
################################################################################

import os
import pandas as pd

configfile: 'config.yaml'

#Import in the phenotypes from the config_phenos.csv
config_phenos = pd.read_csv('config_phenos.csv', header=0)
phenos = config_phenos['pheno']
sumstats_in = config_phenos['sumstats_in']
gwas_n = config_phenos['gwas_n']

# Turn the phenotype table into plain Python lists (safe for expand/zip usage)
PHENOS = phenos.astype(str).tolist()
SUMSTATS_IN = sumstats_in.astype(str).tolist()
GWAS_N = gwas_n.tolist()

print(PHENOS)

rule all:
    input:
        expand(config["output_dir"] + "1_processed_sumstats/" + "sumstats_prscs_{p}.tsv",
        p=PHENOS),
        expand(config["output_dir"] + "2_PRSCS_results/" + "prscs_weights_allchr_{p}.txt",
        p=PHENOS),
        expand(config["output_dir"] + "3_hcmr_scores/prs_cs_scores_{p}.sscore",
        p=PHENOS)

# Run the pipeline per-phenotype, wiring sumstats and n_gwas via zip lists
rule preprocess_sumstats:
    input:
        sumstats=lambda wc: dict(zip(PHENOS, SUMSTATS_IN))[wc.p]
    output:
        prscs=config["output_dir"] + "1_processed_sumstats/" + "sumstats_prscs_{p}.tsv"
    resources:
        mem_mb=16000
    log:
        config["output_dir"] + "logs/" + "preprocess_sumstats_{p}.log"
    shell: """
        module load R/4.3.2-gfbf-2023a
        Rscript aux_scripts/1_summstat_processor.R {input.sumstats} {output.prscs} > {log} 2>&1
    """

rule train_prscs:
    input:
        sst=rules.preprocess_sumstats.output.prscs,
        bim=config["target_dataset_plinkprefix"] + ".bim"
    output:
        done=config["output_dir"] + "2_PRSCS_results/" + "{p}/training.done"
    params:
        prscs_py=config["prscs_py"],
        ref_dir=config["ld_ref_dir"],
        plink_prefix=config["target_dataset_plinkprefix"],
        n_gwas=lambda wc: dict(zip(PHENOS, GWAS_N))[wc.p],
        out_dir=config["output_dir"] + "2_PRSCS_results/" + "{p}/",
        seed=config["seed"],
        n_threads=4
    resources:
        mem_mb=64000
    log:
        config["output_dir"] + "logs/" + "train_prscs_{p}.log"
    conda: "python3.11_ml"
    shell: """
        N_THREADS={params.n_threads}
        export MKL_NUM_THREADS=$N_THREADS
        export NUMEXPR_NUM_THREADS=$N_THREADS
        export OMP_NUM_THREADS=$N_THREADS

        python {params.prscs_py} \
          --ref_dir={params.ref_dir} \
          --bim_prefix={params.plink_prefix} \
          --sst_file={input.sst} \
          --n_gwas={params.n_gwas} \
          --out_dir={params.out_dir} \
          --seed={params.seed} \
          > {log} 2>&1
        touch {output.done}
    """

rule perchr_prscs_weights_merger:
    input:
        done=rules.train_prscs.output.done
    output:
        weights=config["output_dir"] + "2_PRSCS_results/" + "prscs_weights_allchr_{p}.txt"
    params:
        out_dir=config["output_dir"] + "2_PRSCS_results/" + "{p}/"
    resources:
        mem_mb=16000
    log:
        config["output_dir"] + "logs/" + "perchr_prscs_weights_merger_{p}.log"
    shell: """
        cat {params.out_dir}/*.txt > {output.weights} 2> {log}
    """

rule hcmr_evaluator:
    input:
        weights=rules.perchr_prscs_weights_merger.output.weights,
        bedfile=config["target_dataset_plinkprefix"] + ".bed",
        bimfile=config["target_dataset_plinkprefix"] + ".bim",
        famfile=config["target_dataset_plinkprefix"] + ".fam"
    output:
        scores=config["output_dir"] + "3_hcmr_scores/" + "prs_cs_scores_{p}.sscore"
    resources:
        mem_mb=32000
    log:
        config["output_dir"] + "logs/" + "hcmr_evaluator_{p}.log"
    params:
        input_bfile_prefix=config["target_dataset_plinkprefix"],
        out_dir_fileprefix=config["output_dir"] + "3_hcmr_scores/" + "prs_cs_scores_{p}"
    shell: """
        /well/PROCARDIS/jchan/bin/plink2 --bfile {params.input_bfile_prefix} \
          --score {input.weights} 2 4 6 no-mean-imputation \
          --out {params.out_dir_fileprefix} > {log} 2>&1
    """

# Constrain wildcard values to your phenotype list
wildcard_constraints:
    p="|".join(PHENOS)
