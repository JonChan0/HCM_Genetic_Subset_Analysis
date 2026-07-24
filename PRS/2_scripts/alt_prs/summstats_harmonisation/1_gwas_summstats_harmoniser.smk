'''
This outlines the Snakemake pipeline to harmonise all input GWAS summary statistics using EBI GWAS Catalog's GWAS Summstats Harmoniser.
This version includes a conditional rule to reformat specific GWAS outputs instead of harmonising them.
Author: Jonathan Chan
Date: 2025-02-10 (Modified 2026-02-02)
'''

import os
import pandas as pd
import hashlib
import yaml

configfile: 'config_gwas_summstats_harmoniser.yaml'

# --- Input File Partitioning ---
# This section separates input files into two groups based on their path.

all_input_paths = config['input_tsv_gz_files']

# Standard files to be harmonised which have the format of phenotype.tsv.gz
standard_paths = [p for p in all_input_paths]
standard_folders = [os.path.dirname(p) for p in standard_paths]
standard_basenames = [os.path.splitext(os.path.splitext(os.path.basename(p))[0])[0] for p in standard_paths]

# Create a complete list of all basenames and a map to find their original folder
all_basenames = standard_basenames
all_folders = standard_folders
basename_to_folder_map = dict(zip(all_basenames, all_folders))

common_output = [
    expand(config['base_output_folder']+'{basename}/final/{basename}.h.tsv.gz', basename=all_basenames),
    config['desired_output_folder']+'harmonisation_summary.tsv'
]

rule all:
    input: 
         common_output

rule gwas_summstats_harmoniser:
    input:
        # The lambda function uses the map to find the correct input path for any given basename.
        gwas_summstats=lambda wildcards: f"{basename_to_folder_map[wildcards.basename]}/{wildcards.basename}.tsv.gz"
    output:
        harmonised_output_folder = directory(config['base_output_folder']+'{basename}/'),
        harmonised_gwas_summstats = config['base_output_folder']+'{basename}/final/{basename}.h.tsv.gz',
        output_logfiles = config['base_output_folder']+'{basename}/final/{basename}.running.log'
    resources:
        mem_mb = 40000
    params:
        ref_folder = config['gwas_summstats_harmoniser_ref_folder'],
        output_folder = config['desired_output_folder']
    shell:
        '''
        module load Nextflow/24.04.2

        nextflow run EBISPOT/gwas-sumstats-harmoniser -r v1.1.10 \
        --ref {params.ref_folder} \
        --file {input.gwas_summstats} \
        --harm \
        --chromlist 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 \
        --profile executor,conda \
        --to_build '38'

        echo Successfully completed harmonisation of {input.gwas_summstats}
        '''

# Rule = Parse all '.running.log' files to generate a single summary report
rule summarise_harmonisation:
    input:
        # Collect log files from ALL harmonisation jobs.
        lHCMR_files = expand(config['base_output_folder']+'{basename}/final/{basename}.running.log', basename=all_basenames)
    output:
        summary_file = config['desired_output_folder']+'harmonisation_summary.tsv'
    run:
        summary_data = []
        for lHCMR_file in input.lHCMR_files:
            basename = os.path.basename(lHCMR_file).split('_gwas_ssf')[0]
            with open(lHCMR_file, 'r') as f:
                lines = f.readlines()
                success = any("Result\tSUCCESS_HARMONIZATION" in line for line in lines)
                if success:
                    percentage_line = [line for line in lines if 'sites successfully harmonised' in line]
                    clean_percentage = percentage_line[0].strip() if percentage_line else "N/A"
                    summary_data.append(f"{basename}\tYes\t{clean_percentage}")
                else:
                    summary_data.append(f"{basename}\tNo\t(N/A)")
        
        with open(output.summary_file, 'w') as f:
            f.write("Phenotype\tHarmonised\tPercentage Harmonised\n")
            for line in sorted(summary_data): # Sort for consistent output
                f.write(line + "\n")
