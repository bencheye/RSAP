# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 15:59:04 2021

@author: bench
"""
import os
import time
import pandas as pd

# 
configfile: "project_run_config.yaml"
proj_path = config['ResourcePath']['projPath']
script_path = config['ResourcePath']['ScriptPath']
end = config['SampleInfo']['END']
samples = config['SampleInfo']['Sample']
organ = config['GlobalParas']['Organ']
core_num = config['GlobalParas']['CpuCore']
ref_path = config['ResourcePath']['refPath']
reserved_item = config['GlobalParas']['reservedItem']
now_time = config['RunningTime']
alignment_dir = config['OutputDir']['ALIGNMENT']
Count_dir = config['OutputDir']['Count']

sys.path.append(script_path)
import inputBuild

rule run_countMatrix:
    input:
        inputBuild.feature_output(config)

if end == 'pair': 
    rule featureCounts_paired:
        input:
            gtf = os.path.join(ref_path, 'genes.gtf'),
            bam = expand(alignment_dir + "/{sample}.bam", sample=samples)
        output:
            os.path.join(Count_dir, 'featureCounts_countMatrix.txt')
        message:
            "Starting paired-end featureCounts analysis!"
        shell:
            "featureCounts -T {core_num} -F GTF "
            "-a {input.gtf} -o {output} "
            "-p -B -t gene -g gene_id {input.bam} "
            " 2>> Log/featureCounts_{now_time}.log"
else :
    rule featureCounts_single:
        input:
            gtf = os.path.join(ref_path, 'genes.gtf'),
            bam = expand(alignment_dir + "/{sample}.bam", sample=samples)
        output:
            os.path.join(Count_dir, 'featureCounts_countMatrix.txt')
        message:
            "Starting single-end featureCounts analysis!"
        shell:
            "featureCounts -T {core_num} -F GTF "
            "-a {input.gtf} -o {output} "
            "-t gene -g gene_id {input.bam} "
            " 2>> Log/featureCounts_{now_time}.log"

rule preprocess_countMatrix:
    input:
        os.path.join(Count_dir, 'featureCounts_countMatrix.txt')
    output:
        os.path.join(Count_dir, 'countMatrix_geneID.txt')
    message:
        "Starting preprocess count matrix!"
    shell:
        "less {input} | grep -v '^#' | "
        " cut -f1,7- | sed 's/.bam//g' | "
        " sed 's!{alignment_dir}/!!g' > {output} "
        " 2>> Log/featureCounts_{now_time}.log"

rule build_gtf_geneMapping:
    input:
        os.path.join(ref_path, 'genes.gtf')
    output:
        protected(os.path.join(ref_path, 'geneID_mapping_gene_name.txt'))
    message:
        "Starting build geneID mapping gene_name according to the gtf!"
    shell:
        "sh {script_path}/build_geneID_mapping_geneName.sh {input} {output}"

rule generate_expMatrix:
    input:
        mapping = os.path.join(ref_path, 'geneID_mapping_gene_name.txt'),
        count_matrix = os.path.join(Count_dir, 'countMatrix_geneID.txt')
    output:
        os.path.join(Count_dir, 'join_gene_expMatrix.txt')
    message:
        "join the gene expression matrix!"
    shell:
        "join -j 1 {input.mapping} {input.count_matrix} > {output} "
        " 2>> Log/featureCounts_{now_time}.log"

rule preprocess_expMatrix:
    input:
        joined_Exp_matrix = os.path.join(Count_dir, 'join_gene_expMatrix.txt'),
        group_map = os.path.join(proj_path, config['SampleInfo']['MetaData'])
    output:
        exp_matrix = os.path.join(Count_dir, 'gene_expMatrix_final.csv')
    message:
        "Build the Rdata of downstream analysis!"
    shell:
        "Rscript {script_path}/Rcode/geneCountMatrix_preprocess.R "
        "{input.joined_Exp_matrix} {input.group_map} "
        "{reserved_item} {output.exp_matrix} "
        " 2>> Log/featureCounts_{now_time}.log"
