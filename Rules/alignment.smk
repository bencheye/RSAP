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
script_path = config['ResourcePath']['ScriptPath']
end = config['SampleInfo']['END']
organ = config['GlobalParas']['Organ']
core_num = config['GlobalParas']['CpuCore']
ref_path = config['ResourcePath']['refPath']
align_method = config['GlobalParas']['aligMethod']
now_time = config['RunningTime']
fastp_dir = config['OutputDir']['QC']['fastp']
alignment_dir = config['OutputDir']['ALIGNMENT']
gtf_link = config['ResourcePath']['gtfLink']
method_link = config['ResourcePath']['methodLink']

sys.path.append(script_path)
import inputBuild

rule run_aligment:
    input:
        inputBuild.alignment_output(config)

rule download_gtf:
    input:
        
    output:
        protected(os.path.join(ref_path, 'genes.gtf'))
    message:
        "Starting download {organ} gtf file!"
    shell:
        "wget -O {output}.gz {gtf_link};"
        "gunzip -c {output}.gz > {output};"

rule download_genome:
    input:
        
    output:
        os.path.join(ref_path, 'genome.fasta.gz')
    message:
        "Starting download {organ} genome file!"
    shell:
        "wget -O {output} {method_link}"
        

rule gunzip_genome:
    input:
        os.path.join(ref_path, 'genome.fasta.gz')
    output:
        protected(os.path.join(ref_path, 'genome.fasta'))
    message:
        "Starting gunzip {organ} genome file!"
    shell:
       "gunzip -c {input} > {output}" 

rule download_hisat2:
    input:

    output:
        os.path.join(ref_path, 'hisat2.tar.gz') 
    message:
        "Starting download hisat2 index!"
    shell:
        "wget -O {output} {method_link}"

rule build_hisat2_index:
    input:
        gtf = os.path.join(ref_path, 'genes.gtf'),
        index = os.path.join(ref_path, 'hisat2.tar.gz')
    output:
        protected(directory(os.path.join(ref_path, 'hisat2')))
    message:
        "Starting download {organ} hisat2 index!"
    shell:
        "tar -zxvf {input.index} -C {ref_path};"
        "mv {ref_path}/grch38 {output}"

rule build_star_index:
    input:
        fasta = os.path.join(ref_path, 'genome.fasta'),
        gtf = os.path.join(ref_path, 'genes.gtf')
    output:
        protected(directory(os.path.join(ref_path, 'STAR')))
    message:
        "Starting build {organ} {align_method} alignment index"
    shell:
        "STAR  --runMode genomeGenerate "
        "--runThreadN {core_num} "
        "--genomeDir {output} "
        "sjdbGTFfile {input.gtf} "
        "--genomeFastaFiles {input.fasta}"

if align_method == 'STAR':
    if end == 'pair':  
        rule star:
            input:
                R1 = os.path.join(fastp_dir, "{sample}_R1_clean.fq.gz"), 
                R2 = os.path.join(fastp_dir, "{sample}_R2_clean.fq.gz"),
                gtf = os.path.join(ref_path, 'genes.gtf'),
                ref_genome = os.path.join(ref_path, 'STAR')
            output:
                os.path.join(alignment_dir, "{sample}.bam") 
            message:
                "{wildcards.sample} : alignment-star starting running!"
            params:
                "{alignment_dir}/{sample}"
            shell:
                "STAR --genomeDir {input.ref_genome} "
                "--runThreadN {core_num} "
                "--readFilesIn {input.R1} {input.R2} "
                "--readFilesCommand zcat "
                "--outSAMtype BAM SortedByCoordinate "
                "--outBAMsortingThreadN {core_num} "
                "--outFileNamePrefix {params} 2>> Log/alignment_star_{now_time}.log"
    else:            
        rule star:
            input:
                seq_fq = os.path.join(fastp_dir, "{sample}_clean.fq.gz"),
                gtf = os.path.join(ref_path, 'genes.gtf'),
                ref_genome = os.path.join(ref_path, 'STAR')                        
            output:
                os.path.join(alignment_dir, "{sample}.bam")
            message:
                "{wildcards.sample}: alignment-star starting running!"
            params:
                "{alignment_dir}/{sample}"
            shell:
                "STAR --genomeDir {input.ref_genome} "
                "--runThreadN {core_num} "
                "--readFilesIn {input.seq_fq} "
                "--readFilesCommand zcat "
                "--outSAMtype BAM SortedByCoordinate "
                "--outBAMsortingThreadN {core_num} "
                "--outFileNamePrefix {params} 2>> Log/alignment_star_{now_time}.log"
elif align_method == 'hisat2':
    if end == 'pair':  
        rule hisat2:
            input:
                R1 = os.path.join(fastp_dir, "{sample}_R1_clean.fq.gz"), 
                R2 = os.path.join(fastp_dir, "{sample}_R2_clean.fq.gz"),
                index_path = os.path.join(ref_path, 'hisat2')
            output:
                os.path.join(alignment_dir, "{sample}.sam")
            message:
                "{wildcards.sample} : alignment-star starting running!"
            shell:
                "hisat2 -p {core_num} -x {input.index_path}/genome "
                "-1 {input.R1} "
                "-2 {input.R2} "
                "-S {output} 2>> Log/alignment_hisat2_{now_time}.log"
    else:            
        rule hisat2:
            input:
                data = os.path.join(fastp_dir, "{sample}_clean.fq.gz"),
                index_path = os.path.join(ref_path, 'hisat2')          
            output:
                os.path.join(alignment_dir, "{sample}.sam")
            message:
                "{wildcards.sample}: alignment-hisat2 starting running!"
            shell:
                "hisat2 -p {core_num} -x {input.index_path}/genome "
                "-U {input.data} -S {output} 2>> Log/alignment_hisat2_{now_time}.log"
    rule samtools:
        input:
            os.path.join(alignment_dir, "{sample}.sam")
        output:
            bam = os.path.join(alignment_dir, "{sample}.bam"),
            view_bam = temp(os.path.join(alignment_dir, "{sample}_view.bam"))
        message:
            "{wildcards.sample}: samtools starting running!"
        shell:
            "samtools view -bS {input} > {output.view_bam};"
            "samtools sort -@ {core_num} {output.view_bam} -o {output.bam};"
            "samtools index {output.bam} 2>> Log/alignment_samtools_{now_time}.log"
