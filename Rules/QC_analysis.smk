# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 15:59:04 2021

@author: bench
"""
import os
import time

configfile: "project_run_config.yaml"
# sampleInfo
script_path = config['ResourcePath']['ScriptPath']
proj_path = config['ResourcePath']['projPath']
organ = config['GlobalParas']['Organ']
end = config['SampleInfo']['END']
samples = config['SampleInfo']['Sample']
end_seps = config['SampleInfo']['PairedSep']
postfix = config['SampleInfo']['Postfix']
core_num = config['GlobalParas']['CpuCore']
meta_name = config['SampleInfo']['MetaData']
# path
fastqc_dir = config['OutputDir']['QC']['fastqc']
fastp_dir = config['OutputDir']['QC']['fastp']
rmrRNA_dir = config['OutputDir']['QC']['rmrRNA']
now_time = config['RunningTime']
rRNA_ref = config['ResourcePath']['rRNARef']  
db_path = config['ResourcePath']['rRNADBPath']

sys.path.append(script_path)
import inputBuild

rule run_fastqc:
    input:
        inputBuild.QC_output(config)  

# fastqc
if end == 'pair':  
    rule fastqc:
        input:
            "Data/{{sample}}_{{end_sep}}.{postfix}".format(postfix=postfix)
        output:
            os.path.join(fastqc_dir, "{sample}_{end_sep}_fastqc.html")
        message:
            " {wildcards.end_sep}: fastqc starting running!"
        params:
            output_path = {fastqc_dir}
        shell:
            "fastqc --extract -t {core_num} -o {params.output_path} {input} 2>> Log/QC_fastqc_{now_time}.log"
            
    rule summaryReport:
        input:
            expand("Result/QC/Fastqc/{sample}_{end_sep}_fastqc.html", sample=samples, end_sep=end_seps)
        output:
            os.path.join(fastqc_dir, "fastqc_report.html")
        message:
            "multiqc starting running!"          
        params:
            path = {fastqc_dir}
        shell:
            "multiqc {params.path} --filename {output} 2>> Log/QC_fastqc_{now_time}.log"

    rule summaryFastqc:
        input:
            os.path.join(fastqc_dir, "fastqc_report.html")
        output:
            os.path.join(fastqc_dir, "fastqc_result_summary.xls")
        message:
            "fastqc result summary!"            
        params:
            path = {fastqc_dir},
            group_map = os.path.join(proj_path, meta_name) 
        shell:
            "Rscript {script_path}/Rcode/fastqc_result_summary.R "
            "{params.group_map} {params.path} {output} {end} '{end_seps}' "
            " 2>> Log/QC_fastqc_{now_time}.log"            
 
else:            
    rule fastqc:
        input:
            "Data/{{sample}}.{postfix}".format(postfix=postfix)
        output:
            os.path.join(fastqc_dir, "{sample}_fastqc.html")
        message:
            "{wildcards.sample}: fastqc starting running!"
        params:
            output_path = {fastqc_dir}
        shell:
            "fastqc --extract -t {core_num} -o {params.output_path} {input} 2>> Log/QC_fastqc_{now_time}.log"
            
    rule summaryReport:
        input:
            expand("Result/QC/Fastqc/{sample}_fastqc.html", sample=samples)
        output:
            os.path.join(fastqc_dir, "fastqc_report.html")
        message:
            "multiqc starting running!"            
        params:
            path = {fastqc_dir}          
        shell:
            "multiqc {params.path} --filename {output} 2>> Log/QC_fastqc_{now_time}.log"

    rule summaryFastqc:
        input:
            os.path.join(fastqc_dir, "fastqc_report.html")
        output:
            os.path.join(fastqc_dir, "fastqc_result_summary.xls")
        message:
            "fastqc result summary!"            
        params:
            path = {fastqc_dir},
            group_map = os.path.join(proj_path, meta_name)
        shell:
            "Rscript {script_path}/Rcode/fastqc_result_summary.R "
            "{params.group_map} {params.path} {output} {end} 'none' "
            " 2>> Log/QC_fastqc_{now_time}.log"

# remove rRNA       
rule build_rRNA_db:
    input:
        os.path.join(rRNA_ref, 'rRNA.fasta')        
    output:
        db = os.path.join(rRNA_ref, 'rRNA', 'rRNA')+'.1.bt2'
    message:
        "Starting building the {organ} rRNA database!"
    params:
        os.path.join(rRNA_ref, 'rRNA', 'rRNA')
    shell:
        "bowtie2-build {input} {params} "
if end == 'pair': 
    end_R1 = end_seps[0]
    end_R2 = end_seps[1]

if end == 'pair':  
    rule remove_rRNA_pair:
        input:
            R1 = "Data/{{sample}}_{}.{}".format(end_R1, postfix), 
            R2 = "Data/{{sample}}_{}.{}".format(end_R2, postfix),
            db = os.path.join(rRNA_ref, 'rRNA', 'rRNA')+'.1.bt2'
        output:
            oR1 = os.path.join(rmrRNA_dir, '{sample}_rRNA_removed_R1.fq.gz'),
            oR2 = os.path.join(rmrRNA_dir, '{sample}_rRNA_removed_R2.fq.gz'),
            temp = temp(os.path.join(rmrRNA_dir, '{sample}_temp.txt'))             
        shell:
            "sh {script_path}/remove_rRNA.sh {db_path} {core_num} {input.R1} {input.R2} "
            "{rmrRNA_dir} {wildcards.sample} {output.oR1} {output.oR2} > {output.temp}"   
else:            
    rule remove_rRNA_single:
        input:
            data = "Data/{{sample}}.{}".format(postfix),
            db = os.path.join(rRNA_ref, 'rRNA', 'rRNA')+'.1.bt2'          
        output:
            out = os.path.join(rmrRNA_dir, '{sample}_rRNA_removed.fq.gz'),
            temp = temp(os.path.join(rmrRNA_dir, '{sample}_temp.txt'))            
        message:
            "{wildcards.sample}: fastp starting running!"
        shell:
            "bowtie2 -x {db_path} -U {input.data} -p {core_num} "
            "--un-gz {output.out} > {output.temp}"
               
# fastp
if end == 'pair':  
    rule fastp:
        input:
            R1 = os.path.join(rmrRNA_dir, '{sample}_rRNA_removed_R1.fq.gz'), 
            R2 = os.path.join(rmrRNA_dir, '{sample}_rRNA_removed_R2.fq.gz')
        output:
            or1 = os.path.join(fastp_dir, "{sample}_R1_clean.fq.gz"), 
            or2 = os.path.join(fastp_dir, "{sample}_R2_clean.fq.gz"), 
            report = os.path.join(fastp_dir, "HTML_Report/{sample}_fastp_report.html")
        message:
            "{wildcards.sample} : fastp starting running!"
        shell:
            "fastp -w {core_num} -i {input.R1} -I {input.R2} -h {output.report} "
            "-o {output.or1} -O {output.or2} 2>> Log/QC_fastp_{now_time}.log"
else:            
    rule fastp:
        input:
            os.path.join(rmrRNA_dir, '{sample}_rRNA_removed.fq.gz')         
        output:
            r1 = os.path.join(fastp_dir, "{sample}_clean.fq.gz"), 
            report = os.path.join(fastp_dir, "HTML_Report/{sample}_fastp_report.html")
        message:
            "{wildcards.sample}: fastp starting running!"
        shell:
            "fastp -w {core_num} -i {input} -h {output.report} "
            "-o {output.r1} 2>> Log/QC_fastp_{now_time}.log"

rule fastp_summary:
    input:
        expand("Result/QC/Fastp/HTML_Report/{sample}_fastp_report.html", sample=samples)       
    output:
        os.path.join(fastp_dir, "fastp_summary_report.xls")
    message:
        "fastp results summary starting running!"
    params:
        path = os.path.join(fastp_dir, "HTML_Report"),
        group_map = os.path.join(proj_path, meta_name)          
    shell:
        "python {script_path}/fastp_summary.py {params.group_map} {params.path} "
        "{output} 2>> Log/QC_fastp_{now_time}.log"

