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

rule run_report:
    input:
        inputBuild.report_output(config)  

rule report:
    input:
        inputBuild.build_input_running(config)
    output:
        inputBuild.report_output(config)
    message:
        "Starting generating analysis report!"
    shell:
        "python {script_path}/reportGenerate.py {proj_path} {output} "


