# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 15:59:04 2021

@author: bench
"""
import os
import sys

# 
configfile: "project_run_config.yaml"
proj_path = config['ResourcePath']['projPath']
meta_name = config['SampleInfo']['MetaData']
script_path = config['ResourcePath']['ScriptPath']
now_time = config['RunningTime']
Count_dir = config['OutputDir']['Count']

sys.path.append(script_path)
import inputBuild

rule run_build_Rdata:
    input:
        inputBuild.build_Rdata_output(config)

rule build_downstreamAnalysis:
    input:
        exp_matrix = os.path.join(Count_dir, 'gene_expMatrix_final.csv'),
        group_map = os.path.join(proj_path, meta_name)
    output:
        Rdata = os.path.join(Count_dir, 'geneExpressionMatrix.Rdata')
    message:
        "Build the Rdata of downstream analysis!"
    shell:
        "Rscript {script_path}/Rcode/build_downstream_Rdata.R "
        "{input.exp_matrix} {input.group_map} {output.Rdata} "


