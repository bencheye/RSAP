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
samples = config['SampleInfo']['Sample']
group_compare_combined = ' '.join(config['GroupCompare']['Group'])
now_time = now_time = config['RunningTime']
GEM_dir = config['OutputDir']['Count']
PCA_dir = config['OutputDir']['GroupCompare']['PCA']
NMDS_dir = config['OutputDir']['GroupCompare']['NMDS']

sys.path.append(script_path)
import inputBuild

rule run_GroupCompare:
    input:
        inputBuild.GroupCompare_output(config)

rule PCA_analysis:
    input:
        Rdata = os.path.join(GEM_dir, 'geneExpressionMatrix.Rdata')
    output:
        pcadata = os.path.join(PCA_dir, 'PCA_data.xls'),
        PCA_pdf = os.path.join(PCA_dir, 'PCA_plot.pdf'),
        PCA_png = os.path.join(PCA_dir, 'PCA_plot.png')
    message:
        "Starting PCA analysis!"
    shell:
        "Rscript {script_path}/Rcode/PCA_code.R {input} '{group_compare_combined}' "
        "{output.pcadata} {output.PCA_pdf} {output.PCA_png} "
        " 2>> Log/GroupCompare_PCA_{now_time}.log"

rule MDS_analysis:
    input:
        Rdata = os.path.join(GEM_dir, 'geneExpressionMatrix.Rdata')
    output:
        NMDSdata = os.path.join(NMDS_dir, 'NMDS_data.xls'),
        NMDS_pdf = os.path.join(NMDS_dir, 'NMDS_plot.pdf'),
        NMDS_png = os.path.join(NMDS_dir, 'NMDS_plot.png')
    message:
        "Starting NMDS analysis!"
    shell:
        "Rscript {script_path}/Rcode/MDS_code.R {input} '{group_compare_combined}' "
        "{output.NMDSdata} {output.NMDS_pdf} {output.NMDS_png} "
        " 2>> Log/GroupCompare_NMDS_{now_time}.log"