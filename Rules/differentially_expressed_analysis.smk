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
DEA_method = config['GlobalParas']['DEAMethod']
now_time = now_time = config['RunningTime']
Count_dir = config['OutputDir']['Count']
DEA_dir = config['OutputDir']['DEA']['DEA']
volcanoPlot_dir = config['OutputDir']['DEA']['volcanoPlot']
heatmap_dir = config['OutputDir']['DEA']['Heatmap']

sys.path.append(script_path)
import inputBuild

rule run_DEA:
    input:
        inputBuild.DEA_output(config)

if DEA_method == 'edgeR':
    rule run_edgeR_DEA:
        input:
            Rdata = os.path.join(Count_dir, 'geneExpressionMatrix.Rdata')
        output:
            DEA_resFile = os.path.join(DEA_dir, 'DEA_{group_compare}_result.xls'),
            deg_file =  os.path.join(DEA_dir, 'DEG_{group_compare}_result.xls')
        message:
            "Starting DEA-edgeR analysis!"
        shell:
            "Rscript {script_path}/Rcode/edgeR_DEG_analysis.R {input.Rdata} "
            "'{wildcards.group_compare}' "
            "{output.DEA_resFile} {output.deg_file} "
            "2>> Log/DEA_edgeR_{now_time}.log"
elif DEA_method == 'DESeq2':
    rule run_DESeq2_DEA:
        input:
            Rdata = os.path.join(Count_dir, 'geneExpressionMatrix.Rdata')
        output:
            DEA_resFile = os.path.join(DEA_dir, 'DEA_{group_compare}_result.xls'),
            deg_file =  os.path.join(DEA_dir, 'DEG_{group_compare}_result.xls')
        message:
            "Starting DEA-edgeR analysis!"
        shell:
            "Rscript {script_path}/Rcode/DESeq2.R {input.Rdata} "
            "'{wildcards.group_compare}' "
            "{output.DEA_resFile} {output.deg_file} "
            "2>> log/DEA_DESeq2_{now_time}.log"

rule run_VolcanoPlot:
    input:
        Rdata = os.path.join(Count_dir, 'geneExpressionMatrix.Rdata'),
        deg_file =  os.path.join(DEA_dir, 'DEG_{group_compare}_result.xls')
    output:
        VolcanoPlot_data = os.path.join(volcanoPlot_dir, 'VolcanoPlot_{group_compare}_result.xls'),
        pdf = os.path.join(volcanoPlot_dir, 'Volcanoplot_{group_compare}.pdf'),
        png = os.path.join(volcanoPlot_dir, 'Volcanoplot_{group_compare}.png')
    message:
        "Starting DEA-edgeR analysis!"
    shell:
        "Rscript {script_path}/Rcode/VolcanoPlot.R {input.Rdata} "
        "'{wildcards.group_compare}' "
        "{output.VolcanoPlot_data} {output.pdf} {output.png} "
        "2>> Log/DEA_volcanoPlot_{now_time}.log"

rule Heatmap_analysis:
    input:
        Rdata = os.path.join(Count_dir, 'geneExpressionMatrix.Rdata'),
        deg_file =  os.path.join(DEA_dir, 'DEG_{group_compare}_result.xls')
    output:
        pdf = os.path.join(heatmap_dir, 'HeatmapPlot_{group_compare}.pdf')
    message:
        "Starting Heatmap analysis!"
    shell:
        "Rscript {script_path}/Rcode/heatmapPlot.R {input.Rdata} "
        "'{wildcards.group_compare}' {output.pdf} "
        "2>> Log/DEA_heatmap_{now_time}.log"

rule Heatmap_pdfToPng:
    input:
        pdf = os.path.join(heatmap_dir, 'HeatmapPlot_{group_compare}.pdf')
    output:
        png = os.path.join(heatmap_dir, 'HeatmapPlot_{group_compare}.png')
    message:
        "Starting transform pdf into png of Heatmap analysis!"
    shell:
        "python {script_path}/pdfToPng.py {input.pdf} {output.png}"