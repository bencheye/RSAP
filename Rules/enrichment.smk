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
now_time = now_time = config['RunningTime']
Count_dir = config['OutputDir']['Count']
DEA_dir = config['OutputDir']['DEA']['DEA']
ERKEGG_dir = config['OutputDir']['ENRICHMENT']['KEGG']
ERGO_dir = config['OutputDir']['ENRICHMENT']['GO']
organ = config['GlobalParas']['Organ']

sys.path.append(script_path)
import inputBuild

rule run_enrichment:
    input:
        inputBuild.enrichment_output(config)

rule enrichment_analysis:
    input:
        Rdata = os.path.join(Count_dir, 'geneExpressionMatrix.Rdata'),
        deg_file =  os.path.join(DEA_dir, 'DEG_{group_compare}_result.xls')
    output:
        KEGG_pdf = os.path.join(ERKEGG_dir, 'DEG_KEGG_enrichment_{group_compare}.pdf'),
        KEGG_png = os.path.join(ERKEGG_dir, 'DEG_KEGG_enrichment_{group_compare}.png'),
        GO_BP_pdf = os.path.join(ERGO_dir, 'DEG_GeneOntology_BP_enrichment_{group_compare}.pdf'),
        GO_BP_png = os.path.join(ERGO_dir, 'DEG_GeneOntology_BP_enrichment_{group_compare}.png')
    message:
        "Starting enrichment analysis!"
    shell:
        "Rscript {script_path}/Rcode/enrichment.R {input.Rdata} "
        "'{wildcards.group_compare}' {output.GO_BP_pdf} "
        "{output.GO_BP_png} {output.KEGG_pdf} {output.KEGG_png} {organ} "
        "2>> Log/enrichment_{now_time}.log"



