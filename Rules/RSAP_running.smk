import os
import sys
import time
import ruamel.yaml
import pandas as pd
from ruamel import yaml
from multiprocessing import cpu_count

# Input
configfile: "project_run_config.yaml"
# path
script_path = config['ResourcePath']['ScriptPath']
rules_path = config['ResourcePath']['RulesPath']
proj_path = config['ResourcePath']['projPath']
ref_path = config['ResourcePath']['refPath']
rRNA_ref = config['ResourcePath']['rRNARef']  
db_path = config['ResourcePath']['rRNADBPath']
# global parameter
gtf_link = config['ResourcePath']['gtfLink']
method_link = config['ResourcePath']['methodLink']
organ = config['GlobalParas']['Organ']
core_num = config['GlobalParas']['CpuCore']
align_method = config['GlobalParas']['aligMethod']
now_time = now_time = config['RunningTime']
reserved_item = config['GlobalParas']['reservedItem']
DEA_method = config['GlobalParas']['DEAMethod']
# sample
samples = config['SampleInfo']['Sample']
meta_name = config['SampleInfo']['MetaData']
group_compare_combined = ' '.join(config['GroupCompare']['Group'])
end = config['SampleInfo']['END']
end_seps = config['SampleInfo']['PairedSep']
postfix = config['SampleInfo']['Postfix']
# output dir
fastqc_dir = config['OutputDir']['QC']['fastqc']
fastp_dir = config['OutputDir']['QC']['fastp']
rmrRNA_dir = config['OutputDir']['QC']['rmrRNA']
alignment_dir = config['OutputDir']['ALIGNMENT']
GEM_dir = config['OutputDir']['GEM']
PCA_dir = config['OutputDir']['GVir']['PCA']
NMDS_dir = config['OutputDir']['GVir']['NMDS']
DEA_dir = config['OutputDir']['DEA']['DEA']
volcanoPlot_dir = config['OutputDir']['DEA']['volcanoPlot']
heatmap_dir = config['OutputDir']['DEA']['Heatmap']
ERKEGG_dir = config['OutputDir']['ENRICHMENT']['KEGG']
ERGO_dir = config['OutputDir']['ENRICHMENT']['GO']

sys.path.append(script_path)
import inputBuild
rule onekey:
    input:
        inputBuild.build_input_all(config)

for smk_file in os.listdir(rules_path):
    include: os.path.join(rules_path, smk_file)




