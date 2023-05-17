# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 09:18:41 2021
function: generate rnaseq pipeline running config 
@author: benchenye
"""
import os
import sys
import ruamel.yaml
from ruamel import yaml

script_path = sys.argv[1]
proj_name = sys.argv[2]
End = sys.argv[3]
rules_path = sys.argv[4]
proj_path = sys.argv[5]
reference_path = sys.argv[6]
globalYaml_Path = sys.argv[7]
yaml_file = os.path.join(proj_path, 'project_run_config.yaml') 

if End == 'pair':
    end_sep = ['R1', 'R2']
elif End == 'single':
    end_sep = ""

yaml_str = """\
# RSAP running parameters
# Please check the parameters, and adjust them according to your circumstance
# There are some parameters you must

# Project name
PROJECT : {}
RunningTime : 
# ==================== Sample Information =====================
SampleInfo:
    MetaData: "metadata.xls"    
    Sample : 
    Group :
    END :  {}                   
    PairedSep : {}              
    Postfix: 'fastq.gz'        

# ==================== Resource Path ==========================
## DO NOT change these parameters
ResourcePath:
    ScriptPath : {}
    ReportPath: Report
    ReportTemplate: /RSAP/ReportTemplate/md
    RulesPath : {}
    projPath : {}
    referencePath : {} 
    globalYamlPath : {}
    aligDBLink : "aligment_database_download_link.xls"
    rRNARef : 
    rRNADBPath :
    gtfLink : 
    methodLink :  

# ==================== Global parameters ======================
GlobalParas :
    Organ : homo_sapiens           
    referenceType : genome
    CpuCore : 4                   
    aligMethod : "hisat2"          
    DEAMethod : "edgeR"            
    reservedItem: 'protein_coding'  

# ====================== Group compare ========================
GroupCompare:

# ================== Control of the workflow ==================
## just usage the value of "yes" or "no"
## 'yes' means running this analysis, and 'no' means not
Running:
    QC: yes  
    ALIGNMENT: yes  
    Count: yes  
    GroupCompare: yes 
    DEA: yes   
    ENRICHMENT: yes  
OutputDir:
    QC:
        fastqc : 'Result/QC/Fastqc'
        fastp : 'Result/QC/Fastp'
        rmrRNA : 'Result/QC/rm_rRNA'
    ALIGNMENT: 
    Count: Result/GeneExpMatrix
    GroupCompare:
        PCA : 'Result/GVir/PCA'
        NMDS : 'Result/GVir/NMDS' 
    DEA: 
        DEA : 
        volcanoPlot : 'Result/DEA/Volcanoplot' 
        Heatmap : 'Result/DEA/Heatmap'
    ENRICHMENT: 
        KEGG : 'Result/Enrichment/KEGG'
        GO : 'Result/Enrichment/GO'
        
""".format(proj_name, End, end_sep, script_path, rules_path, 
    proj_path, reference_path, globalYaml_Path)

# load yaml    
yal = ruamel.yaml.YAML()  
yaml_dict = yal.load(yaml_str)
# writte yaml
with open(yaml_file, 'w', encoding='utf-8') as file:
    yaml.dump(yaml_dict, file, Dumper=yaml.RoundTripDumper)

