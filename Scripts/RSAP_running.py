# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 09:18:41 2021
function: generate rnaseq pipeline running config 
@author: benchenye
"""

import os
import sys
import time
import ruamel.yaml
import pandas as pd
from ruamel import yaml
from multiprocessing import cpu_count

# function
def getYamlConfig(yaml_str):
    with open(yaml_str, 'r', encoding='utf-8') as file:
        yaml_str = file.read()
        yal = ruamel.yaml.YAML()  
        yaml_dict = yal.load(yaml_str)
    return(yaml_dict)

# Input
yaml_file = 'project_run_config.yaml'
yaml_dict = getYamlConfig(yaml_file)
# load parameters
proj_path = yaml_dict['ResourcePath']['projPath']
script_path = yaml_dict['ResourcePath']['ScriptPath']
rules_path = yaml_dict['ResourcePath']['RulesPath']

print("Starting running check!")
check_res = os.system("python {}/run_initional_check_error.py {}".format(
    script_path, proj_path))
if check_res == 0:
    os.system('snakemake -s {}/RSAP_running.smk'.format(rules_path))
else:
    print('check not pass!')
print("RSAP analysis has finished!")