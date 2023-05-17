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

# Input
proj_path = os.path.abspath(sys.argv[1])
yaml_file = os.path.join(proj_path, 'project_run_config.yaml')
now_time = time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime())

def unique(data):
    new_list = []
    for item in data:
        if item not in new_list:
            new_list.append(item)
    return new_list

def getYamlConfig(yaml_str):
    with open(yaml_str, 'r', encoding='utf-8') as file:
        yaml_str = file.read()
        yal = ruamel.yaml.YAML()  
        yaml_dict = yal.load(yaml_str)
    return(yaml_dict)

# read the config file
yaml_dict = getYamlConfig(yaml_file)

# read metaData to get sample and group information
meta_name = yaml_dict['SampleInfo']['MetaData']
try:
    meta_data = pd.read_csv(meta_name, sep = '\t')
except FileNotFoundError:
    print("File is not found!")    
sample = meta_data['Sample'].tolist()
group = meta_data['Group'].tolist()

# load parameters
resource_path = yaml_dict['ResourcePath']['referencePath']
reference_type = yaml_dict['GlobalParas']['referenceType']
organ = yaml_dict['GlobalParas']['Organ']
align_method = yaml_dict['GlobalParas']['aligMethod']
script_path = yaml_dict['ResourcePath']['ScriptPath']
alig_DB_link = yaml_dict['ResourcePath']['aligDBLink']
DEA_method = yaml_dict['GlobalParas']['DEAMethod']

# preprocess
sys.path.append(script_path)
import functionLibrary
rRNA_ref = os.path.join(resource_path, organ)  
db_path = os.path.join(rRNA_ref, 'rRNA', 'rRNA')
ref_path = os.path.join(resource_path, organ, reference_type)
gtf_link, method_link = functionLibrary.getAligmentDBLink(align_method, resource_path,
 alig_DB_link, organ, reference_type)
DEA_dir = os.path.join('Result/DEA', DEA_method)

# update sample and group info
print("Fill parameters' Sample and Group!")
yaml_dict['SampleInfo']['Sample'] = sample
yaml_dict['SampleInfo']['Group'] = group
yaml_dict['RunningTime'] = now_time
yaml_dict['ResourcePath']['refPath'] = ref_path
yaml_dict['ResourcePath']['rRNARef'] = rRNA_ref
yaml_dict['ResourcePath']['rRNADBPath'] = db_path
yaml_dict['OutputDir']['ALIGNMENT'] = os.path.join('Result/Alignment', align_method)
yaml_dict['ResourcePath']['gtfLink'] = gtf_link
yaml_dict['ResourcePath']['methodLink'] = method_link
yaml_dict['OutputDir']['DEA']['DEA'] = DEA_dir

# update group_compare
print("Fill parameters' group_compare!")
uniq_group = unique(group)
if yaml_dict['GroupCompare'] == None:
    yaml_dict.update({"GroupCompare": {"Group":uniq_group}}) 
    num = 1
    for i in range(0, len(uniq_group)):
        for j in range(i+1, len(uniq_group)):
            #group_compare = group[i] + ' ' + group[j]
            group_compare = [uniq_group[i], uniq_group[j]]
            group_name = 'group' + str(num)
            num = num + 1
            yaml_dict['GroupCompare'][group_name] = group_compare

# writte yaml
with open(yaml_file, 'w', encoding='utf-8') as file:
    yaml.dump(yaml_dict, file, Dumper=yaml.RoundTripDumper)

## check input and parameter error, before running
# load parameters
End = yaml_dict['SampleInfo']['END']
CpuCore = yaml_dict['GlobalParas']['CpuCore']
Organ = yaml_dict['GlobalParas']['Organ']
alig_method = yaml_dict['GlobalParas']['aligMethod']
DEA_method = yaml_dict['GlobalParas']['DEAMethod']
pair_end = yaml_dict['SampleInfo']['PairedSep']
Postfix = yaml_dict['SampleInfo']['Postfix']

print("Check running parameters!")
# check the number of cpucore whether is greater than computer's
if CpuCore > cpu_count():
    print('Please change the CpuCore value, it is greater than the cpucore of computer!')
    sys.exit(1)
else:
    print("cpu core check pass!")
# check sample
print("Just check input sample when starting the whole pipeline!")
if yaml_dict['Running']['QC'] == 'yes':
    for num in range(len(sample)):
        if End == 'pair' :
            R1 = os.path.join(proj_path, 'Data/{}_{}.{}'.format(sample[num], pair_end[0], Postfix))
            R2 = os.path.join(proj_path, 'Data/{}_{}.{}'.format(sample[num], pair_end[1], Postfix))
            if not os.path.exists(R1):
                print('There is not the file: {} !'.format(R1))
                sys.exit(1)
            if not os.path.exists(R2):
                print('There is not the file: {} !'.format(R2))
                sys.exit(1)
        elif End == 'single' :
            fq_file = os.path.join(proj_path, 'Data/{}.{}'.format(sample[num], Postfix))
            if not os.path.exists(fq_file) :
                print('There is not the file: {} !'.format(fq_file))
                sys.exit(1)
        else:
            print('sample check pass!')

# load global yaml
global_yaml = yaml_dict['ResourcePath']['globalYamlPath']
globalYaml_dict = getYamlConfig(global_yaml)
aligMethod_list = globalYaml_dict['aligMethod']
DEAMethod_list = globalYaml_dict['DEAMethod']
Organ_list = globalYaml_dict['Organ']
# check organ
if (Organ in Organ_list):
    print ("organ check pass")
    print("Please reconfirm the analysis organ whether is {}!".format(Organ))
    time.sleep(1)
else:
    print('The organ: {} you filled could not find in the "global_env.yaml"!'.format(Organ))
    sys.exit(1)

if (alig_method in aligMethod_list):
    print ("alignment method check pass!")
else:
    print('The alignment method: {} you filled could not find in the "global_env.yaml"!'.format(alig_method))
    sys.exit(1)

if (DEA_method in DEAMethod_list):
    print ("DEA method check pass!")
else:
    print('The DEA method: {} you filled could not find in the "global_env.yaml"!'.format(DEA_method))
    sys.exit(1)

