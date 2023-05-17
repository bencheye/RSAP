
import os
import sys
import time
import ruamel.yaml
from ruamel import yaml

if(len(sys.argv)!=2):
	print("Please input one parameters!")
	print("Please input the parameter: the workspce directory!")
	print("The workspace directory is used to place all the projcet analysis dircetories")  
	exit()

workspace_dir = sys.argv[1]
global_yaml = '/RSAP/global_env.yaml'
Resource_dir = '/RSAP/Resource/Ref_genome/'
aligDBLink = '/RSAP/aligment_database_download_link.xls'

# function
def getYamlConfig(global_yaml):
	with open(global_yaml, 'r', encoding='utf-8') as file:
		yaml_str = file.read()
		yal = ruamel.yaml.YAML()  
		yaml_dict = yal.load(yaml_str)
	return(yaml_dict)

# params
yaml_dict = getYamlConfig(global_yaml)
yaml_dict['workPath'] = workspace_dir

if not os.path.isdir(workspace_dir):  
    os.makedirs(workspace_dir)

if not os.path.isdir(Resource_dir):  
    os.makedirs(Resource_dir)

os.system('cp {} {}'.format(aligDBLink, Resource_dir))

# writte yaml
with open(global_yaml, 'w', encoding='utf-8') as file:
    yaml.dump(yaml_dict, file, Dumper=yaml.RoundTripDumper)