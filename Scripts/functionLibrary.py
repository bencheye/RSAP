import pandas as pd
import os

# get aligment database download link
def getAligmentDBLink(align_method, resource_path, 
    alig_DB_link, organ, reference_type):
    DB_link_mappingFile = os.path.join(resource_path, alig_DB_link)
    DB_link = pd.read_csv(DB_link_mappingFile, sep = '\t')     
    method_row = (DB_link['Organ'] == organ) & (DB_link['Method'] == align_method) & (DB_link['Type'] == reference_type)
    gtf_row = (DB_link['Organ'] == organ) & (DB_link['Method'] == 'gtf')
    gtf_link = DB_link.loc[gtf_row, 'Link'].tolist()[0]
    method_link = DB_link.loc[method_row, 'Link'].tolist()[0]
    return gtf_link, method_link

def get_group_compareList(group_compares):
    items=group_compares.items()
    group_comparesList = []
    for item in items:
        group = item[0]
        groupCompare = item[1]
        if group != 'Group':
            group_comparesList.append('_vs_'.join(groupCompare))  
    return(group_comparesList) 
    
def buildGroupCompareText(config):
    group_compareList = config['GroupCompare']
    out_list = []
    for item in group_compareList:
        if item != 'Group':
            group_text = '_vs_'.join(group_compareList[item])
            out_list.append(group_text)
    return out_list

def buildRunningModuleText(config):
    running_dict = config['Running']
    out_list = []
    for item in running_dict:
        if running_dict[item] == 'yes':
            out_list.append(item) 
    out_text = ' '.join(out_list)  
    return out_text 

def transformTable2md(input_file, title_text):
    with open(input_file, 'r') as f:  
        lines = f.readlines() 
        count = 0 
        out_text = '<div align="center"> TableNum ' + title_text + '</div>\n'
        for line in lines:
            rec = line.split('\t')
            new_rec = [i.strip() for i in rec]
            new_rec = '| '+ ' | '.join(new_rec) +' |'+'\n'
            out_text = out_text + new_rec
            #md.write(new_rec)
            if count == 0:
                head_rec = '|'
                for i in rec:
                    head_rec += ' --- |'
                head_rec += '\n'
                # md.write(head_rec)
                out_text = out_text + head_rec
            count = count + 1
    return(out_text)

def buildTableRow(item_list):
    out_text = ''
    for item in item_list:
        one_item = '<td>' + item + '</td>'
        out_text = out_text + one_item
    return(out_text)
     
def transformTable2html(input_file, title_text):
    with open(input_file, 'r') as f:  
        lines = f.readlines()  
        out_text = '<div align="center"> TableNum ' + title_text + '</div><br>\n' 
        out_text = out_text + '<table width="100%">'
        for line in lines:
            rec = line.split('\t')
            new_rec = [i.strip() for i in rec]
            new_rec = '<tr>'+ buildTableRow(new_rec)+'</tr>'
            out_text = out_text + new_rec
    return(out_text) 


# prepare data
def prepareReportData(config):
    report_dir = config['ResourcePath']['ReportPath']
    table_dir = os.path.join(report_dir, 'Table')
    figure_dir = os.path.join(report_dir, 'Figure')
    compare_text = buildGroupCompareText(config)
    if not os.path.isdir(table_dir):  
        os.makedirs(table_dir)
    if not os.path.isdir(figure_dir):  
        os.makedirs(figure_dir)
    if config['Running']['QC'] == 'yes':
        os.system('cp {}/fastqc_result_summary.xls {}/'.format(
                config['OutputDir']['QC']['fastqc'], table_dir))
        os.system('cp {}/fastp_summary_report.xls {}/'.format(
                config['OutputDir']['QC']['fastp'], table_dir)) 
    if config['Running']['Count'] == 'yes':
        os.system("head -n 11 {}/gene_expMatrix_final.csv | tr ',' '\t' > {}/geneExpressMatrix.xls".format(
                config['OutputDir']['Count'], table_dir))
    if config['Running']['DEA'] == 'yes':
        for item in compare_text:
            deg_name = "DEG_{}_result.xls".format(item)
            heatmap_name = "HeatmapPlot_{}.png".format(item)
            volcanoPlot_name = "Volcanoplot_{}.png".format(item)
            os.system("head -n 11 {}/{} > {}/DEG_{}.xls".format(
                    config['OutputDir']['DEA']['DEA'], deg_name, table_dir, item))
            os.system("cp {}/{} {}/heatmap_{}.png".format(
                    config['OutputDir']['DEA']['Heatmap'], heatmap_name, figure_dir, item))
            os.system("cp {}/{} {}/volcanoPlot_{}.png".format(
                    config['OutputDir']['DEA']['volcanoPlot'], volcanoPlot_name, figure_dir, item))
    if config['Running']['GroupCompare'] == 'yes':
        os.system("cp {}/PCA_plot.png {}/".format(
                config['OutputDir']['GroupCompare']['PCA'], figure_dir))
        os.system("cp {}/NMDS_plot.png {}/".format(
                config['OutputDir']['GroupCompare']['NMDS'], figure_dir))
    if config['Running']['ENRICHMENT'] == 'yes':
        for item in compare_text:
            kegg_name = "DEG_KEGG_enrichment_{}.png".format(item)
            go_name = "DEG_GeneOntology_BP_enrichment_{}.png".format(item)
            os.system("cp {}/{} {}/DEG_KEGG_enrichment_{}.png".format(
                    config['OutputDir']['ENRICHMENT']['KEGG'], kegg_name, figure_dir, item))
            os.system("cp {}/{} {}/DEG_GeneOntology_BP_enrichment_{}.png".format(
                    config['OutputDir']['ENRICHMENT']['GO'], go_name, figure_dir, item))

