import os
import sys
import ruamel.yaml
from ruamel import yaml

def getYamlConfig(yaml_str):
    with open(yaml_str, 'r', encoding='utf-8') as file:
        yaml_str = file.read()
        yal = ruamel.yaml.YAML()  
        yaml_dict = yal.load(yaml_str)
    return(yaml_dict)

def get_figure_text(img_dir, img_title):
    out_text = '<div align="center"> FigureNum {}</div>\n'.format(img_title)
    out_text = out_text + '![{}]({})'.format(img_title, img_dir)
    return out_text

def headReport(config, report_temp):
    project_name = config['PROJECT'].split('_RNAseq_')[0]
    compare_text = buildGroupCompareText(config)
    compare_text = ' '.join(compare_text)
    running_text = buildRunningModuleText(config)
    head_dir = os.path.join(report_temp, 'head.md')
    with open(head_dir, 'r') as f:
        head_text = f.read()
    out_text = head_text.format(
            project_name, 
            config['RunningTime'], 
            config['SampleInfo']['END'],
            config['GlobalParas']['Organ'],
            config['GlobalParas']['referenceType'],
            config['GlobalParas']['aligMethod'],
            config['GlobalParas']['DEAMethod'],
            config['GlobalParas']['reservedItem'],
            running_text,
            compare_text
            )
    return out_text

def qcReport(config, report_temp):
    fastqcMD_dir = 'Report/Table/fastqc_result_summary.xls'
    fastqc_table = transformTable2md(fastqcMD_dir, 
                                       'Quality statistics of raw sequence data')
    fastpMD_dir = 'Report/Table/fastp_summary_report.xls'
    fastp_table = transformTable2md(fastpMD_dir, 
                                       'Quality statistics after quality control')
    head_dir = os.path.join(report_temp, 'QC.md')
    with open(head_dir, 'r') as f:
        head_text = f.read()
    out_text = head_text.format(fastqc_table, fastp_table)
    return out_text

def alignmentReport(config, report_temp):
    alignment_method = config['GlobalParas']['aligMethod']
    referenceDB = config['GlobalParas']['referenceType']
    head_dir = os.path.join(report_temp, 'Alignment.md')
    with open(head_dir, 'r') as f:
        head_text = f.read()
    out_text = head_text.format(alignment_method, referenceDB)
    return out_text

def countReport(config, report_temp):
    reserved_item = config['GlobalParas']['reservedItem']
    countMat_dir = 'Report/Table/geneExpressMatrix.xls'
    countMat_table = transformTable2md(countMat_dir, 
                                       'Example of gene expression matrix results')
    head_dir = os.path.join(report_temp, 'Count.md')
    with open(head_dir, 'r') as f:
        head_text = f.read()
    out_text = head_text.format(reserved_item, countMat_table)
    return out_text

def gVirReport(config, report_temp):
    pca_dir = 'Figure/PCA_plot.png'
    pca_fig = get_figure_text(pca_dir, 'PCA result diagram')
    nmds_dir = 'Figure/NMDS_plot.png'
    nmds_fig = get_figure_text(nmds_dir, 'NMDS result diagram')    
    head_dir = os.path.join(report_temp, 'GVir.md')
    with open(head_dir, 'r') as f:
        head_text = f.read()
    out_text = head_text.format(pca_fig, nmds_fig)
    return out_text

def deaReport(config, report_temp):
    DEA_method = config['GlobalParas']['DEAMethod']
    compare_text = buildGroupCompareText(config)
    degRes_text = ''
    heatmapRes_text = ''
    volcanoPlotRes_text = ''
    for item in compare_text:
        degMat_dir = 'Report/Table/DEG_{}.xls'.format(item)
        deg_title = 'Example of differentially expressed analysis result of the {}'.format(item)
        degMat_table = transformTable2md(degMat_dir, deg_title)
        degRes_text = degRes_text + degMat_table
        heatmap_dir = 'Figure/heatmap_{}.png'.format(item)
        heatmap_title = 'Heatmap of differentially expressed analysis result of the {}'.format(item)
        heatmap_fig = get_figure_text(heatmap_dir, heatmap_title)
        heatmapRes_text = heatmapRes_text + heatmap_fig
        volcanoPlot_dir = 'Figure/volcanoPlot_{}.png'.format(item)
        volcanoPlot_title = 'VolcanoPlot of differentially expressed analysis result of the {}'.format(item)
        volcanoPlot_fig = get_figure_text(volcanoPlot_dir, volcanoPlot_title)
        volcanoPlotRes_text = volcanoPlotRes_text + volcanoPlot_fig
    head_dir = os.path.join(report_temp, 'DEA.md')
    with open(head_dir, 'r') as f:
        head_text = f.read()
    out_text = head_text.format(DEA_method, degRes_text, heatmapRes_text, volcanoPlotRes_text)
    return out_text

def enrichmentReport(config, report_temp):
    compare_text = buildGroupCompareText(config)
    keggRes_text = ''
    goRes_text = ''
    for item in compare_text:
        kegg_dir = 'Figure/DEG_KEGG_enrichment_{}.png'.format(item)
        kegg_title = 'Enrichment analysis result of KEGG of the {}'.format(item)
        kegg_fig = get_figure_text(kegg_dir, kegg_title)
        keggRes_text = keggRes_text + kegg_fig
        go_dir = 'Figure/DEG_GeneOntology_BP_enrichment_{}.png'.format(item)
        go_title = 'Enrichment analysis result of biological process (BP) of GO of the {}'.format(item)
        go_fig = get_figure_text(go_dir, go_title)
        goRes_text = goRes_text + go_fig
    head_dir = os.path.join(report_temp, 'enrichment.md')
    with open(head_dir, 'r') as f:
        head_text = f.read()
    out_text = head_text.format(goRes_text, keggRes_text)
    return out_text

def generateTableNum(report_text, id_text = 'TableNum'):
    table_num = report_text.count(id_text)
    for num in range(1, table_num+1):
        id_num = 'Table' + str(num)
        report_text = report_text.replace(id_text, id_num, 1)
    return report_text

def generateFigNum(report_text, id_text = 'FigureNum'):
    fig_num = report_text.count(id_text)
    for num in range(1, fig_num+1):
        id_num = 'Figure' + str(num)
        report_text = report_text.replace(id_text, id_num, 1)
    return report_text

def generateFinalReport(config, report_temp):
    report_text = headReport(config, report_temp)
    if config['Running']['QC'] == 'yes':
        report_text = report_text + qcReport(config, report_temp)
    if config['Running']['ALIGNMENT'] == 'yes':
        report_text = report_text + alignmentReport(config, report_temp)
    if config['Running']['Count'] == 'yes':
        report_text = report_text + countReport(config, report_temp)
    if config['Running']['GroupCompare'] == 'yes':
        report_text = report_text + gVirReport(config, report_temp)
    if config['Running']['DEA'] == 'yes':
        report_text = report_text + deaReport(config, report_temp)
    if config['Running']['ENRICHMENT'] == 'yes':
        report_text = report_text + enrichmentReport(config, report_temp)
    report_text = generateTableNum(report_text, id_text = 'TableNum')
    report_text = generateFigNum(report_text, id_text = 'FigureNum')
    return report_text

# Input
proj_path = os.path.abspath(sys.argv[1])
md_output = os.path.abspath(sys.argv[2])
yaml_file = os.path.join(proj_path, 'project_run_config.yaml')

# read the config file
config = getYamlConfig(yaml_file)
script_path = config['ResourcePath']['ScriptPath']
sys.path.append(script_path)
from functionLibrary import *

report_temp = config['ResourcePath']['ReportTemplate']
prepareReportData(config)
md_text = generateFinalReport(config, report_temp)
with open(md_output, 'w') as f:
    f.write(md_text)




