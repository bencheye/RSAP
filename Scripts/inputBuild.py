# function 
from snakemake.io import expand
import os

def get_group_compareList(group_compares):
    items=group_compares.items()
    group_comparesList = []
    for item in items:
        group = item[0]
        groupCompare = item[1]
        if group != 'Group':
            group_comparesList.append('_vs_'.join(groupCompare))  
    return(group_comparesList) 

# QC input
def QC_output(config):
    samples = config['SampleInfo']['Sample']
    end_seps = config['SampleInfo']['PairedSep']
    end = config['SampleInfo']['END']
    input_list = []
    input_list.append(fastqc_output(config))
    input_list.append(remove_rRNA_output(config))
    input_list.append(fastp_output(config))
    return(input_list)

def fastqc_output(config):
    samples = config['SampleInfo']['Sample']
    end_seps = config['SampleInfo']['PairedSep']
    end = config['SampleInfo']['END']   
    fastqc_dir = config['OutputDir']['QC']['fastqc'] 
    output_items = []
    if end == 'pair':  
        html_report = expand(fastqc_dir + "/{sample}_{end_sep}_fastqc.html", 
            sample=samples, end_sep=end_seps)
        output_items.append(html_report)
    elif end == 'single':
        html_report = expand(fastqc_dir + "/{sample}_fastqc.html", sample=samples)
        output_items.append(html_report)
    # summary_report = os.path.join(fastqc_dir, "fastqc_report.html")
    fastqc_summary = os.path.join(fastqc_dir, "fastqc_result_summary.xls")  
    # output_items.append(summary_report)
    output_items.append(fastqc_summary)  
    return(output_items)

def remove_rRNA_output(config):
    samples = config['SampleInfo']['Sample']
    end = config['SampleInfo']['END'] 
    rmrRNA_dir = config['OutputDir']['QC']['rmrRNA']
    output_items = []
    if end == 'pair':  
        rmrRNA_R1 = expand(rmrRNA_dir + '/{sample}_rRNA_removed_R1.fq.gz', sample=samples)
        rmrRNA_R2 = expand(rmrRNA_dir + '/{sample}_rRNA_removed_R2.fq.gz', sample=samples)
        output_items.append(rmrRNA_R1)
        output_items.append(rmrRNA_R2)
    elif end == 'single':
        rmrRNA = expand(rmrRNA_dir + "/{sample}_rRNA_removed.fq.gz", sample=samples) 
        output_items.append(rmrRNA)
    return(output_items)

def fastp_output(config):
    samples = config['SampleInfo']['Sample']
    end = config['SampleInfo']['END'] 
    fastp_dir = config['OutputDir']['QC']['fastp']
    output_items = []
    if end == 'pair':  
        R1 = expand(fastp_dir + "/{sample}_R1_clean.fq.gz", sample=samples)
        R2 = expand(fastp_dir + "/{sample}_R2_clean.fq.gz", sample=samples)
        output_items.append(R1)
        output_items.append(R2)
    elif end == 'single':
        data = expand(fastp_dir + "/{sample}_clean.fq.gz", sample=samples) 
        output_items.append(data)
    html_report = expand(fastp_dir + "/HTML_Report/{sample}_fastp_report.html", sample=samples)
    fastp_summary = os.path.join(fastp_dir, "fastp_summary_report.xls")
    output_items.append(html_report)
    output_items.append(fastp_summary)
    return(output_items)

def alignment_output(config):
    samples = config['SampleInfo']['Sample']
    alignment_dir = config['OutputDir']['ALIGNMENT']
    output_items = []
    output = expand(alignment_dir + "/{sample}.bam", sample=samples)
    output_items.append(output)
    return(output_items)

def feature_output(config):
    Count_dir = config['OutputDir']['Count']
    output_items = [] 
    output = os.path.join(Count_dir, 'gene_expMatrix_final.csv')
    output_items.append(output)
    return(output_items) 
                    
def build_Rdata_output(config):
    Count_dir = config['OutputDir']['Count']
    output_items = [] 
    output = os.path.join(Count_dir, 'geneExpressionMatrix.Rdata')
    output_items.append(output)
    return(output_items)

def DEA_output(config):
    DEA_method = config['GlobalParas']['DEAMethod']
    group_compares = config['GroupCompare']
    group_comparesList = get_group_compareList(group_compares)
    DEA_dir = config['OutputDir']['DEA']['DEA']
    volcanoPlot_dir = config['OutputDir']['DEA']['volcanoPlot']
    heatmap_dir = config['OutputDir']['DEA']['Heatmap']
    output_items = []
    de_res = expand(DEA_dir + '/DEA_{group_compare}_result.xls', group_compare = group_comparesList)
    deg = expand(DEA_dir + '/DEG_{group_compare}_result.xls', group_compare = group_comparesList)
    volcanoplotData = expand(volcanoPlot_dir + '/VolcanoPlot_{group_compare}_result.xls',
     group_compare = group_comparesList)
    pdf = expand(volcanoPlot_dir + '/Volcanoplot_{group_compare}.pdf', group_compare = group_comparesList)
    png = expand(volcanoPlot_dir + '/Volcanoplot_{group_compare}.png', group_compare = group_comparesList)
    htpdf = expand(heatmap_dir + '/HeatmapPlot_{group_compare}.pdf', 
        group_compare = group_comparesList)
    htpng = expand(heatmap_dir + '/HeatmapPlot_{group_compare}.png', 
        group_compare = group_comparesList)  
    output_items.append(de_res)
    output_items.append(deg)
    output_items.append(volcanoplotData)
    output_items.append(pdf)
    output_items.append(png)
    output_items.append(htpdf)
    output_items.append(htpng)             
    return(output_items)

def NMDS_output(config):
    output_items = [] 
    NMDS_dir = config['OutputDir']['GroupCompare']['NMDS']
    NMDSdata = os.path.join(NMDS_dir, 'NMDS_data.xls')
    NMDS_pdf = os.path.join(NMDS_dir, 'NMDS_plot.pdf')
    NMDS_png = os.path.join(NMDS_dir, 'NMDS_plot.png')
    output_items.append(NMDSdata)
    output_items.append(NMDS_pdf)
    output_items.append(NMDS_png)
    return(output_items) 

def PCA_output(config):
    output_items = [] 
    PCA_dir = config['OutputDir']['GroupCompare']['PCA']
    PCAdata = os.path.join(PCA_dir, 'PCA_data.xls')
    PCA_pdf = os.path.join(PCA_dir, 'PCA_plot.pdf')
    PCA_png = os.path.join(PCA_dir, 'PCA_plot.png')
    output_items.append(PCAdata)
    output_items.append(PCA_pdf)
    output_items.append(PCA_png)
    return(output_items) 

def GroupCompare_output(config):
    input_list = []
    input_list.append(NMDS_output(config))
    input_list.append(PCA_output(config))
    return(input_list)

def enrichment_output(config):
    group_compares = config['GroupCompare']
    group_comparesList = get_group_compareList(group_compares)
    ERKEGG_dir = config['OutputDir']['ENRICHMENT']['KEGG']
    ERGO_dir = config['OutputDir']['ENRICHMENT']['GO']
    output_items = [] 
    KEGG_pdf = expand(ERKEGG_dir + '/DEG_KEGG_enrichment_{group_compare}.pdf',
     group_compare = group_comparesList)
    KEGG_png = expand(ERKEGG_dir + '/DEG_KEGG_enrichment_{group_compare}.png',
     group_compare = group_comparesList)
    GO_pdf = expand(ERGO_dir + '/DEG_GeneOntology_BP_enrichment_{group_compare}.pdf',
     group_compare = group_comparesList)
    GO_png = expand(ERGO_dir + '/DEG_GeneOntology_BP_enrichment_{group_compare}.png',
     group_compare = group_comparesList)    
    output_items.append(KEGG_pdf)
    output_items.append(KEGG_png)
    output_items.append(GO_pdf)
    output_items.append(GO_png)
    return(output_items) 

# def GSEA_output(config):
#     group_compares = config['GroupCompare']
#     group_comparesList = get_group_compareList(group_compares)
#     GSEAKEGG_dir = config['OutputDir']['GSEA']['KEGG']
#     output_items = [] 
#     KEGG_data = expand(GSEAKEGG_dir + '/KEGG_GSEA_data_{group_compare}.xls', group_compare = group_comparesList)
#     KEGG_pdf = expand(GSEAKEGG_dir + '/KEGG_GSEA_plot_{group_compare}.pdf', group_compare = group_comparesList)
#     KEGG_png = expand(GSEAKEGG_dir + '/KEGG_GSEA_plot_{group_compare}.png', group_compare = group_comparesList)  
#     output_items.append(KEGG_pdf)
#     output_items.append(KEGG_png)
#     output_items.append(KEGG_data)    
#     return(output_items)

def build_input_running(config):
    input_list = []
    if config['Running']['QC']:
        input_list.append(QC_output(config))
    if config['Running']['ALIGNMENT']:
        input_list.append(alignment_output(config))
    if config['Running']['Count']:
        input_list.append(feature_output(config))
    if config['Running']['GroupCompare']:
        input_list.append(GroupCompare_output(config))
    if config['Running']['DEA']:
        input_list.append(DEA_output(config))
    if config['Running']['ENRICHMENT']:
        input_list.append(enrichment_output(config))
    return(input_list)

def report_output(config):
    report_path = config['ResourcePath']['ReportPath']
    output_item = os.path.join(report_path, 'RSAP_analysis_result_report.md')
    return output_item

def build_input_all(config):
    input_list = []
    input_list.append(build_input_running(config))
    input_list.append(report_output(config))
    return(input_list)
