# Categorize different elements by the distance between them
# Pre-screened differentially expressed genes were found, with the horizontal axis as samples and the vertical axis as genes
# Clinical information of samples (HCC liver cancer, NOR, normal group)
# row cluster tree group color number, no group set to 0, default 0
# group split by column
# lengend name--value mean
####################################################################
# dataset
# groupInfo, groupMap, group-sample mapping info
# group_colName, groupName, e.g.'Group'
# group_compare, group compare
# group_annotation
# mark_gene, need to mark gene symbol
heatmapPlot <- function(dataset, 
                        groupInfo = NULL, 
                        group_colName = NULL,
                        group_compare = NULL,
                        group_annotation = NULL, 
                        mark_gene = NULL,
                        pdf_file = pdf_file,
                        parms
){
  # library dependance
  library(ComplexHeatmap)
  library(dendextend)
  library(dplyr)
  # plot parameter
  row_dend_colorDivNum = params$row_dend_colorDivNum
  col_dend_colorDivNum = params$col_dend_colorDivNum
  lengend_name <- params$lengend_name
  show_column_names = params$show_column_names
  show_row_names = params$show_row_names
  Plotwidth <- params$Plotwidth 
  Plotheight <- params$Plotheight
  # filter data
  if(!is.null(group_compare) & length(group_compare)>0){
    chosed_samples <- unlist(groupInfo[group_colName]) %in% group_compare
    if(sum(chosed_samples) != dim(dataset)[2]){
      dataset <- dataset[ ,chosed_samples]
      groupInfo <- groupInfo[chosed_samples, ]
    }
  }
  # filter show genes 
  show_genes_num <- params$show_gene_num
  if((!is.null(mark_gene)) & (length(mark_gene)>show_genes_num)){
    mark_gene <- mark_gene[1:show_genes_num]
  }
  colors_list <- c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#33A02C",
    "#A65628", "#F781BF", "#E5C494", "#B3B3B3", "#B3DE69", "#FCCDE5",
    "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", 
    "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",   
    "#A6CEE3", "#B2DF8A", "#1F78B4", "#E31A1C", "#6A3D9A", "#BC80BD",
    "#CCEBC5", "#FFED6F", "#B15928", "#999999", "#D9D9D9", "#FF7F00", 
    "#FB9A99", "#FDBF6F", "#CAB2D6", "#FFFF99", "#FFFF33"
  )
  color_usedNum <- 0
  # group annotation color for samples--columns
  if(length(group_annotation)>0){
    columns_annotation_color <- list()
    for(len in 1:length(group_annotation)){
      uniq_var <- groupInfo[group_annotation[len]] %>% unique %>% unlist
      var_color <- colors_list[(color_usedNum+1) : (length(uniq_var)+color_usedNum)]
      names(var_color) <- uniq_var
      color_usedNum <- color_usedNum + length(uniq_var)
      columns_annotation_color[[group_annotation[len]]] = var_color
      # assign the group elements into groupName
      #assign(group_annotation[len], var_list)
    }
    # group annotation dataFrame
    annotation_df <- groupInfo[colnames(groupInfo) %in% group_annotation] %>% data.frame
    top_annotation = HeatmapAnnotation(df = annotation_df,
                                       col = columns_annotation_color,
                                       annotation_name_side = "right"
    )
  }else{
    top_annotation = NULL
  }
  # 突出重要元素(基因)
  # gene_pos 需要与 mark_gene 对应
  if(length(mark_gene)>0){
    gene_pos <- c()
    for(i in 1:length(mark_gene)){
      gene_pos[i] <- which(rownames(dataset) == mark_gene[i])
    }
    row_anno <-  rowAnnotation(mark_gene = anno_mark(at = gene_pos, labels = mark_gene))
  }else{
    row_anno <- NULL
  }
  # 根据特定分组来分割热�?  if(length(group_colName)>0){
    group_split <- groupInfo[group_colName] %>% unlist
    if(col_dend_colorDivNum>=0){
      group_split <- length(group_compare)
  }else{
    group_split <- NULL
  }
  # 修改聚类树不同分组的颜色
  if(row_dend_colorDivNum > 0){
    row_dend <- as.dendrogram(hclust(dist(dataset)))
    # `color_branches()` returns a dendrogram object
    row_dend = dendextend::color_branches(row_dend, k = row_dend_colorDivNum)
  }else{
    row_dend <- FALSE
  }
  if(col_dend_colorDivNum > 0){
    col_dend <- as.dendrogram(hclust(dist(t(dataset))))
    col_dend = dendextend::color_branches(col_dend, k = col_dend_colorDivNum)
  }else if(col_dend_colorDivNum < 0){
    col_dend <- FALSE
  }else{
    if(length(group_compare)>0){
      col_dend <- as.dendrogram(hclust(dist(t(dataset))))
      col_dend = dendextend::color_branches(col_dend, k = length(group_compare)) 
    }else{
      col_dend <- FALSE
    }
  }
  # plot heatmap global parameter
  ht_opt(
    #heatmap_column_names_gp = gpar(fontface = "italic"), 
    #heatmap_column_title_gp = gpar(fontsize = 10),
    legend_border = "black",
    heatmap_border = TRUE,
    annotation_border = TRUE
  )
  if(is.null(lengend_name)){
    lengend_name <- 'Value'
  }
  # plot
  ht_list <- Heatmap(dataset, 
                     name = lengend_name, 
                     column_split = group_split, 
                     top_annotation = top_annotation,
                     show_column_names = show_column_names,
                     show_row_names = show_row_names,
                     column_title = NULL,
                     cluster_rows = row_dend,
                     cluster_columns = col_dend,
                     right_annotation = row_anno,
                     heatmap_legend_param = list(title_position = "leftcenter-rot")
  )
  print(pdf_file)
  # Output
  pdf(file = pdf_file, width = Plotwidth, height = Plotheight)
  draw(ht_list, merge_legend = TRUE)
  dev.off()

}

###############################################################################################
# Input
args <- commandArgs(T)
Rdata = args[1]
group_compare = args[2]
pdf_file = args[3]
compare_name <- group_compare
group_compare <- strsplit(group_compare, '_vs_')[[1]]
load(Rdata)
# plot parameter
params <- list(
        "lengend_name" = 'GeneExpression',
        "row_dend_colorDivNum" = 10,
        "col_dend_colorDivNum" = 0,
        "show_row_names" = FALSE,
        "show_column_names" = TRUE,
        "Plotwidth" = 20,
        "Plotheight" = 16,
        "show_gene_num" = 20
)


de_resList <- de_ResList[[compare_name]]
# get deg df 
degDf <- de_resList[['degDf']]
log_norm <- de_resList[['logNorm']]
# plot
heatmapPlot(dataset = log_norm, 
            groupInfo = groupMap,
            group_colName = 'Group', # 'Group'
            group_compare = group_compare,
            group_annotation = 'Group', 
            mark_gene = rownames(degDf),
            pdf_file = pdf_file,
            parms = params
)


