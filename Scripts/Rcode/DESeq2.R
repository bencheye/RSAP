
options(stringsAsFactors = F)
options(digits = 4)

library(DESeq2)
library(dplyr)
 
# DEA ------------------------------------------------------------------
DESeq2_DEA <- function(count_mat, group_compare, groupMap, DEA_resFile, deg_file){
  # deg filter parm
  threshold_logFC <- 1 
  threshold_adj.PVal <- 0.05
  condition <- factor(groupMap$Group, levels = group_compare)
  col_data <- data.frame(row.names = colnames(count_mat), condition)
  dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = col_data,
                                design = ~ condition)
  dds_filter <- dds[ rowSums(counts(dds)) > 50 ]
  
  dds_filter.sizefactor <- estimateSizeFactors(dds_filter)  
  size_factor <- sizeFactors(dds_filter.sizefactor)
  count_mat <- dds_filter@assays@data@listData$counts
  for(i in 1:dim(count_mat)[2]){
    count_mat[ ,i] <- count_mat[ ,i]/size_factor[i]
  }
  log_norm <- log2(count_mat+1)
  dds_out <- DESeq(dds_filter)
  res <- results(dds_out)
  res <- res[order(res$log2FoldChange, decreasing = T), ]
  res <- data.frame(res)
  res <- res[!is.na(res$padj), ]
  colnames(res) <- c("baseMean", 'logFC', 'lfcSE', 'stat', 'PValue', 'adj.PVal')
  res <- res %>% dplyr::select(logFC, adj.PVal, PValue, baseMean, lfcSE)
  deg_df <-subset(res,(adj.PVal < 0.05) & (abs(logFC) > threshold_logFC))
  out_deg <- cbind(rownames(deg_df), deg_df)
  colnames(out_deg)[1] <- 'GeneSymbol' 
  # save data
  write.table(res, file = DEA_resFile, quote = F, row.names = T, sep = '\t')
  write.table(out_deg, file = deg_file, quote = F, row.names = F, sep = '\t')
  out_res <- list(
    'logNorm' = log_norm,
    'deAnalysisRes' = res,
    'degDf' = deg_df
  )
  return(out_res)
}

# Input -------------------------------------------------------------------
args <- commandArgs(T)
Rdata <- args[1]
group_compare <- args[2]
DEA_resFile <- args[3]
deg_file <- args[4]

compare_name <- group_compare
group_compare <- strsplit(group_compare, '_vs_')[[1]]
load(Rdata)

# de analysis
if(!exists('de_ResList')){
  de_ResList <- list()
}
if(!exists('degDfList')){
  degDfList <- list()
}

# de analysis return result list
de_ResList[[compare_name]] <- DESeq2_DEA(count_mat, group_compare, groupMap, DEA_resFile, deg_file)
# get deg df 
degDfList[[compare_name]] <- de_ResList[[compare_name]]$degDf

save(count_mat, groupMap, de_ResList, degDfList, file = Rdata)
