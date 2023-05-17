# function using edgeR performs TPM normaliztion and DEG analysis
# date 2021-08-26
# author by benchenye
# mail-814643229@qq.com
# history
library(edgeR)
library(dplyr)
# reverseCompare, reverse control and treat compare
edgeRDEA <- function(exp_mat, group_map, group_compare, DEA_resFile, deg_file, reverseCompare = FALSE){
  # deg filter parm
  threshold_logFC <- 1 
  threshold_adj.PVal <- 0.05
  # today
  today <- format(Sys.Date(), format = '%Y%m%d')  
  # get special group samples' expMatrix
  group_map <- group_map[(group_map$Group %in% group_compare), ]
  expr_mat <- exp_mat[ ,(colnames(exp_mat) %in% group_map$Sample)]
  group <- factor(group_map$Group, levels=group_compare)
  # TPM normalization
  dgelist <- DGEList(counts = as.matrix(expr_mat), group = group)
  # keep <- rowSums(dgelist$counts) >= 50
  keep <- rowSums(cpm(dgelist) > 1 ) >= 2
  dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]  
  # TMM normalization
  dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
  #plotMDS(dgelist_norm)
  norm_res <- dgelist_norm@.Data[[1]] %>% data.frame()
  log_norm <- log2(norm_res+1)
  # dispersion
  # design <- model.matrix(~group)
  design <- model.matrix(~0+group)
  colnames(design) <- levels(group)
  rownames(design) <- colnames(expr_mat)
  dgelist_norm <- estimateDisp(dgelist_norm, design, robust=TRUE)
  # DE analysis
  fit <- glmFit(dgelist_norm, design)
  # res <- glmLRT(fit)
  cntr.vs.KD <- makeContrasts(paste0(group_compare,collapse = "-"), levels=design)
  if(reverseCompare == TRUE){
    cntr.vs.KD <- cntr.vs.KD*(-1) 
  }
  res <- glmLRT(fit, contrast=cntr.vs.KD)
  # fit2 <- glmLRT(fit, contrast=c(-1,1))
  res <- topTags(res, n=nrow(expr_mat)) %>% data.frame()
  colnames(res) <- c('logFC', 'logCPM', 'LR', 'PValue', 'adj.PVal')

  res <- res %>% dplyr::select(logFC, adj.PVal, PValue, logCPM, LR) 
  # filter differentially expressed genes
  deg_df <- res[(abs(res$logFC) >= threshold_logFC & res$adj.PVal<threshold_adj.PVal), ] 
  out_deg <- cbind(rownames(deg_df), deg_df)
  colnames(out_deg)[1] <- 'GeneSymbol' 
  # save data
  write.table(res, file = DEA_resFile, quote = F, row.names = T, sep = '\t')
  write.table(out_deg, file = deg_file, quote = F, row.names = F, sep = '\t')
  #return data
  out_res <- list(
    'logNorm' = log_norm,
    'deAnalysisRes' = res,
    'degDf' = deg_df
  )
  return(out_res)
}

# Input
args <- commandArgs(T)
Rdata <- args[1]
group_compare <- args[2]
DEA_resFile <- args[3]
deg_file <- args[4]
# DEA_Rdata <- args[5]
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
de_ResList[[compare_name]] <- edgeRDEA(count_mat, groupMap, group_compare, DEA_resFile, deg_file)
# get deg df 
degDfList[[compare_name]] <- de_ResList[[compare_name]]$degDf

save(count_mat, groupMap, de_ResList, degDfList, file = Rdata)







