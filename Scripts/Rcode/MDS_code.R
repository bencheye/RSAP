
MDSPlot <- function(count_mat, group_map, group_compare, MDSdata_output, pdf_output, png_output){
  library(ggplot2)
  library(edgeR)
  library(dplyr)
  today <- format(Sys.Date(), format = '%Y%m%d')
  # get special group samples' expMatrix
  group_map <- group_map[(group_map$Group %in% group_compare), ]
  count_mat <- count_mat[ ,(colnames(count_mat) %in% group_map$Sample)]
  group <- group_map$Group
  sample <- group_map$Sample
  # log_normal_count
  dgelist <- DGEList(counts = as.matrix(count_mat))
  keep <- rowSums(cpm(dgelist) > 1 ) >= 2
  dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
  # TMM normalization
  dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
  # log-norm+1
  norm_res <- dgelist_norm@.Data[[1]] %>% data.frame()
  log_norm <- log2(norm_res+1)
  # mds
  mds <- plotMDS(dgelist_norm)
  var_explained <- (mds$var.explained*100) %>% round(., digits = 0) %>% as.character
  x_lab <- paste0(mds$axislabel, '1', '(', var_explained[1], '%', ')')
  y_lab <- paste0(mds$axislabel, '2', '(', var_explained[2], '%', ')')
  toplot <- data.frame(Dim1 = mds$x, Dim2 = mds$y, Sample = sample, Group = group)
  write.table(toplot, file=MDSdata_output, quote = F, row.names = F, sep = '\t')
  # theme -------------------------------------------------------------------
  mytheme<-theme_bw()+theme(
    panel.grid=element_blank(),
    text=element_text(colour="black"),
    axis.text=element_text(size=rel(1)),
    axis.title=element_text(size=rel(1.2),lineheight=0.4,hjust=0.5),
    panel.border=element_rect(size=0.3,colour="black"),
    plot.title=element_text(size=rel(1.5),lineheight=0.7,hjust=0.5),
    legend.title=element_text(colour="black", size=rel(1.2)),
    legend.text = element_text(size = rel(1.2)),
    legend.key=element_rect(colour = "white"),
    axis.line=element_blank()
  )
  p = ggplot(toplot, aes(Dim1, Dim2, colour = Group)) + geom_point(size = 3) + mytheme +
  xlab(x_lab) + ylab(y_lab) + 
  geom_text(aes(label = Sample, vjust = 1.1, hjust = -0.5), show_guide = FALSE)
  group_name <- ''
  for (i in 1:length(group_compare)) {
    group_name <- paste(group_name, group_compare[i], sep = '_')
  }
  ggsave(file = pdf_output, plot = p, width=20, height=16, unit="cm")
  ggsave(file = png_output, plot = p, width=20, height=16, unit="cm")
  #return(log_norm)
}

# Input
args <- commandArgs(T)
Rdata = args[1]
group_compare = args[2]
MDSdata_output = args[3]
pdf_output = args[4]
png_output = args[5]
# Rdata
# Rdata <- paste0(projPath, '/Result/ExpressionMatrix/geneExpressionMatrix.Rdata')
# Output <- paste0(projPath, '/Result/MDSPlot')
# dir.create(Output)
load(Rdata)
group_compare <- strsplit(group_compare, ' ')[[1]]

# plot mds
MDSPlot(count_mat, groupMap, group_compare, MDSdata_output, pdf_output, png_output)


