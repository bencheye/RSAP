# x, is a list of all the gsea analysis result
# geneSetID, is/are the geneset(s) will be plot, can be one value or a vector,
# e.g. 1:10, is to plot top 10 genesets' enrichment result
# compare_text, is group compare text, must be consistent with Fold change
# plotFile output dir

gseaplot4 <- function(x, geneSetID, compare_text, file_name, 
                      color = "green", base_size = 11, ES_geom = "line",
                      rel_heights = c(1.5, 0.5, 1), subplots = 1:3, 
                      Pwidth = 20, Pheight = 16
) 
{
  ES_geom <- match.arg(ES_geom, c("line", "dot"))
  geneList <- position <- NULL
  if(length(geneSetID) == 1) {
    gsdata <- gsInfo(x, geneSetID)
  }else{
    gsdata <- do.call(rbind, lapply(geneSetID, gsInfo, object = x))
  }
  p <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) + theme_classic(base_size) + 
    theme(panel.grid.major = element_line(colour = "grey92"), 
          panel.grid.minor = element_line(colour = "grey92"), 
          panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) + 
    scale_x_continuous(expand = c(0, 0))
  if(ES_geom == "line") {
    es_layer <- geom_line(aes_(y = ~runningScore, color = ~Description), 
                          size = 1) 
  }else{
    es_layer <- geom_point(aes_(y = ~runningScore, color = ~Description), 
                           size = 1, data = subset(gsdata, position == 1)) 
  }
  #########################################################
  p.res <- p + es_layer + theme(legend.position = 'none')  + scale_color_brewer(palette="Paired")
  p.res <- p.res + ylab("Running Enrichment Score") + 
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
          axis.line.x = element_blank(), plot.margin = margin(t = 0.2, r = 0.2, b = 0, l = 0.2, unit = "cm"))
  p.res <- p.res + theme(panel.background = element_rect(fill='white', colour='black'))
  if(length(geneSetID) == 1) {
    xx <- dim(gsdata)[1]*0.8
    yy <- max(gsdata$runningScore)
    nes <- round(x@result$NES[geneSetID], 2)
    nes <- paste0("NES = ", as.character(nes))
    padj <- signif(x@result$p.adjust[geneSetID], 2)
    padj <- paste0("padj = ", as.character(padj))
    p.res <- p.res + 
      annotate("text", x = xx, y = yy*0.95, label = nes, size = 6) +
      annotate("text", x = xx, y = yy*0.88, label = padj, size = 6)
  }
  i <- 0
  for (term in unique(gsdata$Description)) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == 
                   term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1
  }
  p2 <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) + ylab(NULL) + 
    geom_linerange(aes_(ymin = ~ymin, ymax = ~ymax, color = ~Description))  +  
    scale_color_brewer(palette="Paired") +
    theme_classic(base_size) + 
    theme(legend.position = "none", 
          plot.margin = margin(t = -0.1, b = 0, unit = "cm"), 
          axis.ticks = element_blank(), 
          axis.text = element_blank(), 
          axis.line.x = element_blank()) + 
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
  if(length(geneSetID) == 1) {
    v <- seq(1, sum(gsdata$position), length.out = 9)
    inv <- findInterval(rev(cumsum(gsdata$position)), v)
    if (min(inv) == 0) 
      inv <- inv + 1
    col <- c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))
    ymin <- min(p2$data$ymin)
    yy <- max(p2$data$ymax - p2$data$ymin) * 0.3
    xmin <- which(!duplicated(inv))
    xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
    d <- data.frame(ymin = ymin, ymax = yy, xmin = xmin, 
                    xmax = xmax, col = col[unique(inv)])
    p2 <- p2 + geom_rect(aes_(xmin = ~xmin, xmax = ~xmax, 
                              ymin = ~ymin, ymax = ~ymax, fill = ~I(col)), data = d, 
                         alpha = 0.9, inherit.aes = FALSE)
  }
  df2 <- p$data
  df2$y <- p$data$geneList[df2$x]
  p.pos <- p + geom_segment(data = df2, aes_(x = ~x, xend = ~x, 
                                             y = ~y, yend = 0), color = "grey")
  p.pos <- p.pos + ylab("Ranked List Metric") + xlab(compare_text) + 
    theme(plot.margin = margin(t = -0.1, r = 0.2, b = 0.2, l = 0.2, unit = "cm"))
  p.pos <- p.pos + theme(panel.background = element_rect(fill='white', colour='black'))
  p.res <- p.res + ggtitle(compare_text) +  theme(plot.title = element_text(hjust = 0.5))
  if (length(color) == length(geneSetID)) {
    p.res <- p.res + scale_color_manual(values = color)
    if (length(color) == 1) {
      p.res <- p.res + theme(legend.position = "none")
      p2 <- p2 + scale_color_manual(values = "black")
    }
    else {
      p2 <- p2 + scale_color_manual(values = color)
    }
  }
  plotlist <- list(p.res, p2, p.pos)[subplots]
  n <- length(plotlist)
  plotlist[[n]] <- plotlist[[n]] + theme(axis.line.x = element_line(), 
                                         axis.ticks.x = element_line(), axis.text.x = element_text())
  if(length(subplots) == 1) 
    return(plotlist[[1]] + theme(plot.margin = margin(t = 0.2, 
                                                      r = 0.2, b = 0.2, l = 0.2, unit = "cm")))
  if(length(rel_heights) > length(subplots)) 
    rel_heights <- rel_heights[subplots]
  pp <- plot_grid(plotlist = plotlist, ncol = 1, align = "v", rel_heights = rel_heights)
  if(length(geneSetID)>1){
    mycolor <- RColorBrewer::brewer.pal(length(geneSetID), "Set2")
    len_var <- length(geneSetID)
    if(len_var < 3){
      mycolor <- mycolor[1:len_var]
    }
    legend <- ggplot() +
      annotate("point", x=1,y=1:length(geneSetID)/5,shape=15, color=mycolor, size = 3) +
      annotate("text", x=1.01, y=1:length(geneSetID)/5, label=x@result$ID[geneSetID], hjust=0) +
      xlim(0.99, 1.2) + ylim(0, (length(geneSetID)/5+1)) + theme_void()
    pp <- plot_grid(plotlist = list(pp, legend)[1:2], ncol = 2)
  }
  pp
  ggsave(file = file_name, width = Pwidth, height = Pheight)
  #return(pp)
}



################################################################################

# MSigDb analysis
# H: hallmark gene sets
# C1: positional gene sets
# C2: curated gene sets
# C3: motif gene sets
# C4: computational gene sets
# C5: GO gene sets
# C6: oncogenic signatures
# C7: immunologic signatures

#@ 'analysisType' 'mSigdbr, kegg, go'
#@ 'geneNameType' 'gene_symbol, entrez_gene'
# geneList must be numeric, its value is logFC or FC; its names is genesymbol
# group_compare, group compare information
# reverseCompare, reverse group_compare to match Foldchange group
# gs_cat, choose geneset 
# gs_subcat, choose geneset's subcat
gseaAnalysis <- function(geneList = NULL, 
                         group_compare,
                         GSEA_res_data,
                         pdf_file,
                         png_file,
                         reverseCompare = TRUE,
                         annotationSpecies = 'Homo sapiens',
                         gs_cat = NULL,
                         gs_subcat = NULL,
                         symbolToID = FALSE,
                         geneNameType = 'gene_symbol',
                         parms = GSEA_params
){
  # library dependance ------------------------------------------------------
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(msigdbr)
 # library(pathview)
  library(cowplot)
  library(ggridges)
  # Input
  export_plot <- list()
  group_compare <- group_compare
  if(reverseCompare == TRUE){
    compare_text <- paste0(group_compare[2], ' vs ', group_compare[1])
  }else{
    compare_text <- paste0(group_compare[2], ' vs ', group_compare[1])
  }

  if(!is.numeric(geneList)){
    stop('"geneList" must be numeric!')
  }
  # omit any NA values 
  geneList<-na.omit(geneList)
  # sort the list in decreasing order 
  geneList = sort(geneList, decreasing = TRUE)
  if(symbolToID == TRUE){
    geneList <- geneSymbolToGeneID(geneList = geneList)
    geneNameType <- 'entrezID'
  }
  if(geneNameType != 'gene_symbol' & geneNameType != 'entrez_gene'){
    stop("geneNameType must be 'gene_symbol' OR 'entrez_gene'!")
  }
  if(is.null(gs_cat)){
    m_t2g <- msigdbr(species = annotationSpecies)
  }else{
    m_t2g <- msigdbr(species = annotationSpecies, category = gs_cat) %>% 
      dplyr::select(gs_name, gs_subcat, geneNameType)
  }
  if(!is.null(gs_subcat)){
    uniq_subtype <- table(m_t2g$gs_subcat) 
    uniq_subtype <- names(uniq_subtype) %>% as.character
    filter_subcat <- stringr::str_detect(uniq_subtype, sprintf("%s$", gs_subcat))
    if(sum(filter_subcat) == 0){
      stop(sprintf('Subcategory must be those: %s!', paste(uniq_subtype, collapse = ' ')))
    }
    # filter the subcategory data
    m_t2g <- m_t2g[m_t2g$gs_subcat ==  uniq_subtype[filter_subcat], ]
  }
  m_t2g <- m_t2g[ ,-2]
  # analysis GSEA
  print(paste0(gs_cat, '_', gs_subcat,' : ', 'Starting GSEA analysis!'))
  gsea_res <- GSEA(geneList = geneList, 
                   TERM2GENE = m_t2g, 
                   minGSSize = parms$minGSSize,
                   maxGSSize = parms$maxGSSize,
                   eps = parms$eps,
                   pvalueCutoff = parms$pvalueCutoff,
                   pAdjustMethod = parms$pAdjustMethod,
                   verbose = parms$verbose,
                   seed = parms$seed,
                   by = parms$by
  )
  if(dim(gsea_res@result)[1] == 0){
    print(paste0(gs_cat, '_', gs_subcat,' : ',
                 'There is no term enriched under specific pvalueCutoff'))
    return(0)
  }else{
    print("Starting plot GSEA result picture!")
  }
  # plot
  write.table(gsea_res@result, file = GSEA_res_data, 
              row.names = F, quote = F, sep = '\t')
  # theme
  mytheme <- theme_bw() + theme(
    panel.grid=element_blank(),
    text=element_text(colour="black"),
    axis.text=element_text(size=rel(1)),
    axis.title=element_text(size=rel(1.25),lineheight=0.4,hjust=0.5),
    panel.border=element_rect(size=0.3,colour="black"),
    plot.title=element_text(size=rel(1.5),lineheight=0.7,hjust=0.5),
    legend.title=element_text(colour="white"),
    legend.key=element_rect(colour = "white"),
    axis.line=element_blank()
  )
  if(parms$showCategory > 8){
    showCategory <- 8
  }else if(parms$showCategory > dim(gsea_res@result)[1]){
    showCategory <- dim(gsea_res@result)[1]
  }else{
    showCategory <- parms$showCategory
  }
  gseaplot4(gsea_res, 1:showCategory, 
            file_name = pdf_file, compare_text = compare_text, base_size = 15) 
  gseaplot4(gsea_res, 1:showCategory, 
          file_name = png_file, compare_text = compare_text, base_size = 15) 
}

geneSymbolToGeneID <- function(geneList){
  fromType <- "SYMBOL"
  toType <- "ENTREZID"
  tranformed_res <- bitr(names(geneList), fromType = fromType,
                         toType = toType, OrgDb="org.Hs.eg.db", drop = FALSE)
  # multi-mapping
  dup_name <- tranformed_res[fromType] %>% duplicated(.) %>% tranformed_res[., fromType]
  if(length(dup_name)>0){
    dup_name_collapse <- paste(dup_name, collapse = ' ')
    if(length(dup_name)>1){
      warning(paste0(
        sprintf("There are %d genes with duplicated mapping ID: ", 
                length(dup_name)), dup_name_collapse))
    }
    warning(paste0("There is a gene with duplicated mapping ID: ", dup_name_collapse))
    warning(dup_name_collapse, ' will be removed!')
    tranformed_res <- tranformed_res[!(unlist(tranformed_res[fromType]) %in% dup_name), ]
  }
  # na 
  na_site <- tranformed_res[toType] %>% is.na %>% which
  if(length(na_site) > 0){
    geneList <- geneList[-na_site]
    names(geneList) <- tranformed_res[-na_site, toType]
  }else{
    names(geneList) <- tranformed_res[toType]
  }
  return(geneList)
}

gsInfo <- function(object, geneSetID) {
  geneList <- object@geneList
  if (is.numeric(geneSetID))
    geneSetID <- object@result[geneSetID, "ID"]
  
  geneSet <- object@geneSets[[geneSetID]]
  exponent <- object@params[["exponent"]]
  df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore))/20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList
  df$Description <- object@result[geneSetID, "Description"]
  return(df)
}

gseaScores <- getFromNamespace("gseaScores", "DOSE")

################################################################################
library(dplyr)
# Input
args <- commandArgs(T)
Rdata = args[1]
organ = args[2]
group_compare = args[3]
KEGG_GSEA_res_data = args[4]
KEGG_pdf = args[5]
KEGG_png = args[6]
compare_name <- group_compare
group_compare <- strsplit(group_compare, '_vs_')[[1]]

GSEA_params <- list(
  exponent = 1,
  minGSSize = 30,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea",
  save_as = 'pdf',
  showCategory = 8,
  base_height = 10,
  base_width = 20
)
load(Rdata)
de_resList <- de_ResList[[compare_name]]
# get deg df 
de_results <- de_ResList[[compare_name]][[2]]
# we want the log2 fold change 
geneList <- de_results['logFC'] %>% unlist
# name the vector
names(geneList) <- rownames(de_results)
# omit any NA values 
geneList<-na.omit(geneList)
# sort the list in decreasing order 
geneList = sort(geneList, decreasing = TRUE)

if(organ == 'homo_sapiens'){
  annotationSpecies = 'Homo sapiens'
}else if(organ == 'musculus'){
  annotationSpecies = "Mus musculus"
  }else{
    print('Please check organ name!')
  }
  
symbolToID = FALSE
geneNameType = 'gene_symbol'
# KEGG
gs_cat <- 'C2'
gs_subcat <- 'KEGG'
GSEA_res_data = KEGG_GSEA_res_data
pdf_file = KEGG_pdf
png_file = KEGG_png
gsea_res <- gseaAnalysis(geneList = geneList,
                         group_compare = group_compare,
                         GSEA_res_data,
                         pdf_file,
                         png_file,  
                         reverseCompare = FALSE,
                         annotationSpecies = annotationSpecies,
                         gs_cat = gs_cat,
                         gs_subcat = gs_subcat,
                         symbolToID = symbolToID,
                         geneNameType = geneNameType,
                         parms = GSEA_params
) 
# # BP
# gs_cat <- 'C5'
# gs_subcat <- 'BP'
# GSEA_res_data = GO_GSEA_res_data
# pdf_file = GO_BP_pdf
# png_file = GO_BP_png
# gsea_res <- gseaAnalysis(geneList = geneList,
#                          group_compare = group_compare,
#                          GSEA_res_data,
#                          pdf_file,
#                          png_file,
#                          reverseCompare = FALSE,
#                          annotationSpecies = annotationSpecies,
#                          gs_cat = gs_cat,
#                          gs_subcat = gs_subcat,
#                          symbolToID = symbolToID,
#                          geneNameType = geneNameType,
#                          parms = GSEA_params
# ) 