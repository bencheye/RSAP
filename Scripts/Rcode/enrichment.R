library(DOSE)
library(org.Hs.eg.db)
library(clusterProfiler)
library(topGO)
library(Rgraphviz)
library(dplyr)
library(gridExtra)
library(ggplot2)
options(stringsAsFactors = F)

# Function ----------------------------------------------------------------
# BP ----------------------------------------------------------------------
BpPlot <- function(eg, organ, GO_BP_pdf, GO_BP_png){
  # theme -------------------------------------------------------------------
  mytheme<-theme_bw()+theme(
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

  if(organ == 'homo_sapiens'){
      organDB = 'org.Hs.eg.db'
    }else if(organ == 'musculus'){
      organDB = 'org.Mm.eg.db'
      }else{
        print('Please check the organ!')
      }
  ego_bp <- enrichGO(gene          = eg$ENTREZID,
                     #universe      = eg$SYMBOL,
                     OrgDb         = organDB,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.2,
                     readable      = TRUE)
  res <- data.frame(ego_bp@result)
  #res <- data.frame(ego_bp1@result)
  numA <- sapply(res$GeneRatio, function(x){return(strsplit(x, '/')[[1]][1])}) %>% as.numeric()
  numB <- sapply(res$GeneRatio, function(x){return(strsplit(x, '/')[[1]][2])}) %>% as.numeric()
  res$GeneRatio <- round(numA/numB, 2)
  res <- res[order(res$Count[1:10]), ]
  res$Description <- factor(res$Description, levels = res$Description)
  p <- ggplot(data = res[1:10, ], 
              mapping = aes(x = GeneRatio, y = Description, colour = p.adjust, size = Count)) + 
    geom_point() + mytheme + labs(y='') + 
    scale_color_gradient(low="#1f77b4", high="#cde64c")
  ggsave(file = GO_BP_pdf, plot = p, width=20, height=16, unit="cm")
  ggsave(file = GO_BP_png, plot = p, width=20, height=16, unit="cm")
  #return(p)
}

# KEGG --------------------------------------------------------------------
KeggPlot <- function(eg, organ, KEGG_pdf, KEGG_png){
  # theme -------------------------------------------------------------------
  mytheme<-theme_bw()+theme(
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
  if(organ == 'homo_sapiens'){
      organDB = 'hsa'
    }else if(organ == 'musculus'){
      organDB = 'mmu'
      }else{
        print('Please check the organ!')
      }
  kk_bp <- enrichKEGG(eg$ENTREZID, organism=organDB,
                      keyType = "ncbi-geneid",
                      pvalueCutoff=0.05, pAdjustMethod="BH",
                      qvalueCutoff=0.1)
  res <- data.frame(kk_bp@result)
  #res <- data.frame(ego_bp1@result)
  numA <- sapply(res$GeneRatio, function(x){return(strsplit(x, '/')[[1]][1])}) %>% as.numeric()
  numB <- sapply(res$GeneRatio, function(x){return(strsplit(x, '/')[[1]][2])}) %>% as.numeric()
  res$GeneRatio <- round(numA/numB, 2)
  res <- res[order(res$Count[1:10]), ]
  res$Description <- factor(res$Description, levels = res$Description) 
  p <- ggplot(data = res[1:10, ], 
              mapping = aes(x = GeneRatio, y = Description, colour = p.adjust, size = Count)) + 
    geom_point() + mytheme + labs(y='') +
    scale_color_gradient(low="#1f77b4", high="#cde64c")
  ggsave(file = KEGG_pdf, plot = p, width=20, height=16, unit="cm")
  ggsave(file = KEGG_png, plot = p, width=20, height=16, unit="cm")    
  #return(p)
}

plotEnrichment <- function(deGeneCharacter, group_compare, output_dir, organ){
  if(!is.character(deGeneCharacter)){
    stop('deGeneCharacter must be the vector!')
  }
  if(organ == 'homo_sapiens'){
    organDB = org.Hs.eg.db
  }else if(organ == 'musculus'){
    organDB = org.Mm.eg.db
    }else{
      print('Please check the organ!')
    }
  # transform the gene id
  eg <- bitr(deGeneCharacter, fromType="SYMBOL", 
              toType="ENTREZID", OrgDb=organDB)
  today <- format(Sys.Date(), format = '%Y%m%d')
  # plot GO BP
  BpPlot(eg, organ, GO_BP_pdf, GO_BP_png)
  #ggsave(bp_name, width = 8, height = 8)
  # plot kegg
  KeggPlot(eg, organ, KEGG_pdf, KEGG_png)
  #ggsave(kegg_name, width = 8, height = 8)  
}

###########################################################################################
# Input
args <- commandArgs(T)
Rdata = args[1]
group_compare = args[2]
GO_BP_pdf = args[3]
GO_BP_png = args[4]
KEGG_pdf = args[5]
KEGG_png = args[6]
organ = args[7]
compare_name <- group_compare
group_compare <- strsplit(group_compare, '_vs_')[[1]]
load(Rdata)
de_resList <- de_ResList[[compare_name]]
# get deg df 
degDf <- de_resList[['degDf']]
# plot
plotEnrichment(rownames(degDf), group_compare, output_dir, organ)

