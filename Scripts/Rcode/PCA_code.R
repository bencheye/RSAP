# dimReduction ------------------------------------------------------------
# plot parameters
# Options: tSNE, UMAP, PCA
# options: pdf, jpeg, tiff, png, bmp, or svg
# color_levels = NULL to plot all group

# plot PCA

library(stats)
library(ggplot2)
library(gridExtra)
library(edgeR)  
library(dplyr)
args <- commandArgs(T)

# Input
Rdata = args[1]
group_compare = args[2]
pcadata_output = args[3]
pdf_output = args[4]
png_output = args[5]

load(Rdata)
group_compare <- strsplit(group_compare, ' ')[[1]]

##########original table###################
  groupMap <- groupMap[(groupMap$Group %in% group_compare), ]
  count_mat <- count_mat[ ,(colnames(count_mat) %in% groupMap$Sample)]
  group <- groupMap$Group
  sample <- groupMap$Sample
  # log_normal_count
  dgelist <- DGEList(counts = as.matrix(count_mat))
  keep <- rowSums(cpm(dgelist) > 1 ) >= 2
  dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
  # TMM normalization
  dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
  # log-norm+1
  norm_res <- dgelist_norm@.Data[[1]] %>% data.frame()
  log_norm <- log2(norm_res+1)

#############PCA######################
clean_table <-t(log_norm)
pc_table <-prcomp(clean_table,scale. = TRUE)
gg_table <- as.data.frame(pc_table$x)[1:3]

#pcadata_output <- paste(Output,"PCA_data.xls",sep="/")
write.table(gg_table,file = pcadata_output,sep = "\t",row.names = TRUE,col.names = TRUE)
################map################
#groupMap <-read.csv("b",header = TRUE,sep ="\t",check.names = FALSE)
gg_table['sample'] <- sample
gg_table['group'] <- group
################theme################

mytheme<-theme_bw()+theme(
  panel.grid=element_blank(),
  text=element_text(colour="black"),
  axis.text=element_text(size=rel(1)),
  axis.title=element_text(size=rel(1.2),lineheight=0.4,hjust=0.5),
  panel.border=element_rect(size=0.3,colour="black"),
  plot.title=element_text(size=rel(1.2),lineheight=0.7,hjust=0.5),
  legend.title=element_text(colour="white"),
  legend.text = element_text(size = rel(1.2)),
  legend.key=element_rect(colour = "white"),
  axis.line=element_blank()
)

###################judge#######################
Group_tmp <- unique(groupMap[,2])
Group_num <- length(Group_tmp)
if (Group_num <= 4){legend_position <- c(0.98,0.01);widths_value <- c(4,2);heights_value <- c(4,2)} else if (Group_num < 7 ){legend_position <- "none";widths_value <- c(5,2);heights_value <- c(5,2)} else {legend_position <- "none";widths_value <- c(5,3);heights_value <- c(5,3)}

###################PC1-PC2###########################
pc1_ratio <- summary(pc_table)$importance[2,1]*100
pc2_ratio <- summary(pc_table)$importance[2,2]*100
#pc3_ratio <- summary(pc_table)$importance[3,3]*100 -pc1_ratio-pc2_ratio

pc1_max <-max(gg_table$PC1)
pc1_min <-min(gg_table$PC1)
pc2_max <-max(gg_table$PC2)
pc2_min <-min(gg_table$PC2)

p <- ggplot(data = gg_table,aes(x = PC1,y = PC2,colour = group))
p <- p+geom_point(size=3)+mytheme+labs(x=paste("PC1-Percent variant explained ",round(pc1_ratio,2),"%",sep = ""),
                                 y =paste("PC2-Percent variant explained ",round(pc2_ratio,2),"%",sep = ""))+
  stat_ellipse(aes(color = group),lwd =0.75)+ggtitle("PCA analysis")+
# xlim(pc1_min-(pc1_max-pc1_min)*0.8,pc1_max+(pc1_max-pc1_min)*0.3)+ylim(pc2_min-(pc2_max-pc2_min)*0.8,pc2_max+(pc2_max-pc2_min)*0.5)+
  geom_text(aes(label = sample, vjust = 1.1, hjust = -0.5), show_guide = FALSE)
p
ggsave(file = pdf_output, plot = p, width=20, height=16, unit="cm")
ggsave(file = png_output, plot = p, width=20, height=16, unit="cm")