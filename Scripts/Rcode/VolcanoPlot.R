# Export LabeledVolcanoPlot -----------------------------------------------
# Target_tag, gene (target) is which group
# expLevel_mark, gene (target) is positive (up) or negative (down)
prepareDataVolcanoPlot <- function(de_results, compare_element, 
                                   Target_group = 'NULL', reverseCompare){
  if(reverseCompare == TRUE){
    Target_tag <- rep(compare_element[2], dim(de_results)[1])
    Target_tag[de_results$logFC < 0] <- compare_element[1]
  }else{
    Target_tag <- rep(compare_element[1], dim(de_results)[1])
    Target_tag[de_results$logFC < 0] <- compare_element[2]
  }
  expLevel_mark <- rep('positive_label', dim(de_results)[1])
  expLevel_mark[de_results$logFC < 0] <- 'negative_label'
  Target_group_membership.s <- Target_group
  VolcanoPlot_data <- data.frame(Target.tag = Target_tag,
                                 Target.expLevel = expLevel_mark,
                                 Target.group.membership.s = Target_group_membership.s,
                                 Target.Name = rownames(de_results),
                                 Log2 = de_results$logFC,
                                 Pvalue = de_results$PValue,
                                 Adjusted.pvalue = de_results$adj.PVal,
                                 X.log10.pvalue = log10(de_results$PValue),
                                 X.log10.adjusted.pvalue = log10(de_results$adj.PVal)
  )
  return(VolcanoPlot_data)
}

#' mainVolcanoPlot
#'
#' create volcano plot figure using user input variables and DE results from DSPDA
#' @param de data frame of DE results from DSPDA
#' @return ggFigure ggplot of volcano plot
#' @export
mainVolcanoPlot <- function(de_results, 
                            plot_title = '',
                            pdf_file = pdf_file, 
                            png_file = png_file,
                            VolcanoPlot_data = VolcanoPlot_data,
                            para_group = para_VolcanoPlot){
  # Libraries:
  library(ggplot2)
  library(ggrepel)
  library(testthat)
  library(scales)
  library(stats)
  library(stringr)
  plot_title <- plot_title 
  if(plot_title == ''){
    group_item <- de_results$Target.tag[!duplicated(de_results$Target.tag)]
    plot_title <- paste(group_item[1], group_item[2], sep = '_')
  }
  plot_width <- para_group$plot_width
  plot_height <- para_group$plot_height
  # main function 
  # test for valid input variables
  #testVariableFormats(de=de_results)
  # calculate FDR and add it as a column
  de_results$FDR <- p.adjust(de_results$Pvalue, method="fdr")
  # create volcano plot
  volcanoPlot_results <- plotVolcano(de=de_results)
  ggsave(file = pdf_file, plot = volcanoPlot_results$plot, 
    width=plot_width, height=plot_height, unit="cm")
  ggsave(file = png_file, plot = volcanoPlot_results$plot, 
    width=plot_width, height=plot_height, unit="cm")
  write.table(x=volcanoPlot_results$gene_labels, file=VolcanoPlot_data,
              sep="\t", quote=FALSE, row.names=FALSE)
}

#' plotVolcano
#'
#' create volcano plot figure using user input variables and DE results from DSPDA
#' @param de data frame of DE results from DSPDA
#' @return ggFigure ggplot of volcano plot
#' @export
plotVolcano <- function(de, para_group = para_VolcanoPlot){
  # para
  # User Options
  n_genes <- para_VolcanoPlot$n_genes
  pval_thresh <- para_VolcanoPlot$pval_thresh
  fdr_thresh <- para_VolcanoPlot$fdr_thresh
  fc_thresh <- para_VolcanoPlot$fc_thresh
  plot_title <- para_VolcanoPlot$plot_title 
  gene_list <- para_VolcanoPlot$gene_list 
  show_legend <- para_VolcanoPlot$show_legend
  label_fc <- para_VolcanoPlot$label_fc
  target_groups <- para_VolcanoPlot$target_groups 
  # stable parameter
  # number of target groups 
  color_options <-  c("#3A6CA1", "#FFD861", "#CF4244", "#47BAB4", "#474747", "#EB739E", 
                      "#318026", "#A66293", "#F28E2B", "#8F6954")
  default_color <- "grey65"
  fc_color <- "grey30"
  font_size <- 8
  label_size <- 2
  font_family <- "sans"
  # x axis label
  negative_label <- de$Target.tag[de$Target.expLevel == 'negative_label'][1]
  positive_label <- de$Target.tag[de$Target.expLevel == 'positive_label'][1]
  
  # determine highest x point to make volcano plot equal on both sides
  maxFC <- max(abs(de$Log2))
  maxPval <- min(de$Pvalue)
  # find closest FDR value to fdr_thresh and use that pvalue to add y axis cutoff line
  if(!is.null(fdr_thresh)){
    fdr_pval <- mean(de$Pvalue[which(abs(de$FDR - fdr_thresh) == 
                                       min(abs(de$FDR - fdr_thresh)))])
  }else{
    fdr_pval <- NULL
  }
  # create basic volcano plot with correct formatting
  ggFigure <- ggplot(de, aes(x=Log2, y=Pvalue))+
    geom_point(color=default_color)+
    labs(y="adj.Pvalue",
         x=paste(negative_label, "<-", "log2(FC)", "->", positive_label, sep=" "),
         title=plot_title)+
    theme_bw(base_size = font_size) +
    theme(text = element_text(family = font_family))+
    scale_x_continuous(limits=c(-maxFC, maxFC))
  # this makes for easier testing if not running in DSPDA, will flip yaxis of graph
  # scale_y_continuous(trans=change_axis_revlog_trans(base=10),
  #                    labels=function(x) format(x, trim=TRUE, digits=4,
  #                                              scientific=ifelse(maxPval < 0.0001, TRUE, FALSE), 
  #                                              drop0trailing=TRUE)) 
  if(show_legend == FALSE){
    ggFigure <- ggFigure + theme(legend.position="none")
  }
  # subset de to only include genes either in specified target groups or above pval/fdr threshold
  if(!is.null(target_groups)){
    probe_groups <- strsplit(de$Target.group.membership.s, split=", ")
    de_subset_list <- lapply(X=target_groups, 
                             FUN=function(x){
                               #return de rows with specified target group named in column 
                               return(cbind(de[grep(x=de$Target.group.membership.s, pattern=x),], x))
                             })
    names(de_subset_list) <- target_groups
    gene_coloring <- as.data.frame(do.call(rbind, de_subset_list))
    names(gene_coloring)[names(gene_coloring) == "x"] <- "Target_coloring"
    gene_coloring$Target_coloring <- str_wrap(gene_coloring$Target_coloring, width=45)
    color_label <- "Target Group\nMembership"
  }else{
    color_options <- color_options[1:2]
    if(!is.null(fdr_thresh) & !is.null(pval_thresh)){
      if(fdr_pval < pval_thresh){
        gene_coloring <- de[which(de$Pvalue < pval_thresh),]
        label_thresh_low <- paste("pval <", pval_thresh)
        label_thresh_high <- paste("adj.Pvalue <", fdr_thresh)
        high_thresh <- fdr_pval
      }else{
        gene_coloring <- de[which(de$Pvalue < fdr_pval),]
        label_thresh_low <- paste("adj.Pvalue <", fdr_thresh)
        label_thresh_high <- paste("pval <", pval_thresh)
        high_thresh <- pval_thresh
      }
      # label points as positive or negative FC
      gene_coloring$Target_coloring <- ifelse(test=gene_coloring$Log2 < 0, 
                                              yes=negative_label, 
                                              no=positive_label)
      # label points as above pval or FDR threshold
      gene_coloring$Target_coloring <- ifelse(test=gene_coloring$Pvalue >= high_thresh, 
                                              yes=paste(label_thresh_low, gene_coloring$Target_coloring), 
                                              no=paste(label_thresh_high, gene_coloring$Target_coloring))
      # make color options have muted colors for higher threshold
      color_options <- c(color_options, muted(color_options, l=80))
      names(color_options) <- c(paste(label_thresh_high, negative_label), 
                                paste(label_thresh_high, positive_label),
                                paste(label_thresh_low, negative_label), 
                                paste(label_thresh_low, positive_label))
    }else{
      if(is.null(fdr_thresh)){
        gene_coloring <- de[which(de$Pvalue < pval_thresh),]
        label_thresh <- paste("pval <", pval_thresh)
      }else{
        gene_coloring <- de[which(de$FDR < fdr_thresh),]
        label_thresh <- paste("adj.Pvalue <", fdr_thresh)
      }
      gene_coloring$Target_coloring <- ifelse(test=gene_coloring$Log2 < 0, 
                                              yes=paste(label_thresh, negative_label), 
                                              no=paste(label_thresh, positive_label))
      
      names(color_options) <- c(paste(label_thresh, negative_label), 
                                paste(label_thresh, positive_label))
    }
    color_label <- "Significance:"
    # color by fc_thresh if fc_thresh is not NULL
    if(!is.null(fc_thresh)){
      gene_coloring$Target_coloring[abs(gene_coloring$Log2) < fc_thresh] <- paste("FC <", fc_thresh)
      color_options <- c(color_options, fc_color)
      names(color_options)[length(color_options)] <- paste("log2FC <", fc_thresh)
    }
  }
  color_options <- c(color_options, default_color)
  names(color_options)[length(color_options)] <- "Not Specified"
  # add coloring to ggplot
  ggFigure <- ggFigure + geom_point(data=gene_coloring, aes(x=Log2, y=Pvalue, color=Target_coloring))+
    labs(color=color_label)+
    scale_color_manual(values=color_options)
  # add threshold line values to y axis
  yaxis <- data.frame(brk=as.numeric(pretty_breaks(n=4)(0:max(de$X.log10.pvalue))))
  yaxis$brk <- 10^-(yaxis$brk)
  # keep scientific notation if small enough pvalues when changing to character
  yaxis$label <- format(yaxis$brk, trim=TRUE, digits=4,
                        scientific=ifelse(maxPval < 0.0001, TRUE, FALSE), 
                        drop0trailing=TRUE)
  # add threshold lines if thresholds are not NULL
  if(!is.null(fc_thresh)){
    ggFigure <- ggFigure + geom_vline(xintercept=fc_thresh, linetype="dotted")+
      geom_vline(xintercept=-fc_thresh, linetype="dotted")+
      annotate("text", x=fc_thresh+0.35, y=1,family=font_family, size=font_size/3,
               label=paste0("FC=", round(2^fc_thresh, digits=2)))
  }
  if(!is.null(pval_thresh)){
    ggFigure <- ggFigure + geom_hline(yintercept=pval_thresh, linetype="dotted")
    yaxis <- rbind(yaxis, c(pval_thresh, paste0('pval=',pval_thresh)))
  }
  if(!is.null(fdr_thresh)){
    ggFigure <- ggFigure + geom_hline(yintercept=fdr_pval, linetype="dotted")
    yaxis <- rbind(yaxis, c(fdr_pval, paste0('FDR=',fdr_thresh)))
  }
  # order yaxis in increasing value 
  yaxis$brk <- as.numeric(yaxis$brk)
  yaxis <- yaxis[order(yaxis$brk, decreasing=F),]
  # subset de to only contain genes to label on plot, either by user specified genes or top n_genes by pval
  if(!is.null(gene_list)){
    gene_labels <- subset(de, subset=Target.Name %in% gene_list)
  }else{
    # remove gene labels for genes below lowest pvalue cutoff
    gene_labels <- de[which(de$Pvalue < min(fdr_pval, pval_thresh, na.rm = T)),]
    # only label genes above fc_thresh if set by user else only look at pvalue
    if(!is.null(fc_thresh) & label_fc == FALSE){
      gene_labels <- gene_labels[which(abs(gene_labels$Log2) > fc_thresh),]
    }
    # only keep top # of genes by pvalue
    gene_labels <- gene_labels[head(order(gene_labels$Pvalue, decreasing=FALSE), n=n_genes),]
  }  
  # add gene labels to ggplot
  ggFigure <- ggFigure + geom_text_repel(data=gene_labels, aes(x=Log2, y=Pvalue, label=Target.Name), 
                                         family=font_family, force=5, fontface="bold", min.segment.length=0.1,
                                         size=label_size)
  # add y axis labels 
  ggFigure <- ggFigure + scale_y_continuous(trans=change_axis_revlog_trans(base=10), breaks=as.numeric(yaxis$brk),
                                            labels=yaxis$label)
  colnames(gene_labels) <- c("Target tag", "Target group memership/s", "Target Name", "Log2", 
                             "Pvalue", "Adjusted pvalue", "-log10 pvalue", "-log10 adjusted pvalue",
                             "FDR")
  volcanoPlot <- list(ggFigure, gene_labels)
  names(volcanoPlot) <- c("plot", "gene_labels")
  return(volcanoPlot)
}

#' testAreColors
#'
#' checks if all colors in a vector are valid color names
#' @param colors color names 
#' @return TRUE/FALSE statement on valid color status
#' @export
testAreColors <- function(colors) {
  return(sapply(colors, function(X) {
    tryCatch(is.matrix(col2rgb(X)), 
             error=function(e) FALSE)
  }))
}

#' testIdenticalClass
#'
#' raises error if given object is not the assumed variable class
#' @param object object to determine variable class
#' @param object_name name of object for error message
#' @param class_name expected variable class
#' @return None, errors out if class is not expected
#' @export
testIdenticalClass <- function(object, object_name, class_name){
  expect_identical(object=class(object), expected=class_name, 
                   label=paste(object_name, "variable must be",
                               class_name,"\n You've supplied a", class(object)))
}

#' testVariableFormats
#'
#' check all user input variables for class and other validity markers
#' @param de data frame of DE results from DSPDA
#' @return None, errors out if input is not valid
#' @export
testVariableFormats <- function(de=de_results){
  
  # check user input de file for number of columns and column name matching
  expected_col_names <- c("Target.tag", "Target.group.membership.s", "Target.Name", "Log2", 
                          "Pvalue", "Adjusted.pvalue", "X.log10.pvalue", "X.log10.adjusted.pvalue")
  
  expect_identical(object=ncol(de), expected=length(expected_col_names),
                   label="Number of columns in given volcano plot tab delimited file do not match expected. 
                   Make sure file is tab delimited")
  
  expect_identical(object=colnames(de), expected=expected_col_names, 
                   label="Column names in given volcano plot tab delimited file do not match expected.")
  
  # check that all target_groups have at least one gene if not NULL 
  if(!is.null(target_groups)){
    invisible(lapply(X=target_groups, FUN=function(x){
      genes_in_group <- grep(x=de$Target.group.membership.s, pattern=x)
      if(length(genes_in_group) == 0){
        fail(message=paste(x, "is not a valid probe group. Please check spelling or remove from target_groups before continuing"))
      }
    }))
  }
  
  ############################### USER DEFINED VARIABLE CLASS CHECKS ###############################
  numeric_variables <- c("n_genes", "pval_thresh", "fdr_thresh", "fc_thresh", 
                         "font_size", "label_size", "plot_width", "plot_height")
  
  character_variables <- c("plot_title", "gene_list", "target_groups", 
                           "negative_label", "positive_label")
  
  logical_variables <- c("show_legend", "label_fc")
  
  for(v in 1:length(numeric_variables)){
    if(!is.null(eval(parse(text=numeric_variables[v])))){
      testIdenticalClass(object=eval(parse(text=numeric_variables[v])), 
                         object_name=numeric_variables[v], 
                         class_name="numeric")
    }
    
  }
  
  for(v in 1:length(character_variables)){
    if(!is.null(eval(parse(text=character_variables[v])))){
      testIdenticalClass(object=eval(parse(text=character_variables[v])), 
                         object_name=character_variables[v], 
                         class_name="character")
    }
  }
  
  for(v in 1:length(logical_variables)){
    if(!is.null(eval(parse(text=logical_variables[v])))){
      testIdenticalClass(object=eval(parse(text=logical_variables[v])), 
                         object_name=logical_variables[v], 
                         class_name="logical")
    }
  }
  
  # check that either n_genes or gene_list is not NULL
  if(is.null(gene_list) & is.null(n_genes)){
    fail(message="Either n_genes or gene_list must not be NULL
         both n_genes and gene_list are currently NULL")
  }
  
  # check that either pval_thresh or fdr_thresh is not NULL
  if(is.null(pval_thresh) & is.null(fdr_thresh)){
    fail(message="Either fdr_thresh or pval_thresh must not be NULL
         both fdr_thresh and pval_thresh are currently NULL")
  }
  
  # check that gene_list only contains genes in de results
  if(!is.null(gene_list)){
    expect_true(object=all(gene_list %in% de$Target.Name), 
                label="At least one gene in gene_list does not match genes in volcano plot file results")
  }
  
  ################################# USER DEFINED VARIABLE  CHECKS ################################## 
  
  allowed_fonts <- c("serif", "sans", "mono")
  
  # check that font_family is an allowable font
  if(!font_family %in% allowed_fonts){
    fail(message=paste(font_family, "is not a valid font. Allowed fonts are",
                       paste(allowed_fonts, collapse=", ")))
  }
  
  # check that default_color is an allowable color
  if(!testAreColors(colors=default_color)){
    fail(message=paste(paste(default_color, collapse=", "), "is not a valid color"))
  }
  
  # check that fc_color is an allowable color
  if(!testAreColors(colors=fc_color)){
    fail(message=paste(paste(fc_color, collapse=", "), "is not a valid color"))
  }
  
  # check that color_options are all allowable colors
  if(!all(testAreColors(colors=color_options))){
    error_colors <- color_options[which(!testAreColors(colors=color_options))]
    fail(message=paste(paste(error_colors, collapse=", "), "is/are not valid color(s)"))
  }
  
  # check that enough colors are given for labeled target groups
  expect_gte(object=length(color_options), expected=max(length(target_groups), 2), 
             label=paste("Not enough colors were given./n", max(length(target_groups), 2), 
                         "expected, only", length(color_options), "given"))
  
  expected_output_format <- c("png", "jpg", "tiff", "svg", "pdf", "bmp")
  
  if(!output_format %in% expected_output_format){
    fail(message=paste("Output format not in expected list of formats.\n", output_format, 
                       "given\n expected", paste(expected_output_format, collapse=", ")))
  }
}

#' change_axis_revlog_trans
#'
#' reverse log transform axis; used to return pvalue rather than -log10(pvalue) on yaxis
#' @param base base in which logs are computed
#' @return revlog_trans reverse log transformation
#' @export
change_axis_revlog_trans <- function(base=exp(1)){
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  revlog_trans <- trans_new(name=paste0("revlog-", base), 
                            transform=trans, 
                            inverse=inv, 
                            breaks=log_breaks(base=base),
                            domain=c(1e-100, Inf))
  return(revlog_trans)
}

###########################################################################################
# Input
args <- commandArgs(T)
Rdata <- args[1]
group_compare <- args[2]
VolcanoPlot_data <- args[3]
pdf_file <- args[4]
png_file <- args[5]
compare_name <- group_compare
group_compare <- strsplit(group_compare, '_vs_')[[1]]
reverseCompare <- FALSE

load(Rdata)
# de analysis
# plot parameter
para_VolcanoPlot <- list("pval_thresh" = NULL,
                         "fdr_thresh" = 0.05,
                         "fc_thresh" = 1,
                         "n_genes" = 20,
                         "plot_width" = 20,
                         "plot_height" = 16,
                         "plot_title" = "", 
                         "gene_list" = NULL, 
                         "show_legend" = TRUE,
                         "label_fc" = FALSE,
                         "target_groups" = NULL 
)

de_resList <- de_ResList[[compare_name]]
VolcanoPlot <- prepareDataVolcanoPlot(de_resList[[2]], group_compare, 
                Target_group = 'NULL', reverseCompare = reverseCompare)       
mainVolcanoPlot(de_results = VolcanoPlot, para_group = para_VolcanoPlot,
                pdf_file = pdf_file, png_file = png_file, VolcanoPlot_data = VolcanoPlot_data)

