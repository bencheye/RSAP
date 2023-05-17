### Differential expression analysis (DEA)

Gene differential expression analysis is based on the quantitative expression value obtained by transcriptome sequence to identify the genes with differential expression changes. For RNASeq data, there are two common methods, `edgeR` (v3.34.1) and `DESeq2` (v1.34.0) , are provided. {} is used here with `adjusted pvalue < 0.05` and `|log2foldchange| > 1` to identify differentially expressed genes. 

The results of differentially expressed analysis of the first ten rows is shown in the table below:

{}

Then, `Volcano plot` and `heatmap` are used to visualize the results of differential expression analysis. 

For `heatmap`, the column corresponds to the samples and the row corresponds to the genes. The color block above is used to distinguish the group information, and the names of the top 20 differentially expressed genes are marked on the right. The color depth of each block in the heatmap is corresponds to its expression value. The heatmap is shown in the figure below.

{}

For `volcano plot`, the abscissa is the value of `log2foldchange`, and the ordinate is the `adjusted p value`. Yellow indicates up-regulated genes and blue shows down regulated genes. The nodes of the top 20 differentially expressed genes in the picture are marked. The volcano plot is shown in the figure below.

{}

