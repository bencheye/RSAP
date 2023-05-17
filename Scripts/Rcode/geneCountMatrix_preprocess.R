# function using edgeR performs TPM normaliztion and DEG analysis
# date 2021-08-26
# author by benchenye
# mail-814643229@qq.com
# history

library(edgeR)
library(dplyr)
library(readr)
options(stringsAsFactors = F)
args=commandArgs(T)

count_file <- args[1] 
group_map <- args[2]
reserved_item <- args[3]
out_name <- args[4]

# input
count_mat <- read.csv(count_file, header = T, sep = ' ')
# filter gene_type
count_mat <- count_mat[(count_mat$gene_type %in% reserved_item), ]
##################
count_mat <- count_mat[ ,3:(dim(count_mat)[2])]
if((is.na(count_mat[[dim(count_mat)[2]]]) %>% sum) == dim(count_mat)[1]){
	count_mat <- count_mat[ ,-dim(count_mat)[2]]
	}
# sum counts according to the gene_name
count_mat <- aggregate(x = count_mat[,2:ncol(count_mat)],
                          by = list(count_mat$gene_name),
                          FUN = sum)
                          dim(count_mat)
rownames(count_mat) <- count_mat[ ,1]
count_mat <- count_mat[ ,-1]
# GroupMapMat
groupMap <- read_delim(group_map, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
write.csv(count_mat, file = out_name, quote = F)

