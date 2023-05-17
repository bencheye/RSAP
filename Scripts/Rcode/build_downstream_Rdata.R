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
RData_name <- args[3]

# input
count_mat <- read.csv(count_file, header = T, row.names=1)
# GroupMapMat
groupMap <- read_delim(group_map, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
save(count_mat, groupMap, file = RData_name)

