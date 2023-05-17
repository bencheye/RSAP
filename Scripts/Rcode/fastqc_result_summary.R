
# function
getOutput <- function(sample, file_name){
  out_list <- list('Sample' = sample)
  content <- readLines(file_name)
  # Sequence length
  site <- grepl('^Sequence length', content) %>% which
  item <- content[site]
  out_list['Length(bp)'] <- str_split(item, "\t", simplify = TRUE)[1,2] %>% as.numeric()
  # %GC
  site <- grepl('^%GC', content) %>% which
  item <- content[site]
  out_list['GC(%)'] <- str_split(item, "\t", simplify = TRUE)[1,2] %>% as.numeric()
  # Per sequence quality scores
  start_site <- grepl('^#Quality', content) %>% which
  stop_site <- grepl('^>>END_MODULE', content) %>% which
  stop_site <- stop_site[stop_site > start_site][1]
  # filter
  quality_score <- content[(start_site+1):(stop_site-1)]
  quality_score <- str_split(quality_score, "\t", simplify = TRUE) %>% data.frame()
  colnames(quality_score) <- c('Score', 'Reads')
  quality_score$Score <- as.numeric(quality_score$Score)
  quality_score$Reads <- as.numeric(quality_score$Reads)
  # reads number
  out_list['Reads'] <- sum(quality_score$Reads)
  # Q20
  out_list['Q20(%)']<- sum(quality_score$Reads[quality_score$Score>=20])  #/as.numeric(out_list[['reads']])*100
  # Q30
  out_list['Q30(%)'] <- sum(quality_score$Reads[quality_score$Score>=30])  # /as.numeric(out_list[['reads']])*100
  return(out_list)
}
# 
qcSummary <- function(R1_res, R2_res){
  res_list <- list()
  res_list['Sample'] <- R1_res[['Sample']]
  res_list['Length(bp)'] <- (R1_res[['Length(bp)']]+R2_res[['Length(bp)']])/2
  res_list['Reads'] <- R1_res[['Reads']]+R2_res[['Reads']]
  #res_list['Bases(bp)'] <- R1_res[['Bases(bp)']]+R2_res[['Bases(bp)']]
  res_list['GC(%)'] <- (R1_res[['GC(%)']]+R2_res[['GC(%)']])/2
  Q20 <- (R1_res[['Q20(%)']]+R2_res[['Q20(%)']])/ res_list[['Reads']]*100 
  res_list['Q20(%)'] <- round(Q20, 2)
  Q30 <- (R1_res[['Q30(%)']]+R2_res[['Q30(%)']])/res_list[['Reads']]*100 %>% round(., 2)
  res_list['Q30(%)'] <- round(Q30, 2)
  return(res_list)
}

library(dplyr)
library(stringr)
# Analysis
args <- commandArgs(T)
sample_list_name <- args[1]
fastqc_path <- args[2]
out_name <- args[3]
end <- args[4]
paired_sep <- args[5]

if(paired_sep != 'None'){
  paired_sep <- strsplit(paired_sep, ' ')[[1]]
}

sample_list <- read.table(sample_list_name, sep = '\t', header = T)[ ,1]
out_summary <- data.frame()

if(end == 'pair'){
  for(i in 1:length(sample_list)){
    sample <- sample_list[i]
    R1_name <- paste0(fastqc_path, '/', sample, '_', paired_sep[1], '_fastqc', '/fastqc_data.txt')
    R2_name <- paste0(fastqc_path, '/', sample, '_', paired_sep[2], '_fastqc', '/fastqc_data.txt')
    R1_res <- getOutput(sample, R1_name)
    R2_res <- getOutput(sample, R2_name)
    sample_summary <- qcSummary(R1_res, R2_res) %>% unlist()
    out_summary <- rbind(out_summary, sample_summary)
  }
}else{
  for(i in 1:length(sample_list)){
    sample <- sample_list[i]
    R1_name <- paste0(fastqc_path, '/', sample, '_fastqc', '/fastqc_data.txt')
    sample_summary <- getOutput(sample, R1_name) %>% unlist()
    out_summary <- rbind(out_summary, sample_summary)
}
}

colnames(out_summary) <- c('Sample', 'Length(bp)', 'Reads',	 'GC(%)', 'Q20(%)', 'Q30(%)')
write.table(out_summary, file = out_name, sep = '\t', quote = F, row.names = F)



