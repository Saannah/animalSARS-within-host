#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
###
library(stringr)
filename = args[1]
print(filename)

#check if file is empty
if(file.size(filename) == 0){file.create(paste0(filename, ".parsed.NoSegSites"))
  stop()}

tab = read.delim(filename, header = FALSE)
colnames(tab) = c("CHROM", "POS", "REF", "ALT", "INFO")


### PARSE INFO INTO SEPARATE COLUMNNs
info_obj = str_split(tab$INFO, ";")
DEPTH = vector()
VAF = vector()
ConsensusTag = vector()

for (i in 1:length(info_obj)){
  DEPTH = c(DEPTH, as.numeric(substr(unlist(info_obj[i])[1], (unlist(gregexpr('=', unlist(info_obj[i])[1]))[1]+ 1), nchar(unlist(info_obj[i])[1]))))
  VAF = c(VAF, as.numeric(substr(unlist(info_obj[i])[2], (unlist(gregexpr('=', unlist(info_obj[i])[2]))[1]+ 1), nchar(unlist(info_obj[i])[2]))))
  ConsensusTag = c(ConsensusTag, substr(unlist(info_obj[i])[3], (unlist(gregexpr('=', unlist(info_obj[i])[3]))[1]+ 1), nchar(unlist(info_obj[i])[3])))
}
tab$DEPTH = DEPTH
tab$VAF = VAF
tab$ConsensusTag = ConsensusTag
tab$ALT_COUNT = as.integer(tab$DEPTH * tab$VAF)

write.table(tab, file = paste0(filename, ".parsed"), row.names = FALSE, sep = "\t")
