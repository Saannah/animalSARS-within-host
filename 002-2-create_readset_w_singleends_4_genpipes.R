## script to create readset files using fastq files, including both single end and paired end
## this script inputs the "realpath" for all fastqs poold together, both single end and paired end
library(stringr)
setwd("/Users/sana/Documents/AIM3/002-readset_files/")
### load Hector's example for help
ex = read.delim("example_readset_from_Hector.txt")

### load the tsv file containing 
paths = read.delim("all_sample_paths.tsv", header = FALSE)
paths = paths[, 1]

#remove paths that are from the test folder
paths = paths[-which(str_detect(paths, "test"))]

##parse paths to extract accessions
obj = str_split(paths, "/")
accessions = vector()
for(i in 1:length(obj)){
  accessions = c(accessions, unlist(obj[i])[10])
}

uniq_accessions = unique(accessions)
run_type = vector()
for(i in uniq_accessions){
  run_type = c(run_type, length(which(accessions == i)))
}
#the order of the array "accession" is the same order as runtype
readset = data.frame(matrix(ncol = length(colnames(ex)), nrow = length(uniq_accessions)))
colnames(readset) <- colnames(ex)
readset$Sample = uniq_accessions
readset$Readset = uniq_accessions
readset$Library = uniq_accessions
readset$Run = unique(ex$Run)
readset$Lane = unique(ex$Lane)
readset$Adapter1 = unique(ex$Adapter1)
readset$Adapter2 = unique(ex$Adapter2)
readset$QualityOffset = unique(ex$QualityOffset)
readset$BED = unique(ex$BED)
readset$BAM = unique(ex$BAM)

## run type
run_type = as.character(run_type)
run_type = str_replace_all(run_type, "1", "SINGLE_END")
run_type = str_replace_all(run_type, "2", "PAIRED_END")
readset$RunType = run_type

## set paths for fastq1 and fastq2
for(a in readset$Sample){
  if(length(which(str_detect(paths, a))) == 1){
    readset$FASTQ1[which(readset$Sample == a)] = paths[which(str_detect(paths, a))]
  }
  if(length(which(str_detect(paths, a))) == 2){
    sub_paths = paths[which(str_detect(paths, a))]
    readset$FASTQ1[which(readset$Sample == a)] = sub_paths[which(str_detect(sub_paths, "1.fastq"))]
    readset$FASTQ2[which(readset$Sample == a)] = sub_paths[which(str_detect(sub_paths, "2.fastq"))]
  }
}

#write the readset file
write.table(readset, file = "readset_for_paired_and_single.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
