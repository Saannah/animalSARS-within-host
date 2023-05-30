#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

##script for generating readset file
#FILE NEEDED: list of accessions that have paired-end reads
species = args[1]

## step 1: load Hector's example readset file:
readset_example = read.delim("/lustre06/project/6006375/sanna/aim3/genpipes/1-readset_file/0-example_from_hector/example_readset_from_Hector.txt")
doubles = read.delim(paste0("/home/sanna/projects/def-shapiro/sanna/aim3/genpipes/1-readset_file/1-paired_end_accessions_list/",species,"_doubles.txt"), header = FALSE)
doubles = doubles[,1]
library(stringr)
doubles = str_split(doubles, "/")
doubles_accessions = vector()
for(i in 1:length(doubles)){
  doubles_accessions = c(doubles_accessions, unlist(doubles[i])[1])
}

readset_final = data.frame(matrix(ncol = length(colnames(readset_example)), nrow = length(doubles_accessions)))
colnames(readset_final) <- colnames(readset_example)
readset_final$Sample = doubles_accessions
readset_final$Readset = doubles_accessions
readset_final$Library = doubles_accessions
readset_final$RunType = unique(readset_example$RunType)
readset_final$Run = unique(readset_example$Run)
readset_final$Lane = unique(readset_example$Lane)
readset_final$Adapter1 = unique(readset_example$Adapter1)
readset_final$Adapter2 = unique(readset_example$Adapter2)
readset_final$QualityOffset = unique(readset_example$QualityOffset)
readset_final$BED = unique(readset_example$BED)
readset_final$BAM = unique(readset_example$BAM)
#create path for fastq1
path_fastq_1 = paste0("/home/sanna/projects/def-shapiro/sanna/aim3/sra/", species, "/", doubles_accessions, "/", 
                      doubles_accessions, "_1.fastq")
path_fastq_2 = paste0("/home/sanna/projects/def-shapiro/sanna/aim3/sra/", species, "/", doubles_accessions, "/", 
                      doubles_accessions, "_2.fastq")

readset_final$FASTQ1 = path_fastq_1
readset_final$FASTQ2 = path_fastq_2

#write table
write.table(readset_final, file = paste0(species, "_readset.txt"), row.names = FALSE, sep = "\t", quote = FALSE)
