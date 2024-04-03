#### script to parse vcf files into a table format, and match with metadata
library(stringr)

setwd("/Users/sana/Documents/AIM3/021-full_pipeline/0-starting_files")
filenames = read.delim("filenames", header = FALSE)
filenames = filenames[, 1]
length(filenames)

#create a vector with sample names in the same order as file names
sample_names = vector()
filenames_obj = str_split(filenames, "\\.")
for(i in 1:length(filenames_obj)){
  sample_names = c(sample_names, unlist(filenames_obj[i])[1])
}

#load metadata
metadata = read.delim("metadata_all_w_human_cleaned.tsv")
filenames = filenames[which(sample_names %in% metadata$Run)]
metadata = metadata[which(metadata$Run %in% sample_names), ]
sample_names = sample_names[which(sample_names %in% metadata$Run)]

#parse the matched file names and metadata into tables and save into a separate directory
for(filename in filenames){
  
  #check if file is empty
  if(file.size(filename) == 0){file.create(paste0("/users/sana/Documents/AIM3/021-full_pipeline/",filename, ".unfiltered.parsed.NoSegSitesCalled"))
    next}
  
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
  
  write.table(tab, file = paste0("/users/sana/Documents/AIM3/021-full_pipeline/", filename, ".unfiltered.parsed"), row.names = FALSE, sep = "\t")
  
}


## write matched metadata
write.table(metadata, "matched_metadata.tsv", sep = "\t", row.names = FALSE)


### load bedfile names:
setwd("/users/sana/Documents/AIM3/021-full_pipeline/0-starting_files/")
bednames = read.delim("/users/sana/Documents/AIM3/021-full_pipeline/0-starting_files/bednames", header = FALSE)
bednames = bednames[, 1]
bedobj = str_split(bednames, "\\.")
bedsamples = vector()
for(i in 1:length(bedobj)){
  bedsamples = c(bedsamples, unlist(bedobj[i])[1])
}

bednames = bednames[which(bedsamples %in% metadata$Run)]
for(i in bednames){
  system(paste0("cp ", i, " ../00-ALL_FILES/"))
  
}



###match freyja files
freyjanames = read.delim("/users/sana/Documents/AIM3/021-full_pipeline/0-freyja_files/freyjanames", header = FALSE)
freyjanames = freyjanames[, 1]

samples = vector()
for(i in 1:length(freyjanames)){
  samples = c(samples,
              unlist(str_split(freyjanames[i], "\\."))[1])
}

### load matched metadata
metadata = read.delim("/users/sana/Documents/AIM3/021-full_pipeline/00-ALL_FILES/matched_metadata.tsv")

matched_freyja_filenames = freyjanames[which(samples %in% metadata$Run)]

setwd("/users/sana/Documents/AIM3/021-full_pipeline/0-freyja_files/")
for(file in matched_freyja_filenames){
  system(command = paste0("cp ", file, " ../00-ALL_FILES/"))
}


