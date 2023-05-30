### load the metadata file

metadata = read.delim("/users/sana/Documents/AIM3/016-glmm_w_human_ref/all_metadata_plus_human.tsv")

dim(metadata)
length(unique(metadata$Run))
unique(metadata$Species)
unique(metadata$lab_infected)
unique(metadata$assay_type)
#change WGA to WGS
metadata$assay_type[which(metadata$assay_type == "WGA")] = "WGS"
#the targeted captures are probably amplicon
metadata$assay_type[which(metadata$assay_type == "Targeted-Capture")] = "AMPLICON"

unique(metadata$LibraryLayout)
unique(metadata$Instrument)
unique(metadata$Platform)
unique(metadata$Isolation_source)

metadata$tissue = "nasal/pharyngeal-swab"
metadata$tissue[which(metadata$Isolation_source == "retropharyngeal_lymph_node")] = "lymph-tissue"
unique(metadata$tissue)

write.table(metadata, file = "metadata_all_w_human_cleaned.tsv", sep = "\t", row.names = FALSE)
