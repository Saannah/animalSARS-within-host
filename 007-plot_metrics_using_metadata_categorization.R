### script to match calculated stats to the metadata downloaded from NCBI
library(stringr)
#load the metrics dataframe
metrics = read.delim("/Users/sana/Documents/AIM3/variant_plots-old/dataframe.tsv")
n = read.delim("/Users/sana/Dropbox/AIM3/variant_plots/NoSegSites_names_list.txt")
n = n[, 1]
nosegsites_names_list = vector()
for(i in 1:length(n)){
  nosegsites_names_list = c(nosegsites_names_list, unlist(str_split(n, "\\.")[i])[1])
}

#load metadata for all 4 species
setwd("/Users/sana/Dropbox/AIM3/004-metadata_categorization/")
cat = read.delim("cat_final.tsv")
dog = read.delim("dog_final.tsv")
deer = read.delim("deer_final.tsv")
mink = read.delim("mink_final.tsv")
all_metadata = rbind(cat[, 1:3], dog[, 1:3], mink[, 1:3], deer[, 1:3])



#add a column to metrics dataframe and specify species
metrics$species = "TBD"
metrics$species[which(metrics$SRA_ID %in% cat$Run)] = "cat"
metrics$species[which(metrics$SRA_ID %in% dog$Run)] = "dog"
metrics$species[which(metrics$SRA_ID %in% mink$Run)] = "mink"
metrics$species[which(metrics$SRA_ID %in% deer$Run)] = "deer"

#add a columnn to metrics dataframe and specify lab infected or not
metrics$lab_infected = "TBD"
#lab infected cats:


species = c("cat", "dog", "mink", "deer")
for(s in species){
  metadata = read.delim(paste0("/Users/sana/Dropbox/AIM3/004-metadata_categorization/", s, "_final.tsv"))
  metrics$lab_infected[match(metadata[, 1], metrics$SRA_ID)] = metadata[, 2] 
}
metrics = metrics[which(metrics$SRA_ID %in% all_metadata$Run), ]
metrics$lab_infected[which(metrics$lab_infected == "probably no")] = "no"
metrics$lab_infected[which(metrics$lab_infected == "no")] = "Natural Infection"
metrics$lab_infected[which(metrics$lab_infected == "yes")] = "Lab Infected"
metrics = metrics[-which(metrics$lab_infected == "Missing data"), ]

library(ggplot)
library(tidyverse)
#plot number of segregating sites
p = ggplot(metrics, aes(x=species, y=n_segregating_sites, fill = lab_infected)) + 
  geom_jitter(aes(color = lab_infected), width = 0.2, size = 1) + 
  scale_color_manual(values = c("#DDA7AB", "#4A5899")) +
  theme_bw()
ggsave(p, filename = "n_segsites.pdf", device = "pdf")

#plot average depth
p = ggplot(metrics, aes(x=species, y=average_depth, fill = lab_infected)) + 
  geom_jitter(aes(color = lab_infected), width = 0.2, size = 1) + 
  scale_color_manual(values = c("#DDA7AB", "#4A5899")) +
  theme_bw()
ggsave(p, filename = "average_depth.pdf", device = "pdf")

numeric_dataframe = metrics
numeric_dataframe$pi = as.numeric(numeric_dataframe$pi)
numeric_dataframe$tajimas_D = as.numeric(numeric_dataframe$tajimas_D)

#plot pi
p = ggplot(numeric_dataframe, aes(x=species, y=pi, fill = lab_infected)) + 
  geom_jitter(aes(color = lab_infected), width = 0.2, size = 1) + 
  scale_color_manual(values = c("#DDA7AB", "#4A5899")) +
  theme_bw()
ggsave(p, filename = "pi.pdf", device = "pdf")

#plot tajimas D
p = ggplot(numeric_dataframe, aes(x=species, y=tajimas_D, fill = lab_infected)) + 
  geom_jitter(aes(color = lab_infected), width = 0.2, size = 1) + 
  scale_color_manual(values = c("#DDA7AB", "#4A5899")) +
  theme_bw()
ggsave(p, filename = "tajimas_D.pdf", device = "pdf")
