### script to match calculated stats to the metadata downloaded from NCBI
library(stringr)
#load the metrics dataframe
metrics = read.delim("/Users/sana/Documents/AIM3/006-all_runs_variant_files/metrics_df_forparsed.tsv")
# n = read.delim("/Users/sana/Dropbox/AIM3/variant_plots/NoSegSites_names_list.txt")
# 
# #load metadata for all 4 species
# setwd("/Users/sana/Documents/AIM3/004-metadata_categorization/paired_single_new")
# cat = read.delim("cat_final.tsv")
# dog = read.delim("dog_final.tsv")
# deer = read.delim("deer_final.tsv")
# mink = read.delim("mink_final.tsv")
# all_metadata = rbind(cat[, 1:3], dog[, 1:3], mink[, 1:3], deer[, 1:3])
# 
# #remove data that do not exist in metadata and vice versa
# metrics = metrics[which(metrics$SRA_ID %in% all_metadata$Run), ]
# all_metadata = all_metadata[which(all_metadata$Run %in% metrics$SRA_ID), ]
# 
# 
# #add a column to metrics dataframe and specify species
# metrics$species = "TBD"
# metrics$species[which(metrics$SRA_ID %in% cat$Run)] = "cat"
# metrics$species[which(metrics$SRA_ID %in% dog$Run)] = "dog"
# metrics$species[which(metrics$SRA_ID %in% mink$Run)] = "mink"
# metrics$species[which(metrics$SRA_ID %in% deer$Run)] = "deer"
# 
# #add a columnn to metrics dataframe and specify lab infected or not
# metrics$lab_infected = "TBD"
# metrics$lab_infected[match(all_metadata$Run, metrics$SRA_ID)] = all_metadata$lab_infected
# 
# 
# #add a column to specify runtype
# readset = read.delim("/Users/sana/Documents/AIM3/006-all_runs_variant_files/readset_for_paired_and_single.tsv")
# readset = readset[which(readset$Sample %in% metrics$SRA_ID), ]
# metrics$runtype = "TBD"
# metrics$runtype[match(readset$Sample, metrics$SRA_ID)] = readset$RunType
# 
# 
# metrics$sample_average_depth = as.numeric(metrics$sample_average_depth)
# metrics$n_segregating_sites = as.numeric(metrics$n_segregating_sites)


library(ggplot2)
library(tidyverse)
#plot number of segregating sites
p = ggplot(metrics, aes(x=species, y=n_segregating_sites, fill = lab_infected, shape=runtype)) + 
  geom_jitter(aes(color = lab_infected), width = 0.2, size = 1) + 
  scale_color_manual(values = c("#DDA7AB", "#4A5899")) +
  scale_shape_manual(values = c(17, 8)) + ylim(0, 100) +
  theme_bw()
ggsave(p, filename = "n_segsites.pdf", device = "pdf")

metrics$vcf_avg_depth_after_filtering_sites = as.numeric(
  metrics$vcf_avg_depth_after_filtering_sites
)
#plot average depth
p = ggplot(metrics, aes(x=species, y=vcf_avg_depth_after_filtering_sites, fill = lab_infected, shape = runtype)) + 
  geom_jitter(aes(color = lab_infected), width = 0.2, size = 1) + 
  scale_color_manual(values = c("#DDA7AB", "#4A5899")) +
  scale_shape_manual(values = c(17, 8)) + ylim(0, 15000) +
  theme_bw()
ggsave(p, filename = "average_depth.pdf", device = "pdf")

numeric_dataframe = metrics
numeric_dataframe$pi = as.numeric(numeric_dataframe$pi)
numeric_dataframe$tajimas_D = as.numeric(numeric_dataframe$tajimas_D)

#plot pi
p = ggplot(numeric_dataframe, aes(x=species, y=pi, fill = lab_infected, shape = runtype)) + 
  geom_jitter(aes(color = lab_infected), width = 0.2, size = 1) + 
  scale_color_manual(values = c("#DDA7AB", "#4A5899")) +
  scale_shape_manual(values = c(17, 8)) + ylim(0, 750) +
  theme_bw()
ggsave(p, filename = "pi.pdf", device = "pdf")

#plot tajimas D
p = ggplot(numeric_dataframe, aes(x=species, y=tajimas_D, fill = lab_infected, shape = runtype)) + 
  geom_jitter(aes(color = lab_infected), width = 0.2, size = 1) + 
  scale_color_manual(values = c("#DDA7AB", "#4A5899")) +
  scale_shape_manual(values = c(17, 8)) + ylim(0, 10) +
  theme_bw()
ggsave(p, filename = "tajimas_D.pdf", device = "pdf")


# 
# ##### n_seg_sites ouyliers:
# sra = metrics$SRA_ID[which(metrics$n_segregating_sites > 200)]
# all_metadata[which(all_metadata$Run %in% sra),]
# mink_out = mink[which(mink$Run %in% sra), ]
# 
# ##### average depth outliers
# sra = metrics$SRA_ID[which(metrics$average_depth_after_filtering_sites < 200)]
# mink_out = mink[which(mink$Run %in% sra), ]
# 
# #### mink from netherlands 
# 
