### script to match calculated stats to the metadata downloaded from NCBI
library(stringr)
#### load coverage data
setwd("/Users/sana/Documents/AIM3/007-all_runs_depth_summary_files")
names = read.delim("names.tsv", header = FALSE)
names = names[,1]
names = names[-c(717, 718)]

coverage_df = data.frame(id = vector(), mean_coverage = vector())

for(i in names){
  t = read.delim(i)
  coverage_df[nrow(coverage_df) + 1, ] = c(t$sample_id[1], t$mean[1])
}


#load the metrics dataframe
metrics = read.delim("/Users/sana/Documents/AIM3/pres_plots/dataframe.tsv")
n = read.delim("/Users/sana/Documents/AIM3/variant_plots-old/NoSegSites_names_list.txt")
n = n[, 1]
nosegsites_names_list = vector()
for(i in 1:length(n)){
  nosegsites_names_list = c(nosegsites_names_list, unlist(str_split(n, "\\.")[i])[1])
}

#load metadata for all 4 species
setwd("/Users/sana/Documents/AIM3/004-metadata_categorization/")
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
  metadata = read.delim(paste0("/Users/sana/Documents/AIM3/004-metadata_categorization/", s, "_final.tsv"))
  metrics$lab_infected[match(metadata[, 1], metrics$SRA_ID)] = metadata[, 2] 
}
metrics = metrics[which(metrics$SRA_ID %in% all_metadata$Run), ]
metrics$lab_infected[which(metrics$lab_infected == "probably no")] = "no"
metrics$lab_infected[which(metrics$lab_infected == "no")] = "Natural Infection"
metrics$lab_infected[which(metrics$lab_infected == "yes")] = "Lab Infected"
metrics = metrics[-which(metrics$lab_infected == "Missing data"), ]


metrics = metrics[which(metrics$enough_avg_depth == "y"),]

##### this sectiona dds the two samples that had no segsites returned but had enough coverage
noseg_coverage_df = coverage_df[which(coverage_df$id %in% nosegsites_names_list),]
to_add_zeros = noseg_coverage_df[which(as.numeric(noseg_coverage_df$mean_coverage) > 50),]

c(to_add_zeros[1,1], "y", to_add_zeros[1,2], 0, "-", "-", "TBD", "TBD")
c(to_add_zeros[2,1], "y", to_add_zeros[2,2], 0, "-", "-", "TBD", "TBD")

## these two rows are samples that had no seg sites but > 50 x coverage
metrics[215,] = c(to_add_zeros[1,1], "y", to_add_zeros[1,2], 0, "-", "-", "mink", "Natural Infection")
metrics[216,] = c(to_add_zeros[2,1], "y", to_add_zeros[2,2], 0, "-", "-", "mink", "Natural Infection")




library(ggplot)
library(tidyverse)
setwd("/Users/sana/Documents/AIM3/pres_plots/")
#plot number of segregating sites
metrics$n_segregating_sites = as.numeric(metrics$n_segregating_sites)
p = ggplot(metrics, aes(x=species, y=n_segregating_sites, fill = lab_infected)) + 
  geom_jitter(aes(color = lab_infected), width = 0.09, size = 1) + 
  scale_color_manual(values = c("#DDA7AB", "#4A5899")) +
  theme_bw()
ggsave(p, filename = "n_segsites.pdf", device = "pdf")


metrics$average_depth = as.numeric(coverage_df$mean_coverage[which(coverage_df$id %in% metrics$SRA_ID)])
#plot average depth
p = ggplot(metrics, aes(x=species, y=average_depth, fill = lab_infected)) + 
  geom_jitter(aes(color = lab_infected), width = 0.09, size = 1) + 
  scale_color_manual(values = c("#DDA7AB", "#4A5899")) +
  theme_bw()
ggsave(p, filename = "average_depth.pdf", device = "pdf")

numeric_dataframe = metrics
numeric_dataframe$pi = as.numeric(numeric_dataframe$pi)
numeric_dataframe$tajimas_D = as.numeric(numeric_dataframe$tajimas_D)

metrics$pi = as.numeric(metrics$pi)
#plot pi
p = ggplot(numeric_dataframe, aes(x=species, y=pi, fill = lab_infected)) + 
  geom_jitter(aes(color = lab_infected), width = 0.09, size = 1) + 
  scale_color_manual(values = c("#DDA7AB", "#4A5899")) +
  theme_bw()
ggsave(p, filename = "pi.pdf", device = "pdf")

metrics$tajimas_D = as.numeric(metrics$tajimas_D)
#plot tajimas D
p = ggplot(numeric_dataframe, aes(x=species, y=tajimas_D, fill = lab_infected)) + 
  geom_jitter(aes(color = lab_infected), width = 0.09, size = 1) + 
  scale_color_manual(values = c("#DDA7AB", "#4A5899")) +
  theme_bw()
ggsave(p, filename = "tajimas_D.pdf", device = "pdf")










### add no seg site files
## create table with coverage summary
save_metrics = metrics

coverage_df$mean_coverage[which(coverage_df$id %in% metrics$SRA_ID)]


#### wilcoxon test
wilcoxon = data.frame(species1 = vector(), species2 = vector(), p_value_n_seg = vector(), p_value_coverage = vector(),
                      p_value_pi = vector(), p_value_tajimasD = vector())
wilcoxon[1,] = c("mink", "deer", -1, -1, -1, -1)
wilcoxon[2,] = c("mink", "cat", -1, -1, -1, -1)
wilcoxon[3,] = c("mink", "dog", -1, -1, -1, -1)
wilcoxon[4,] = c("deer", "cat", -1, -1, -1, -1)
wilcoxon[5,] = c("deer", "dog", -1, -1, -1, -1)
wilcoxon[6,] = c("dog", "cat", -1, -1, -1, -1)

for(i in 1:6){
  species1 = wilcoxon$species1[i]
  species2 = wilcoxon$species2[i]
  wilcoxon$p_value_n_seg[i] = wilcox.test(metrics$n_segregating_sites[which(metrics$species == species1)], metrics$n_segregating_sites[which(metrics$species == species2)])$p.value
  wilcoxon$p_value_coverage[i] = wilcox.test(metrics$average_depth[which(metrics$species == species1)], metrics$average_depth[which(metrics$species == species2)])$p.value
  wilcoxon$p_value_pi[i] = wilcox.test(metrics$pi[which(metrics$species == species1)], metrics$pi[which(metrics$species == species2)])$p.value
  wilcoxon$p_value_tajimasD[i] = wilcox.test(metrics$tajimas_D[which(metrics$species == species1)], metrics$tajimas_D[which(metrics$species == species2)])$p.value
  
  }

write.table(wilcoxon, file = "wilcox.tsv", sep = "\t", row.names = FALSE)
