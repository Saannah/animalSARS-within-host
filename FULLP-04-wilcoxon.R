### wilcoxon test for the metrics nseg, pi, shannon, richness, and tajima's D

rm(list = ls())
library(car)
library(MASS)
library(lme4)




## load previous metric dataset
all_data = read.delim("/users/sana/Documents/AIM3/021-full_pipeline/1-metrics/metrics_strict_thresholds_min50_5percent.tsv")

## load tonkins pi data
metrics_w_tonkin = read.delim("/users/sana/Documents/AIM3/021-full_pipeline/1-metrics/metrics_w_tonkins_pi.tsv")
pi_tonkin = metrics_w_tonkin[, c("sample_id", "pi_tonkin")]

## merge tonkins pi into all_data
all_data = merge(all_data, pi_tonkin, "sample_id")

#remove sites that have lower than 50 mean depth
all_data = all_data[-which(all_data$mean_depth_from_bed < 50), ]
all_data = all_data[-which(all_data$min_50_site_count < 15000), ]


all_data$Species = relevel(as.factor(all_data$Species), ref = "human")

all_data$tissue_type = "swab"
all_data$tissue_type[which(all_data$tissue == "lymph-tissue")] = "lymph"
all_data$tissue_type = relevel(as.factor(all_data$tissue_type), ref = "swab")

all_data$lab_infected = relevel(as.factor(all_data$lab_infected), ref = "natural")

all_data$tissue = relevel(as.factor(all_data$tissue), ref = "nasal/pharyngeal-swab")
all_data$assay_type = relevel(as.factor(all_data$assay_type), ref = "AMPLICON")
all_data$BioProject = as.factor(all_data$BioProject)
all_data$shannon = as.numeric(all_data$filtered_shannon)
all_data$dpth = as.vector(scale(all_data$mean_depth_from_bed))
all_data$brdth = all_data$min_50_site_count/29903
all_data$filtered_richness = as.numeric(all_data$filtered_richness)
all_data$filtered_shannon = as.numeric(all_data$filtered_shannon)

### scale tajima's D
#all_data$tajimasd = scale(all_data$tajimasd)

############################################################################################################################
################################################# wilcoxon test

wilcoxon_df = data.frame(species = vector(),
                         nseg_pvalue = vector(),
                         pi_pvalue = vector(),
                         shannon_pvalue = vector(),
                         richness_pvalue = vector(),
                         tajimasd_pvalue = vector())


for(species in c("cat", "dog", "mink", "deer")){
  #nseg test
  w_nseg = wilcox.test(all_data$n_seg[which(all_data$Species == "human")],
                       all_data$n_seg[which(all_data$Species == species)])
  #pi test
  w_pi = wilcox.test(all_data$pi_tonkin[which(all_data$Species == "human")],
                     all_data$pi_tonkin[which(all_data$Species == species)])
  #shannon test
  w_shannon = wilcox.test(all_data$filtered_shannon[which(all_data$Species == "human")],
                          all_data$filtered_shannon[which(all_data$Species == species)])
  #richness test
  w_richness = wilcox.test(all_data$filtered_richness[which(all_data$Species == "human")],
                           all_data$filtered_richness[which(all_data$Species == species)])
  
  w_tajimasd = wilcox.test(all_data$tajimasd[which(all_data$Species == "human")],
                           all_data$tajimasd[which(all_data$Species == species)])
  
  wilcoxon_df[nrow(wilcoxon_df) + 1, ] = c(paste0(species, " vs human"),
                                           w_nseg$p.value,
                                           w_pi$p.value,
                                           w_shannon$p.value,
                                           w_richness$p.value,
                                           w_tajimasd$p.value)
  
}


write.table(wilcoxon_df, "/users/sana/Documents/AIM3/021-full_pipeline/3-wilcoxon/wilcoxon_strictcutoff.tsv", sep = "\t", row.names = FALSE)



#### repeat wilcoxon test after removal of lymph node samples
## subset dataframe into nasopharyngeal samples only
all_data = all_data[-which(all_data$tissue_type == "lymph"), ]

wilcoxon_df = data.frame(species = vector(),
                         nseg_pvalue = vector(),
                         pi_pvalue = vector(),
                         shannon_pvalue = vector(),
                         richness_pvalue = vector(),
                         tajimasd_pvalue = vector())


for(species in c("cat", "dog", "mink", "deer")){
  #nseg test
  w_nseg = wilcox.test(all_data$n_seg[which(all_data$Species == "human")],
                       all_data$n_seg[which(all_data$Species == species)])
  #pi test
  w_pi = wilcox.test(all_data$pi_tonkin[which(all_data$Species == "human")],
                     all_data$pi_tonkin[which(all_data$Species == species)])
  #shannon test
  w_shannon = wilcox.test(all_data$filtered_shannon[which(all_data$Species == "human")],
                          all_data$filtered_shannon[which(all_data$Species == species)])
  #richness test
  w_richness = wilcox.test(all_data$filtered_richness[which(all_data$Species == "human")],
                           all_data$filtered_richness[which(all_data$Species == species)])
  
  w_tajimasd = wilcox.test(all_data$tajimasd[which(all_data$Species == "human")],
                           all_data$tajimasd[which(all_data$Species == species)])
  
  wilcoxon_df[nrow(wilcoxon_df) + 1, ] = c(paste0(species, " vs human"),
                                           w_nseg$p.value,
                                           w_pi$p.value,
                                           w_shannon$p.value,
                                           w_richness$p.value,
                                           w_tajimasd$p.value)
  
}












#### get the median and mean values for each of the metrics:
## MAKE SURE THE TABLE GOING INTO THIS SECTION IS THE TABLE WITH ALL ENTRIES
rm(list = ls())
mean_table = data.frame(species = vector(),
                         nseg = vector(),
                         pi = vector(),
                         filtered_shannon = vector(),
                         filtered_richness = vector()
                         )


median_table = data.frame(species = vector(),
                        nseg = vector(),
                        pi = vector(),
                        filtered_shannon = vector(),
                        filtered_richness = vector()
)

for(species in c("cat", "dog", "mink", "deer", "human")){
  mean_table = rbind(mean_table,
                     data.frame(species = species,
                                nseg = mean(all_data$n_seg[which(all_data$Species == species)]),
                                pi = mean(all_data$pi[which(all_data$Species == species)]),
                                filtered_shannon = mean(all_data$filtered_shannon[which(all_data$Species == species)]),
                                filtered_richness = mean(all_data$filtered_richness[which(all_data$Species == species)]),
                                tajimasd = mean(all_data$tajimasd[which(all_data$Species == species)])
                     )
                     )
  
  
  median_table = rbind(median_table,
                     data.frame(species = species,
                                nseg = median(all_data$n_seg[which(all_data$Species == species)]),
                                pi = median(all_data$pi[which(all_data$Species == species)]),
                                filtered_shannon = median(all_data$filtered_shannon[which(all_data$Species == species)]),
                                filtered_richness = median(all_data$filtered_richness[which(all_data$Species == species)]),
                                tajimasd = median(all_data$tajimasd[which(all_data$Species == species)])
                     )
  )
  
}
write.table(mean_table, file = "/users/sana/Documents/AIM3/021-full_pipeline/xx-Manuscript/mean_table_strictcutoff_filtered.tsv", sep = "\t", row.names = FALSE)
write.table(median_table, file = "/users/sana/Documents/AIM3/021-full_pipeline/xx-Manuscript/median_table_strictcutoff_filtered.tsv", sep = "\t", row.names = FALSE)




### deer df
deer_data = all_data[which(all_data$Species == "deer"), ]
wilcox.test(deer_data$filtered_richness[which(deer_data$tissue_type == "lymph")],
            deer_data$filtered_richness[which(deer_data$tissue_type == "swab")])
mean(deer_data$filtered_richness[which(deer_data$tissue_type == "lymph")])

mean(deer_data$filtered_richness[which(deer_data$tissue_type == "swab")])

