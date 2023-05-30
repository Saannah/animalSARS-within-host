#### wilcoxon test
wilcoxon = data.frame(species1 = vector(), species2 = vector(), p_value_n_seg = vector(),
                      p_value_pi = vector())
wilcoxon[1,] = c("human", "cat", -1, -1)
wilcoxon[2,] = c("human", "dog", -1, -1)
wilcoxon[3,] = c("human", "mink", -1, -1)
wilcoxon[4,] = c("human", "deer", -1, -1)

#load all_data from script 018
#create pi_data
pi_data = all_data[-which(all_data$pi == "-"), ]
pi_data$pi = as.numeric(pi_data$pi)

for(i in 1:4){
  species1 = wilcoxon$species1[i]
  species2 = wilcoxon$species2[i]
  wilcoxon$p_value_n_seg[i] = wilcox.test(all_data$n_segregating_sites[which(all_data$Species == species1)], all_data$n_segregating_sites[which(all_data$Species == species2)])$p.value
  #wilcoxon$p_value_coverage[i] = wilcox.test(metrics$vcf_avg_depth_after_filtering_sites[which(metrics$species == species1)], metrics$vcf_avg_depth_after_filtering_sites[which(metrics$species == species2)])$p.value
  wilcoxon$p_value_pi[i] = wilcox.test(pi_data$pi[which(pi_data$Species == species1)], pi_data$pi[which(pi_data$Species == species2)])$p.value
  #wilcoxon$p_value_tajimasD[i] = wilcox.test(metrics$tajimas_D[which(metrics$species == species1)], metrics$tajimas_D[which(metrics$species == species2)])$p.value
  
}

write.table(wilcoxon, file = "wilcox.tsv", sep = "\t", row.names = FALSE)
