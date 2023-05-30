#### wilcoxon test
wilcoxon = data.frame(species1 = vector(), species2 = vector(), p_value_n_seg = vector(), p_value_coverage = vector(),
                      p_value_pi = vector(), p_value_tajimasD = vector())
wilcoxon[1,] = c("mink", "deer", -1, -1, -1, -1)
wilcoxon[2,] = c("mink", "cat", -1, -1, -1, -1)
wilcoxon[3,] = c("mink", "dog", -1, -1, -1, -1)
wilcoxon[4,] = c("deer", "cat", -1, -1, -1, -1)
wilcoxon[5,] = c("deer", "dog", -1, -1, -1, -1)
wilcoxon[6,] = c("dog", "cat", -1, -1, -1, -1)

metrics$pi = as.numeric(metrics$pi)
metrics$tajimas_D = as.numeric(metrics$tajimas_D)

for(i in 1:6){
  species1 = wilcoxon$species1[i]
  species2 = wilcoxon$species2[i]
  wilcoxon$p_value_n_seg[i] = wilcox.test(metrics$n_segregating_sites[which(metrics$species == species1)], metrics$n_segregating_sites[which(metrics$species == species2)])$p.value
  wilcoxon$p_value_coverage[i] = wilcox.test(metrics$vcf_avg_depth_after_filtering_sites[which(metrics$species == species1)], metrics$vcf_avg_depth_after_filtering_sites[which(metrics$species == species2)])$p.value
  wilcoxon$p_value_pi[i] = wilcox.test(metrics$pi[which(metrics$species == species1)], metrics$pi[which(metrics$species == species2)])$p.value
  wilcoxon$p_value_tajimasD[i] = wilcox.test(metrics$tajimas_D[which(metrics$species == species1)], metrics$tajimas_D[which(metrics$species == species2)])$p.value
  
}

write.table(wilcoxon, file = "wilcox.tsv", sep = "\t", row.names = FALSE)
