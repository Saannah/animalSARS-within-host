### load table with richness and shannon data
tab = read.delim("/Users/sana/Documents/AIM3/011-freyja/richness_and_shannon.tsv")
setwd("/Users/sana/Documents/AIM3/012-plot_richness_and_shannon")
library(ggplot2)
library(tidyverse)
#plot number of segregating sites
p = ggplot(tab, aes(x=species, y=shannon_diversity, fill = lab_infected)) + 
  geom_jitter(aes(color = lab_infected), width = 0.1, size = 1) + 
  scale_color_manual(values = c("#DDA7AB", "#4A5899")) +
  scale_shape_manual(values = c(17, 8)) +
  theme_bw()
ggsave(p, filename = "shannon_diversity.pdf", device = "pdf", height = 5, width = 5)


p = ggplot(tab, aes(x=species, y=richness, fill = lab_infected)) + 
  geom_jitter(aes(color = lab_infected), height = 0.2, size = 1) + 
  scale_color_manual(values = c("#DDA7AB", "#4A5899")) +
  scale_shape_manual(values = c(17, 8)) +
  theme_bw()
ggsave(p, filename = "richness.pdf", device = "pdf", height = 5, width = 5)
