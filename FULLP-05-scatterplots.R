## scatter plots with the new data
rm(list = ls())
library(ggplot2)
library(ggpubr)
all_data = read.delim("/users/sana/Documents/AIM3/021-full_pipeline/1-metrics/metrics_strict_thresholds_min50_5percent.tsv")


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
all_data$brdth = all_data$min_50_site_count/29930
all_data$filtered_richness = as.numeric(all_data$filtered_richness)
all_data$filtered_shannon = as.numeric(all_data$filtered_shannon)

all_data$Infection = all_data$lab_infected

### load tonkins pi data
tpi = read.delim("/users/sana/Documents/AIM3/021-full_pipeline/1-metrics/metrics_w_tonkins_pi.tsv")
typeof(tpi$pi_tonkin)
tpi = tpi[, c("sample_id", "pi_tonkin")]
all_data = merge(all_data, tpi, "sample_id")

s=22



#### create plots
n_seg_plot = ggplot(all_data, aes(x=factor(Species, level = c("human", "cat", "dog", "mink", "deer")), 
                                  y=n_seg)) + 
  geom_boxplot(aes(x=Species, y=n_seg), fill = "#4A5899", alpha = 0.1, outlier.shape = NA) +
  geom_jitter(aes(color = Infection,
                  shape = tissue_type),
              width = 0.25, height=0, size = 1, stroke = 1) + 
  scale_color_manual(values = c("#4A5899", "#CD5D67")) + 
  scale_shape_manual(values = c(19, 5)) +
  ylim(0,50) +
  theme_bw() + xlab("Species") + ggtitle("a.") +
  ylab("Number of iSNVs") + theme(legend.position = "none",
                                              text = element_text(size=s),
                                              plot.title = element_text(size = 15, face = "bold"))

pi_plot = ggplot(all_data, aes(x=factor(Species, level = c("human", "cat", "dog", "mink", "deer")), 
                               y=pi_tonkin)) + 
  geom_boxplot(aes(x=Species, y=pi_tonkin), fill = "#4A5899", alpha = 0.1, outlier.shape = NA) +
  geom_jitter(aes(color = Infection,
                  shape = tissue_type),
              width = 0.25, height=0, size = 1, stroke = 1) + 
  scale_color_manual(values = c("#4A5899", "#CD5D67")) + 
  scale_shape_manual(values = c(19, 5)) +
  ylim(0,20) +
  theme_bw() + xlab("Species") + ggtitle("b.") +
  ylab("Average number of pairwise differences  ") + theme(legend.position = "none",
                                                         text = element_text(size=s),
                                                         plot.title = element_text(size = 15, face = "bold"))


### scale tajima's D to mean 0 and variance 1
all_data$tajimasd_scaled = scale(all_data$tajimasd)

tajimasd_plot = ggplot(all_data, aes(x=factor(Species, level = c("human", "cat", "dog", "mink", "deer")), 
                                     y=tajimasd_scaled)) + 
  geom_boxplot(aes(x=Species, y=tajimasd_scaled), fill = "#4A5899", alpha = 0.1, outlier.shape = NA) +
  geom_jitter(aes(color = Infection,
                  shape = tissue_type),
              width = 0.25, height=0, size = 1, stroke = 1) + 
  scale_color_manual(values = c("#4A5899", "#CD5D67"), name = "Infection type", labels = c("Natural infection", "Lab infection")) + 
  scale_shape_manual(values = c(19, 5), name = "Tissue type", labels = c("Nasopharyngeal swab", "Lymph node tissue")) +
  ylim(-2,3) +
  theme_bw() + xlab("Species") + ggtitle("c. Tajima's D") +
  ylab("Tajima's D") + theme(legend.position = "none",
                             text = element_text(size=s),
                             plot.title = element_text(size = 15, face = "bold"))

#ggsave(tajimasd_plot, file = "/users/sana/Documents/AIM3/021-full_pipeline/xx-Manuscript/Figure1and2_legend.pdf", device = "pdf", height = 3, width =5)

figure1 = ggarrange(n_seg_plot, pi_plot,
                    nrow =1, ncol = 2, widths = c(6,6),
                    heights = c(10,10)) 
ggsave(figure1, file = "/users/sana/Documents/AIM3/021-full_pipeline/4-scatterplots/figure1_strict_thresholds.pdf",
       device = "pdf", height = 7, width = 14)




richness_plot = ggplot(all_data, aes(x=factor(Species, level = c("human", "cat", "dog", "mink", "deer")), 
                                     y=filtered_richness)) + 
  geom_boxplot(aes(x=Species, y=filtered_richness), fill = "#4A5899", alpha = 0.1, outlier.shape = NA) +
  geom_jitter(aes(color = Infection,
                  shape = tissue_type),
              width = 0.25, height=0, size = 1, stroke = 1) + 
  scale_color_manual(values = c("#4A5899", "#CD5D67")) + 
  scale_shape_manual(values = c(19, 5)) +
  ylim(0,25) +
  theme_bw() + xlab("Species") + ggtitle("a.") +
  ylab("Number of lineages present") + theme(legend.position = "none",
                                             text = element_text(size=s),
                                             plot.title = element_text(size = 15, face = "bold"))

shannon_plot = ggplot(all_data, aes(x=factor(Species, level = c("human", "cat", "dog", "mink", "deer")), 
                                    y=filtered_shannon)) + 
  geom_boxplot(aes(x=Species, y=filtered_shannon), fill = "#4A5899", alpha = 0.1, outlier.shape = NA) +
  geom_jitter(aes(color = Infection,
                  shape = tissue_type),
              width = 0.25, height=0, size = 1, stroke = 1) + 
  scale_color_manual(values = c("#4A5899", "#CD5D67")) + 
  scale_shape_manual(values = c(19, 5)) +
  ylim(0,3.5) +
  theme_bw() + xlab("Species") + ggtitle("b.") +
  ylab("Shannon diversity index") + theme(legend.position = "none",
                                          text = element_text(size=s),
                                          plot.title = element_text(size = 15, face = "bold"))


figure2 = ggarrange(richness_plot, shannon_plot, nrow =1, ncol = 2) 
ggsave(figure2, file = "/users/sana/Documents/AIM3/021-full_pipeline/4-scatterplots/figure2.pdf",
       device = "pdf", height = 7, width = 14)


# all = ggarrange(n_seg_plot, pi_plot, richness_plot, shannon_plot, ncol = 4, nrow = 1, widths=c(5,5.5,5,5))
# ggsave(all, file = "/users/sana/Documents/AIM3/021-full_pipeline/4-scatterplots/all.pdf", device = "pdf", width = 16, height = 9)





### get legend for the plots
legend_plot = ggplot(all_data, aes(x=factor(Species, level = c("human", "cat", "dog", "mink", "deer")), 
                                    y=filtered_shannon)) + 
  geom_boxplot(aes(x=Species, y=filtered_shannon), fill = "#4A5899", alpha = 0.1, outlier.shape = NA) +
  geom_jitter(aes(color = Infection,
                  shape = tissue_type),
              width = 0.25, height=0, size = 1, stroke = 1) + 
  scale_color_manual(values = c("#4A5899", "#CD5D67")) + 
  scale_shape_manual(values = c(19, 5)) +
  ylim(0,3.5) +
  theme_bw() + xlab("Species") + ggtitle("b. Shannon diversity index") +
  ylab("Shannon diversity index") + theme(legend.position = "right",
                                          text = element_text(size=s),
                                          plot.title = element_text(size = 15, face = "bold"),
                                          legend.text=element_text(size=25)) +
  labs(col="Type of infection", shape="Tissue type")


library(ggpubr)
leg = get_legend(legend_plot)
ggsave(leg, file = "/users/sana/Documents/AIM3/021-full_pipeline/xx-Manuscript/Figure1and2_legend_large.pdf", device = "pdf")


