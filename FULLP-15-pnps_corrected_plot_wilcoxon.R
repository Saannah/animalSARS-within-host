######## load and aggregate both genomewide and S gene data
rm(list = ls())
### load metadata 
#metadata = read.delim("/users/sana/Documents/AIM3/021-full_pipeline/00-ALL_FILES/matched_metadata.tsv")
set.seed(1)

pnps_gw = read.delim("/users/sana/Documents/AIM3/021-full_pipeline/10-correction_pnps/all_pnps_corrected.tsv")
pnps_gw$Run = pnps_gw$SRA_id
pnps_gw$pNpS_genomewide = pnps_gw$pNpS

pnps_sgene = read.delim("/users/sana/Documents/AIM3/021-full_pipeline/10-correction_pnps/all_pnps_corrected_sgene.tsv")
pnps_sgene$Run = pnps_sgene$SRA_id


all_pnps = merge(pnps_gw[, c("Run", "pNpS_genomewide")],
                 pnps_sgene[, c("Run", "pNpS_sgene")],
                 "Run")



### filter based on sample quality:
all_data = read.delim("/users/sana/Documents/AIM3/021-full_pipeline/1-metrics/metrics_strict_thresholds_min50_5percent.tsv")
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
all_data$Infection = all_data$lab_infected

### merge pnps with the rest of the data
all_data = merge(all_data, all_pnps, "Run")


s=22
#### plot quality filtered pnps data points
p_gw = ggplot(all_data, aes(x = Species, 
                     y=pNpS_genomewide)) + 
  geom_boxplot(aes(x=Species, y=pNpS_genomewide), fill = "#4A5899", alpha = 0.1, outlier.shape = NA) +
  geom_jitter(aes(color = Infection,
                  shape = tissue_type),
              width = 0.25, height=0, size = 1, stroke = 1) + 
  scale_color_manual(values = c("#4A5899", "#CD5D67")) + 
  scale_shape_manual(values = c(19, 5)) +
  ylim(0,5.5) +
  theme_bw() + xlab("Species") + ggtitle("a. Genome-wide") +
  ylab("pN/pS") + theme(legend.position = "none",
                                              text = element_text(size=s),
                                              plot.title = element_text(size = 15, face = "bold"))

p_sgene = ggplot(all_data, aes(x = Species, 
                            y=pNpS_sgene)) + 
  geom_boxplot(aes(x=Species, y=pNpS_sgene), fill = "#4A5899", alpha = 0.1, outlier.shape = NA) +
  geom_jitter(aes(color = Infection,
                  shape = tissue_type),
              width = 0.25, height=0, size = 1, stroke = 1) + 
  scale_color_manual(values = c("#4A5899", "#CD5D67")) + 
  scale_shape_manual(values = c(19, 5)) +
  ylim(0,5.5) +
  theme_bw() + xlab("Species") + ggtitle("b. Spike gene") +
  ylab("pN/pS") + theme(legend.position = "none",
                       text = element_text(size=s),
                       plot.title = element_text(size = 15, face = "bold"))

both = ggarrange(p_gw, p_sgene, nrow =2, ncol = 1)
ggsave(both, 
       file = "/users/sana/Documents/AIM3/021-full_pipeline/10-correction_pnps/figure3_alldatapoints.pdf", 
       device = "pdf",
       height = 7,
       width = 12)




##### remove lymph node tissue and re-plot
all_data = all_data[which(all_data$tissue_type == "swab"), ]
p_gw = ggplot(all_data, aes(x = Species, 
                            y=pNpS_genomewide)) + 
  geom_boxplot(aes(x=Species, y=pNpS_genomewide), fill = "#4A5899", alpha = 0.1, outlier.shape = NA) +
  geom_jitter(aes(color = Infection,
                  shape = tissue_type),
              width = 0.25, height=0, size = 1, stroke = 1) + 
  scale_color_manual(values = c("#4A5899", "#CD5D67")) + 
  scale_shape_manual(values = c(19, 5)) +
  ylim(0,3) +
  theme_bw() + xlab("Species") + ggtitle("a. pNpS genomewide (swab only samples)") +
  ylab("pNpS") + theme(legend.position = "none",
                       text = element_text(size=s),
                       plot.title = element_text(size = 15, face = "bold"))

p_sgene = ggplot(all_data, aes(x = Species, 
                               y=pNpS_sgene)) + 
  geom_boxplot(aes(x=Species, y=pNpS_sgene), fill = "#4A5899", alpha = 0.1, outlier.shape = NA) +
  geom_jitter(aes(color = Infection,
                  shape = tissue_type),
              width = 0.25, height=0, size = 1, stroke = 1) + 
  scale_color_manual(values = c("#4A5899", "#CD5D67")) + 
  scale_shape_manual(values = c(19, 5)) +
  ylim(0,3) +
  theme_bw() + xlab("Species") + ggtitle("b. pNpS S gene (swab only samples)") +
  ylab("pNpS") + theme(legend.position = "none",
                       text = element_text(size=s),
                       plot.title = element_text(size = 15, face = "bold"))

both = ggarrange(p_gw, p_sgene, nrow =2, ncol = 1)
ggsave(both, 
       file = "/users/sana/Documents/AIM3/021-full_pipeline/10-correction_pnps/figure3_swabonly.pdf", 
       device = "pdf",
       height = 7,
       width = 12)


