### prepare supplementary table 2
## scatter plots with the new data
rm(list = ls())
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
all_data$brdth = all_data$min_50_site_count/29903
all_data$filtered_richness = as.numeric(all_data$filtered_richness)
all_data$filtered_shannon = as.numeric(all_data$filtered_shannon)

all_data$Infection = all_data$lab_infected

### load tonkins pi data
tpi = read.delim("/users/sana/Documents/AIM3/021-full_pipeline/1-metrics/metrics_w_tonkins_pi.tsv")
typeof(tpi$pi_tonkin)
tpi = tpi[, c("sample_id", "pi_tonkin")]
all_data = merge(all_data, tpi, "sample_id")


pnps_gw = read.delim("/users/sana/Documents/AIM3/021-full_pipeline/10-correction_pnps/all_pnps_corrected.tsv")
pnps_gw$Run = pnps_gw$SRA_id
pnps_gw$pNpS_genomewide = pnps_gw$pNpS

pnps_sgene = read.delim("/users/sana/Documents/AIM3/021-full_pipeline/10-correction_pnps/all_pnps_corrected_sgene.tsv")
pnps_sgene$Run = pnps_sgene$SRA_id


all_pnps = merge(pnps_gw[, c("Run", "pNpS_genomewide")],
                 pnps_sgene[, c("Run", "pNpS_sgene")],
                 "Run")

all_data = merge(all_data, all_pnps, "Run")

supp2 = all_data[, c("Run", "Species" ,"n_seg", "pi_tonkin", "filtered_shannon", "filtered_richness", "pNpS_genomewide", "pNpS_sgene","dpth", "brdth",
                     "tissue_type", "Infection", "Collection_Date", "geo_loc_name_country", "assay_type")]

write.table(supp2, file = "/users/sana/Documents/AIM3/021-full_pipeline/xx-Manuscript/supp tables/Supplementary_Table_2_updated_and_corrected.tsv")


