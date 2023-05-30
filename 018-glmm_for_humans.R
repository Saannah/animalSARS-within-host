#load all the metadata
all_metadata = read.delim("/users/sana/Documents/AIM3/016-glmm_w_human_ref/metadata_all_w_human_cleaned.tsv")


### load the nosegsites names
human_no_seg = read.delim("/users/sana/Documents/AIM3/006-all_runs_variant_files/human_runs/nosegnames", header = FALSE)
human_no_seg = strsplit(human_no_seg[, 1], "\\.")

animal_no_seg = read.delim("/users/sana/Documents/AIM3/006-all_runs_variant_files/0-paired_end/nosegnames", header = FALSE)
animal_no_seg = strsplit(animal_no_seg[, 1], "\\.")

noseg = c(human_no_seg, animal_no_seg)
no_seg_calls = vector()
for(i in 1:length(noseg)){
  no_seg_calls = c(no_seg_calls, unlist(noseg[i])[1])
}
metric_df_for_no_seg = data.frame(SRA_ID = no_seg_calls,
                                  enough_avg_depth = "y",
                                  average_depth = "-",
                                  n_segregating_sites = 0,
                                  pi = "-",
                                  tajimas_D = "-")

#load the animal paired-end metrics
animal_paired_end_metrics = read.delim("/Users/sana/Documents/AIM3/006-all_runs_variant_files/0-paired_end/paired_end_df.tsv")
animal_paired_end_metrics = rbind(animal_paired_end_metrics, metric_df_for_no_seg[2:dim(metric_df_for_no_seg)[1], ])
#keep only the ones with cleaned up metadata
animal_paired_end_metrics = animal_paired_end_metrics[which(animal_paired_end_metrics$SRA_ID %in% all_metadata$Run), ]

#load human paired-end metrics
human_paired_end_metrics = read.delim("/Users/sana/Documents/AIM3/006-all_runs_variant_files/human_runs/dataframe.tsv")
human_paired_end_metrics = rbind(human_paired_end_metrics, metric_df_for_no_seg[1, ])
human_paired_end_metrics = human_paired_end_metrics[which(human_paired_end_metrics$SRA_ID %in% all_metadata$Run), ]

#add metadata to the metrics and bind all
d = rbind(animal_paired_end_metrics, human_paired_end_metrics)
c = all_metadata[match(d$SRA_ID, all_metadata$Run), ]
all_data = cbind(d, 
                 c)
write.table(all_data, file = "all_data_w_human.tsv", sep = "\t", row.names = FALSE)

all_data$Species = relevel(as.factor(all_data$Species), ref = "human")
all_data$lab_infected = relevel(as.factor(all_data$lab_infected), ref = "natural")
all_data$tissue = relevel(as.factor(all_data$tissue), ref = "nasal/pharyngeal-swab")
all_data$assay_type = relevel(as.factor(all_data$assay_type), ref = "AMPLICON")

library(car)
library(MASS)
library(lme4)


n_seg = glmer(n_segregating_sites ~  Species + lab_infected + 
                   (1 | tissue) +
                   (1 | average_depth) +
                   #(1 | method_of_collection) +
                   #(1 | runtype) +
                   (1 | assay_type),
                 data = all_data, family = "poisson",
                 control = glmerControl(optimizer = "bobyqa"))
Anova(n_seg)
summary(n_seg)
exp(confint(n_seg))


######## creating a model for pi

pi = glmer(pi ~  Species + lab_infected + 
        (1 | tissue) +
        #(1 | sample_average_depth) +
        (1 | method_of_collection) +
        #(1 | runtype) +
        (1 | assayType), data = tab,
      family = gaussian(link = "log")
)
Anova(pi)
summary(pi)

