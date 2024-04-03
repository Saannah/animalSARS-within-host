rm(list = ls())
library(car)
library(MASS)
library(lme4)





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

############################################################################################################################
################################################# AIC comparisons:



#############################
richness_1 = glmer(filtered_richness ~  tissue_type + Species + lab_infected + 
                  (1 | dpth) +
                  #(1 | tissue) +
                  (1 | brdth), 
                #(1 | vcf_average_depth) +
                #(1 | method_of_collection) +
                #(1 | runtype) +
                #(1 | lab_infected) +
                #(1 | assay_type),
                data = all_data, family = "poisson",
                control=glmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e5)))
AIC(richness_1)  
Anova(richness_1)
summary(richness_1)
#################################################################################################




#############################
richness_11 = glmer(filtered_richness ~  tissue_type + Species + lab_infected + 
                   (1 | dpth)+
                   #(1 | tissue) +
                   (1 | brdth)+
                   #(1 | vcf_average_depth) +
                   #(1 | method_of_collection) +
                   #(1 | runtype) +
                   #(1 | lab_infected) +
                   (1 | assay_type),
                 data = all_data, family = "poisson",
                 control=glmerControl(optimizer="bobyqa",
                                      optCtrl=list(maxfun=2e5)))
AIC(richness_11)  
Anova(richness_11)
summary(richness_11)
#conf1 = exp(confint(n_seg_1))
#################################################################################################



richness_2 = glmer(filtered_richness ~  Species + lab_infected + 
                  (1 | dpth) +
                  (1 | tissue_type) +
                  (1 | brdth),
                #(1 | vcf_average_depth) +
                #(1 | method_of_collection) +
                #(1 | runtype) +
                #(1 | lab_infected) +
                #(1 | assay_type),
                data = all_data, family = "poisson",
                control=glmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e5)))

AIC(richness_2)  
Anova(richness_2)
summary(richness_2)

##########################################################################
richness_3 = glmer(filtered_richness ~  tissue_type + Species + #lab_infected + 
                  (1 | dpth) +
                  #(1 | tissue) +
                  (1 | brdth) +
                  #(1 | vcf_average_depth) +
                  #(1 | method_of_collection) +
                  #(1 | runtype) +
                  (1 | lab_infected),
                #(1 | assay_type),
                data = all_data, family = "poisson",
                control=glmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e5)))
AIC(richness_3)  
Anova(richness_3)
summary(richness_3)
#conf1 = exp(confint(n_seg_1))



##########################################################################
richness_4 = glmer(filtered_richness ~  Species + dpth + brdth +
                  #(1 | dpth) +
                  (1 | tissue) +
                  #(1 | brdth) +
                  #(1 | vcf_average_depth) +
                  #(1 | method_of_collection) +
                  #(1 | runtype) +
                  (1 | lab_infected),
                #(1 | assay_type),
                data = all_data, family = "poisson",
                control=glmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e5)))
AIC(richness_4)  
Anova(richness_4)
summary(richness_4)

##########################################################################
richness_5 = glmer(filtered_richness ~  Species + dpth + brdth + lab_infected +
                  #(1 | dpth) +
                  (1 | tissue),
                #(1 | brdth) +
                #(1 | vcf_average_depth) +
                #(1 | method_of_collection) +
                #(1 | runtype) +
                #(1 | lab_infected),
                #(1 | assay_type),
                data = all_data, family = "poisson",
                control=glmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e5)))
AIC(richness_5)  
Anova(richness_5)
summary(richness_5)

##########################################################################
richness_51 = glmer(filtered_richness ~  Species + dpth + brdth + lab_infected +
                     #(1 | dpth) +
                     (1 | tissue)+
                   #(1 | brdth) +
                   #(1 | vcf_average_depth) +
                   #(1 | method_of_collection) +
                   #(1 | runtype) +
                   #(1 | lab_infected),
                   (1 | assay_type),
                   data = all_data, family = "poisson",
                   control=glmerControl(optimizer="bobyqa",
                                        optCtrl=list(maxfun=2e5)))
AIC(richness_51)  
Anova(richness_51)
summary(richness_51)





##########################################################################
richness_6 = glmer(filtered_richness ~  Species + brdth + lab_infected + tissue_type +
                  (1 | dpth),
                #(1 | tissue),
                #(1 | brdth) +
                #(1 | vcf_average_depth) +
                #(1 | method_of_collection) +
                #(1 | runtype) +
                #(1 | lab_infected),
                #(1 | assay_type),
                data = all_data, family = "poisson",
                control=glmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e5)))
AIC(richness_6)  
Anova(richness_6)
summary(richness_6)


richness_61 = glmer(filtered_richness ~  Species + brdth + lab_infected + tissue_type +
                   (1 | dpth)+
                   #(1 | tissue),
                   #(1 | brdth) +
                   #(1 | vcf_average_depth) +
                   #(1 | method_of_collection) +
                   #(1 | runtype) +
                   #(1 | lab_infected),
                   (1 | assay_type),
                 data = all_data, family = "poisson",
                 control=glmerControl(optimizer="bobyqa",
                                      optCtrl=list(maxfun=2e5)))
AIC(richness_61)  
Anova(richness_61)
summary(richness_61)
 # 
# old_data = read.delim("/users/sana/Documents/AIM3/all_within_host_data.tsv")
# new_data = read.delim("/users/sana/Documents/AIM3/021-full_pipeline/1-metrics/metrics.tsv")
# new_data = new_data[which(new_data$Run %in% old_data$Run), ]
# 
# new_data = new_data[order(new_data$Run), ]
# old_data = old_data[order(old_data$Run), ]
# 
# data.frame(old = old_data$n_segregating_sites,
#            new = new_data$n_seg)
# vcf = read.delim("/users/sana/Documents/AIM3/021-full_pipeline/00-ALL_FILES/ERR5885053.freebayes_calling.consensus.vcf_table.tsv.unfiltered.parsed")
# ggplot(data = new_data) + geom_point(aes(x = Species, y = n_seg))
# 

#### trying the best model 6-1 without the lymph node data points
sub_data = all_data[-which(all_data$tissue_type == "lymph"), ]
dim(sub_data)
richness_612 = glmer(filtered_richness ~  Species + brdth + lab_infected + 
                      (1 | dpth)+
                      #(1 | tissue),
                      #(1 | brdth) +
                      #(1 | vcf_average_depth) +
                      #(1 | method_of_collection) +
                      #(1 | runtype) +
                      #(1 | lab_infected),
                      (1 | assay_type),
                    data = sub_data, family = "poisson",
                    control=glmerControl(optimizer="bobyqa",
                                         optCtrl=list(maxfun=2e5)))
AIC(richness_612)  
Anova(richness_612)
summary(richness_612)

