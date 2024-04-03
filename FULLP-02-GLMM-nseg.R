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
all_data$brdth = all_data$min_50_site_count/29903


############################################################################################################################
################################################# AIC comparisons:



#############################
n_seg_1 = glmer(n_seg ~  tissue_type + Species + lab_infected + 
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
AIC(n_seg_1)  
Anova(n_seg_1)
summary(n_seg_1)
#conf1 = exp(confint(n_seg_1))
#################################################################################################

n_seg_1.nb = glmer.nb(n_seg ~  tissue_type + Species + lab_infected + 
                  (1 | dpth) +
                  #(1 | tissue) +
                  (1 | brdth),
                #(1 | vcf_average_depth) +
                #(1 | method_of_collection) +
                #(1 | runtype) +
                #(1 | lab_infected),
                #(1 | assay_type),
                data = all_data)

allfit = allFit(n_seg_1.nb)
diff_optims_OK <- allfit[sapply(allfit, is, "merMod")]
lapply(allfit, function(x) x@optinfo$conv$lme4$messages)

######
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
overdisp_fun(n_seg_1.nb)

######## running the best model (1) using a quasipoisson distribution:
quasi_table <- function(model,ctab=coef(summary(model))) {
  phi <- sum(residuals(model, type="pearson")^2)/df.residual(model)
  qctab <- within(as.data.frame(ctab),
                  {   `Std. Error` <- `Std. Error`*sqrt(phi)
                  `z value` <- Estimate/`Std. Error`
                  `Pr(>|z|)` <- 2*pnorm(abs(`z value`), lower.tail=FALSE)
                  })
  return(qctab)
}

printCoefmat(quasi_table(n_seg_1),digits=4)



##### testing the best model (model 1) without lymph node tissue samples:
sub_data = all_data[-which(all_data$tissue_type == "lymph"), ]
dim(sub_data)
n_seg_12 = glmer(n_seg ~ Species + lab_infected + 
                  (1 | dpth) +
                  #(1 | tissue) +
                  (1 | brdth), 
                #(1 | vcf_average_depth) +
                #(1 | method_of_collection) +
                #(1 | runtype) +
                #(1 | lab_infected) +
                #(1 | assay_type),
                data = sub_data, family = "poisson",
                control=glmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e5)))
AIC(n_seg_12)  
Anova(n_seg_12)
summary(n_seg_12)


#############################
n_seg_11 = glmer(n_seg ~  tissue_type + Species + lab_infected + 
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
AIC(n_seg_11)  
Anova(n_seg_11)
summary(n_seg_11)
#conf1 = exp(confint(n_seg_1))
#################################################################################################



n_seg_2 = glmer(n_seg ~  Species + lab_infected + 
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

AIC(n_seg_2)  
Anova(n_seg_2)
summary(n_seg_2)

##########################################################################
n_seg_3 = glmer(n_seg ~  tissue_type + Species + #lab_infected + 
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
AIC(n_seg_3)  
Anova(n_seg_3)
summary(n_seg_3)
#conf1 = exp(confint(n_seg_1))



##########################################################################
n_seg_4 = glmer(n_seg ~  Species + dpth + brdth +
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
AIC(n_seg_4)  
Anova(n_seg_4)
summary(n_seg_4)

##########################################################################
n_seg_5 = glmer(n_seg ~  Species + dpth + brdth + lab_infected +
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
AIC(n_seg_5)  
Anova(n_seg_5)
summary(n_seg_5)

##########################################################################
n_seg_6 = glmer(n_seg ~  Species + brdth + lab_infected + tissue_type +
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
AIC(n_seg_6)  
Anova(n_seg_6)
summary(n_seg_6)


n_seg_61 = glmer(n_seg ~  Species + brdth + lab_infected + tissue_type +
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
AIC(n_seg_61)  
Anova(n_seg_61)
summary(n_seg_61)
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
