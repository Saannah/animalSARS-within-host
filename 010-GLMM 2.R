### load the metrics dataframes

tab_150 = read.delim("/Users/Sana/Documents/AIM3/010-metrics_table/metrics_df_forparsed_150.tsv")
tab_50 = read.delim("/Users/Sana/Documents/AIM3/010-metrics_table/metrics_df_forparsed_50.tsv")


################ only paired-ends with a threshold of 50
tab = tab_50[which(tab_50$runtype == "PAIRED_END"), ]
tab$n_segregating_sites = as.numeric(tab$n_segregating_sites)
tab$pi = as.numeric(tab$pi)
tab$tajimas_D = as.numeric(tab$tajimas_D)
tab$sample_average_depth = as.numeric(tab$sample_average_depth)
tab$vcf_avg_depth_before_filtering_sites = as.numeric(tab$vcf_avg_depth_before_filtering_sites)
tab$vcf_avg_depth_after_filtering_sites = as.numeric(tab$vcf_avg_depth_after_filtering_sites)

### convert categorical columns to factor
tab$species = as.factor(tab$species)
tab$lab_infected = as.factor(tab$lab_infected)
tab$SRA_ID = as.factor(tab$SRA_ID)
tab$tissue = as.factor(tab$tissue)
tab$method_of_collection = as.factor(tab$method_of_collection)
tab$assayType = as.factor(tab$assayType)
tab$SRA_ID = as.factor(tab$SRA_ID)

tab$sample_average_depth = scale(tab$sample_average_depth)
### bin coverag elevels
library(dplyr)
tab %>% mutate(new_bin = cut(sample_average_depth, breaks=3)) -> tab
tab$new_bin = as.factor(tab$new_bin)

library(car)
library(MASS)
library(lme4)
qqp(tab$n_segregating_sites, "norm")

qqp(tab$n_segregating_sites, "lnorm")

poisson <- fitdistr(tab$n_segregating_sites, "Poisson")
qqp(tab$n_segregating_sites, "pois",  lambda = poisson$estimate)

gamma <- fitdistr(tab$n_segregating_sites + 1, "gamma")
qqp(tab$n_segregating_sites + 1, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])


nbinom <- fitdistr(tab$n_segregating_sites + 1, "Negative Binomial")
qqp(tab$n_segregating_sites+1, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])

#only normalize variables not response
# try binning average depth, but not recommended: use it as a fixed effect
# lab infected as random effect
n_seg = glmer.nb(n_segregating_sites ~  species + sample_average_depth +
        #(1 | new_bin) +
        (1 | lab_infected) +
        (1 | tissue) +
        #(1 | SRA_ID) +
        (1 | method_of_collection),
        data = tab,
        control = glmerControl(optimizer = "bobyqa"))
Anova(n_seg)
summary(n_seg)

## look at residuals plots
qqp(tab$pi, "lnorm")
pi = glmer(pi ~  species + sample_average_depth +
            (1 | lab_infected) + 
            (1 | tissue) +
            #(1 | new_bin) +
            #(1|SRA_ID),
            (1 | method_of_collection),
           data = tab,
           family = gaussian(link = "log"),
           control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))
           )
Anova(pi)
summary(pi)



### tajima's D
# 
# qqp(tab$tajimas_D, "norm")
# 
# 

# 
# t = tab$tajimas_D[which(!is.na(tab$tajimas_D))]
# poisson <- fitdistr(t, "Poisson")
# qqp(t, "pois",  lambda = poisson$estimate)
# 
# gamma <- fitdistr(t+1, "gamma")
# qqp(t + 1, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])
# 
# 
# nbinom <- fitdistr(t + 1, "Negative Binomial")
# qqp(t +1, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])

## tajima's D looks normal for paired end reads only
tab$sample_average_depth = scale(tab$sample_average_depth)
tajima = lmer(tajimas_D ~  species +  sample_average_depth +
                 (1 | lab_infected) + 
                 (1 | tissue) +
                 #(1 | new_bin) +
                 (1 | method_of_collection),
                 #(1 | SRA_ID),
               data = tab,
               control = lmerControl(optimizer = "bobyqa"))
Anova(tajima)
summary(tajima)





library(DHARMa)
sims_sp <- simulateResiduals(n_seg)
png('nseg_new2.png')
plot(sims_sp,quantreg = T)
dev.off()
