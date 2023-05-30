### load the metrics dataframes

tab_150 = read.delim("/Users/admin/Documents/AIM3/010-metrics_table/metrics_df_forparsed_150.tsv")
tab_50 = read.delim("/Users/admin/Documents/AIM3/010-metrics_table/metrics_df_forparsed_50.tsv")


################ only paired-ends with a threshold of 50
tab = tab_50[which(tab_50$runtype == "PAIRED_END"), ]
tab$n_segregating_sites = as.numeric(tab$n_segregating_sites)
tab$pi = as.numeric(tab$pi)
tab$tajimas_D = as.numeric(tab$tajimas_D)
tab$sample_average_depth = as.numeric(tab$sample_average_depth)
tab$vcf_avg_depth_before_filtering_sites = as.numeric(tab$vcf_avg_depth_before_filtering_sites)
tab$vcf_avg_depth_after_filtering_sites = as.numeric(tab$vcf_avg_depth_after_filtering_sites)


####### normalize columns of interest to mean zero and sd 1
tab$n_segregating_sites = scale(tab$n_segregating_sites)
tab$pi = scale(tab$pi)
tab$tajimas_D = scale(tab$tajimas_D)
tab$sample_average_depth = scale(tab$sample_average_depth)


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

n_seg = glmer.nb(n_segregating_sites ~  species + lab_infected + 
        (1 | tissue) +
        (1 | sample_average_depth) +
        (1 | method_of_collection) +
        #(1 | runtype) +
        (1 | assayType),
        data = tab,
        control = glmerControl(optimizer = "bobyqa"))
Anova(n_seg)
summary(n_seg)

n_seg = glmer(n_segregating_sites ~  species + lab_infected + 
                   (1 | tissue) +
                   (1 | sample_average_depth) +
                   (1 | method_of_collection) +
                   #(1 | runtype) +
                   (1 | assayType),
                 data = tab,
                 control = glmerControl(optimizer = "bobyqa"),
              family = gaussian(link = "log"))






qqp(tab$pi, "lnorm")
qqp(tab$pi, "norm")
pi = glmer(pi ~  species + lab_infected + 
            (1 | tissue) +
            (1 | sample_average_depth) +
            (1 | method_of_collection) +
            #(1 | runtype) +
            (1 | assayType), data = tab,
           family = gaussian(link = "log")
           )
Anova(pi)
summary(pi)



### tajima's D
# 
# qqp(tab$tajimas_D, "norm")
# 
# qqp(tab$tajimas_D, "lnorm")
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
tajima = lmer(tajimas_D ~  species + lab_infected + 
            (1 | sample_average_depth ) +
            (1 | method_of_collection) +
            (1 | tissue) +
            #(1 | runtype) +
            (1 | assayType), 
            data = tab)
Anova(tajima)

# 
# 
# library(lme4)
# m = glmer(n_segregating_sites ~ sample_average_depth + species + lab_infected + 
#             (1 | method_of_collection) +
#             (1 | tissue) +
#             (1 | assayType), data = tab)
# Anova(m)
# 
# m = glmer(pi ~ sample_average_depth + species + lab_infected + 
#             (1 | method_of_collection) +
#             (1 | tissue) +
#             (1 | assayType), data = tab)
# Anova(m)
# 
# 
# m = glmer(tajimas_D ~ sample_average_depth + species + lab_infected + 
#             (1 | method_of_collection) +
#             (1 | tissue) +
#             (1 | assayType), data = tab)
# Anova(m)
