### plot tajimasd for animals
rm(list = ls())
tab = read.delim("/users/sana/Documents/AIM3/021-full_pipeline/1-metrics/metrics_strict_thresholds_min50_5percent.tsv")
tajimasd_plot = ggplot(tab, aes(x=factor(Species, level = c("human", "cat", "dog", "mink", "deer")), 
                                    y=tajimasd)) + 
  geom_boxplot(aes(x=Species, y=tajimasd), fill = "#4A5899", alpha = 0.1, outlier.shape = NA) +
  geom_jitter(aes(color = lab_infected), width = 0.35, height=0, size = 1) + 
  scale_color_manual(values = c("#4A5899", "#CD5D67")) + 
  #ylim(0,3.5) +
  theme_bw() + xlab("Species") + ggtitle("Tajimas D - not scaled") +
  ylab("Tajimas D") + theme(legend.position = "none")

## scale tajimasd
tab$tajimasd_scaled = as.array(scale(tab$tajimasd))
tajimasd_plot_scaled = ggplot(tab, aes(x=factor(Species, level = c("human", "cat", "dog", "mink", "deer")), 
                                y=tajimasd_scaled)) + 
  geom_boxplot(aes(x=Species, y=tajimasd_scaled), fill = "#4A5899", alpha = 0.1, outlier.shape = NA) +
  geom_jitter(aes(color = lab_infected), width = 0.35, height=0, size = 1) + 
  scale_color_manual(values = c("#4A5899", "#CD5D67")) + 
  #ylim(0,3.5) +
  theme_bw() + xlab("Species") + ggtitle("Tajimas D - scaled") +
  ylab("Tajimas D") + theme(legend.position = "none")


##### glmm for tahimasd
library(car)
library(MASS)
library(lme4)

tab$tajimasd_squared = tab$tajimasd ^ 2

qqp(tab$tajimasd_squared, "norm")

qqp(tab$tajimasd_squared, "lnorm")



gamma <- fitdistr(tab$tajimasd, "gamma")
qqp(tab$n_segregating_sites + 1, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])


nbinom <- fitdistr(tab$tajimasd, "Negative Binomial")
qqp(tab$n_segregating_sites+1, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])




#### wilcoxon tests for tajimas D

table = data.frame(species1 = vector(),
                   species2 = vector(),
                   pvalue = vector())

for(species in c("cat", "dog", "mink", "deer")){
  
  
  t = wilcox.test(tab$tajimasd_scaled[which(tab$Species == "human")],
                  tab$tajimasd_scaled[which(tab$Species == species)])
  table = rbind(table,
                data.frame(species1 = "human",
                           species2 = species,
                           pvalue = t$p.value))
}

write.table(table,
            file = "/users/sana/Documents/AIM3/021-full_pipeline/8-tajimasd/wilcoxon.tsv",
            sep = "\t",
            row.names = FALSE)




