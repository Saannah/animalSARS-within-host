### create exploratory plots 
library(ggplot2)
#load all_data from script 018

all_data$colors = all_data$lab_infected
all_data$colors[which(all_data$Species == "human")] = "human"
#plot number of segregating sites
n_seg_plot = ggplot(all_data, aes(x=factor(Species, level = c("human", "cat", "dog", "mink", "deer")), 
                         y=n_segregating_sites, fill = colors)) + 
  geom_jitter(aes(color = colors), width = 0.05, height=0, size = 1) + 
  scale_color_manual(values = c("darkred", "#DDA7AB", "#4A5899"), labels = c("Human", "Lab Infection", "Natural Infection")) + ylim(0,100) +
  theme_bw() + xlab("Species") +
  ylab("Number of Segregating Sites") + theme(legend.position = 'none')

#remove dashed pi's
pi_data = all_data[-which(all_data$pi == "-"), ]
pi_data$pi = as.numeric(pi_data$pi)
pi_plot = ggplot(pi_data, aes(x=factor(Species, level = c("human", "cat", "dog", "mink", "deer")), 
                        y=pi, fill = colors)) + 
  geom_jitter(aes(color = colors), width = 0.05, height=0, size = 1) + 
  scale_color_manual(values = c("darkred", "#DDA7AB", "#4A5899")) + ylim(0,20) +
  theme_bw() + xlab("Species") +
  ylab("Within-host pairwise differences") + theme(legend.position = 'none')
