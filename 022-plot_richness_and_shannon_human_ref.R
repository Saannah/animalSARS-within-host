#run script 21 to get all_freyja_output

richness_plot = ggplot(all_freyja_output, aes(x=factor(Species, level = c("human", "cat", "dog", "mink", "deer")), 
                         y=richness, fill = colors)) + 
  geom_jitter(aes(color = colors), width = 0.05, height=0, size = 1) + 
  scale_color_manual(values = c("darkred", "#DDA7AB", "#4A5899"), labels = c("Human", "Lab Infection", "Natural Infection")) + ylim(0,25) +
  theme_bw() + xlab("Species") +
  ylab("Richness (# of present VOCs)") + theme(legend.position = 'none')

ggarrange(n_seg_plot, pi_plot, richness_plot, nrow = 1)
