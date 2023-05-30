library(stringr)
library(ggplot2)

### tunable variables
abundance_threshold = 0.025

setwd("/Users/sana/Documents/AIM3/011-freyja")
file_names = read.delim("names_list.tsv")
file_names = file_names[, 1]
n = vector()
names = str_split(file_names, "\\.")
for(i in 1:length(names)){
  n = c(n, unlist(names[i])[1])
}

## load table for paired ends
tab_50 = read.delim("/Users/Sana/Documents/AIM3/010-metrics_table/metrics_df_forparsed_50.tsv")
tab = tab_50[which(tab_50$runtype == "PAIRED_END"), ]


#lab_infected = tab[which(tab$lab_infected == "lab"), ]

freyja_df = data.frame(SRA_id = vector(),
                       lineages = vector(),
                       abundance = vector())
for(i in 1:length(file_names)){
  t = read.delim(file_names[i], sep = "\t")
  freyja_df[nrow(freyja_df)+1, ] = c(n[i], t[2,2], t[3,2])
}

## include only SRA's present in table
freyja_df = freyja_df[which(freyja_df$SRA_id %in% tab$SRA_ID), ]
#lab_infected_lineages = freyja_df[match(lab_infected$SRA_ID, freyja_df$SRA_id),]

#### function to find all lineages present and remove low abundances from df
all_lin = vector()
freyja_df$richness = 0
for(i in 1:dim(freyja_df)[1]){
  vec_abund = as.numeric(unlist(str_split(freyja_df$abundance[i], " ")))
  vec_lin = unlist(str_split(freyja_df$lineages[i], " "))
  
  if(length(which(vec_abund <= abundance_threshold)) > 0){
    rm = which(vec_abund <= abundance_threshold)
    vec_abund = vec_abund[-rm]
    vec_lin = vec_lin[-rm]
    
    freyja_df$lineages[i] = paste(vec_lin, collapse = ' ')
    freyja_df$abundance[i] = paste(vec_abund, collapse = ' ')
  }
  freyja_df$richness[i] = length(vec_lin)
  all_lin = unique(c(all_lin, vec_lin))
}

## check if there are 0 richness samples with "other" lineages
freyja_df = freyja_df[order(freyja_df$richness),]
head(freyja_df)
freyja_df[1,] = c("ERR8314858", "other",  0.984811, 1)
freyja_df[2,] = c("SRR18231914", "other",  0.984811, 1)
head(freyja_df)
all_lin = c(all_lin, "other")


### sort freyja_df based on richness (lowest to highet)
freyja_df = freyja_df[order(freyja_df$richness),]


## create dataframe for stackplot
sample_array = vector()
lineage_array = vector()
abundance_array = vector()

for(i in 1:dim(freyja_df)[1]){
  ids = rep(freyja_df$SRA_id[i], length(all_lin))
  abundances = rep(0, length(all_lin))
  abundances[match(unlist(str_split(freyja_df$lineages[i], " ")), all_lin)] = 
    as.numeric(unlist(str_split(freyja_df$abundance[i], " ")))
  sample_array = c(sample_array, ids)
  lineage_array = c(lineage_array, all_lin)
  abundance_array = c(abundance_array, abundances)
}

stackbar_df = data.frame(sample = sample_array,
                         lineage = lineage_array,
                         abundance = abundance_array)
stackbar_df = stackbar_df[-which(stackbar_df$abundance == 0),]

stackbar_df$abundance = as.numeric(stackbar_df$abundance)

species = "mink"
p = ggplot(stackbar_df[which(stackbar_df$sample %in% tab$SRA_ID[which(tab$species == species)]), ], 
       aes(fill=lineage, y=abundance, x=sample)) + 
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle(species)
setwd("/Users/sana/Documents/AIM3/012-plot_richness_and_shannon/")
ggsave(p, file = "mink_abundances.pdf", device = "pdf", width = 30, height = 5)
