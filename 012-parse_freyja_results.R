library(stringr)

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


lab_infected = tab[which(tab$lab_infected == "lab"), ]

freyja_df = data.frame(SRA_id = vector(),
                       lineages = vector(),
                       abundance = vector())
for(i in 1:length(file_names)){
  t = read.delim(file_names[i], sep = "\t")
  freyja_df[nrow(freyja_df)+1, ] = c(n[i], t[2,2], t[3,2])
}

## include only SRA's present in table
freyja_df = freyja_df[which(freyja_df$SRA_id %in% tab$SRA_ID), ]

lab_infected_lineages = freyja_df[match(lab_infected$SRA_ID, freyja_df$SRA_id),]

lineage_df = data.frame(SRA_id = character(),
                        richness = numeric(),
                        shannon_diversity = numeric())

### function to calculate shannon diversity
calculate_shannon = function(vec){
  return(-1 * sum(vec * log(vec)))

}

for(i in 1:dim(freyja_df)[1]){
  vec = as.numeric(unlist(str_split(freyja_df$abundance[i], " ")))
  if(length(which(vec <= abundance_threshold)) > 0){
    vec = vec[-which(vec <= abundance_threshold)]
  }
  lineage_df[nrow(lineage_df) + 1, ] = c(freyja_df$SRA_id[i], length(vec), calculate_shannon(vec))
}


### add column to tab
tab$richness = -1
tab$shannon_diversity = -1
tab$richness[match(lineage_df$SRA_id, tab$SRA_ID)] =lineage_df$richness
tab$shannon_diversity[match(lineage_df$SRA_id, tab$SRA_ID)] =lineage_df$shannon_diversity

### rmove first row of tab because we don't have freyja output for this
### likely because a high quality alignment could not be produced
tab = tab[-1, ]

## write table
write.table(tab, file = "richness_and_shannon.tsv", sep = "\t", row.names = FALSE)
