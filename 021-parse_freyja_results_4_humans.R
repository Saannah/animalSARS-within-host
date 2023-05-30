library(stringr)

### tunable variables
abundance_threshold = 0.025

setwd("/Users/sana/Documents/AIM3/017-freyja_results_4_human/")
file_names = read.delim("names_list.tsv", header = FALSE)
file_names = file_names[, 1]
n = vector()
names = str_split(file_names, "\\.")
for(i in 1:length(names)){
  n = c(n, unlist(names[i])[1])
}

## load all_data from script 18
n = n[which(n %in% all_data$SRA_ID)]


lab_infected = all_data[which(all_data$lab_infected == "lab"), ]

freyja_df = data.frame(SRA_id = vector(),
                       lineages = vector(),
                       abundance = vector())

for(i in 1:(length(file_names)-1)){
  t = read.delim(file_names[i], sep = "\t")
  freyja_df[nrow(freyja_df)+1, ] = c(n[i], t[2,2], t[3,2])
}

## include only SRA's present in table
freyja_df = freyja_df[which(freyja_df$SRA_id %in% all_data$SRA_ID), ]

#lab_infected_lineages = freyja_df[match(lab_infected$SRA_ID, freyja_df$SRA_id),]

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

human_lineages = lineage_df

### run script 12 to get lineages_df for animals before running the rest of the code
animal_lineages = lineage_df

all_freyja_output = rbind(human_lineages, animal_lineages)
all_freyja_output$Species = "tbd"

#check if all freyja outputs exist in the all_data dataframe:
length(which(all_freyja_output$SRA_id %in% all_data$SRA_ID))
all_data = all_data[which(all_data$SRA_ID %in% all_freyja_output$SRA_id), ]
all_freyja_output$Species[match(all_data$SRA_ID, all_freyja_output$SRA_id)] = all_data$Species
all_freyja_output$colors[match(all_data$SRA_ID, all_freyja_output$SRA_id)] = all_data$colors
all_freyja_output$richness = as.numeric(all_freyja_output$richness)
