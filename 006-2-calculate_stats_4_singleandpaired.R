##### customize variables:
min_coverage_thresh = 50
fixation_freq = 0.95

### script for creating metrics for all single and paired-end reads
library(stringr)
#### 1. load depth summary files (from genpipes)
names = read.delim("/Users/sana/Documents/AIM3/007-all_runs_depth_summary_files/depth_file_list.tsv", header = FALSE)
names = names[,1]

#### create dataframe from coverage files
coverage_df = data.frame(id = vector(), mean_coverage = vector())

setwd("/Users/sana/Documents/AIM3/007-all_runs_depth_summary_files")
for(i in names){
  t = read.delim(i)
  coverage_df[nrow(coverage_df) + 1, ] = c(t$sample_id[1], t$mean[1])
}
coverage_df$mean_coverage = as.numeric(coverage_df$mean_coverage) #convert from string to numeric

#### 2. remove average depth coverage outliers annd average < 50
coverage_df = coverage_df[-which(
  coverage_df$mean_coverage %in% boxplot.stats(coverage_df$mean_coverage)$out), ]
coverage_df = coverage_df[-which(coverage_df$mean_coverage < min_coverage_thresh), ]

#### 3. load parsed tables names list (files that have .parsed extention)
pooled_parsed_names = read.delim("/Users/sana/Documents/AIM3/006-all_runs_variant_files/1-pool/parsed_names.tsv", header = FALSE)
pooled_parsed_names = pooled_parsed_names[, 1]

file_name_id = vector()
pooled_parsed_names_obj = str_split(pooled_parsed_names, "\\.")
for(i in 1:length(pooled_parsed_names_obj)){
  file_name_id = c(file_name_id, unlist(pooled_parsed_names_obj[i])[1])
}

#remove file names that are depth outliters (do not exist in coverage_df)
file_names_clean = pooled_parsed_names[which(file_name_id %in% coverage_df$id)]

#create dataframe for metrics
df = data.frame(SRA_ID = vector(),
                sample_average_depth = vector(),
                vcf_avg_depth_before_filtering_sites = vector(),
                vcf_avg_depth_after_filtering_sites = vector(),
                n_segregating_sites = vector(),
                pi = vector(),
                tajimas_D = vector())

#### set working directory to where the variant files are
setwd("/Users/sana/Documents/AIM3/006-all_runs_variant_files/1-pool/")

for (n in file_names_clean){
  print(n)
  ### load the table
  tab = read.delim(n)
  avg_dpth = mean(tab$DEPTH)
  
  #flag for average depth
  flag = "y"
  if(avg_dpth < min_coverage_thresh){flag = "n"}
  id = unlist(str_split(n, "\\."))[1]
  
  #remove sites with low coverage
  if(length(which(tab$DEPTH < min_coverage_thresh)) > 0){tab = tab[-which(tab$DEPTH < min_coverage_thresh),]}
  unfiltered_avg_dpth = avg_dpth
  avg_dpth = mean(tab$DEPTH)
  
  if(length(which(tab$VAF > fixation_freq)) > 0){tab = tab[-which(tab$VAF > fixation_freq),]}
  
  
  if(dim(tab)[1] == 0){
    df[(nrow(df) + 1), ] = c(id, 
                             coverage_df$mean_coverage[which(coverage_df$id == id)],
                             unfiltered_avg_dpth,
                             avg_dpth,
                             0, "-", "-")
    next}
  #print("here")
  n_seg_sites = dim(tab)[1]
  
  tab$pi_count = tab$ALT_COUNT * (tab$DEPTH - tab$ALT_COUNT)
  tab$choose2 = tab$DEPTH * (tab$DEPTH - 1) * 0.5
  pi = sum(tab$pi_count) / sum(tab$choose2) * n_seg_sites
  n = mean(tab$DEPTH)
  x = 1:n
  a1 = sum(1 / x)
  
  a2 = sum((1/x)*(1/x))
  
  b1 = (n + 1) / (3 * (n - 1))
  
  b2 = (2 * (n**2 + n + 3)) / (9 * n * (n - 1))
  
  c1 = b1 - (1 / a1)
  
  c2 = b2 - ((n + 2) / (a1 * n)) + (a2 / a1**2)  
  
  e1 = c1 / a1
  
  e2 = c2 / (a1**2 + a2)
  
  expected_sd = sqrt(e1 * n_seg_sites + e2 * n_seg_sites * (n_seg_sites - 1))
  
  wattersons_theta = n_seg_sites / a1
  
  if (expected_sd == 0){tajimas_d = float("NaN")}
  else
    tajimas_d = (pi - wattersons_theta) / expected_sd
  
  df[(nrow(df) + 1), ] = c(id, 
                           coverage_df$mean_coverage[which(coverage_df$id == id)],
                           unfiltered_avg_dpth, 
                           avg_dpth, 
                           n_seg_sites, 
                           pi, 
                           tajimas_d)
}

metrics = df
### load the final version of metadata with lab and tissue information
cat = read.delim("/Users/sana/Documents/AIM3/009-tissue&method_categorization/cat_tissues.tsv")
dog = read.delim("/Users/sana/Documents/AIM3/009-tissue&method_categorization/dog_tissues.tsv")
deer = read.delim("/Users/sana/Documents/AIM3/009-tissue&method_categorization/deer_tissues.tsv")
mink = read.delim("/Users/sana/Documents/AIM3/009-tissue&method_categorization/mink_tissues.tsv")
all_metadata = rbind(cat[, 1:9],
                     dog[, 1:9],
                     mink[, 1:9],
                     deer[, 1:9])
#remove data that do not exist in metadata and vice versa
metrics = metrics[which(metrics$SRA_ID %in% all_metadata$Run), ]
all_metadata = all_metadata[which(all_metadata$Run %in% metrics$SRA_ID), ]

#add a column to metrics dataframe and specify species
metrics$species = "TBD"
metrics$species[which(metrics$SRA_ID %in% cat$Run)] = "cat"
metrics$species[which(metrics$SRA_ID %in% dog$Run)] = "dog"
metrics$species[which(metrics$SRA_ID %in% mink$Run)] = "mink"
metrics$species[which(metrics$SRA_ID %in% deer$Run)] = "deer"

#add a columnn to metrics dataframe and specify lab infected or not
metrics$lab_infected = "TBD"
metrics$lab_infected[match(all_metadata$Run, metrics$SRA_ID)] = all_metadata$lab_infected

#add a column to specify runtype
readset = read.delim("/Users/sana/Documents/AIM3/006-all_runs_variant_files/readset_for_paired_and_single.tsv")
readset = readset[which(readset$Sample %in% metrics$SRA_ID), ]
metrics$runtype = "TBD"
metrics$runtype[match(readset$Sample, metrics$SRA_ID)] = readset$RunType

#add a column for tissue type and method
metrics$tissue = "TBD"
metrics$tissue[match(all_metadata$Run, metrics$SRA_ID)] = all_metadata$tissue
metrics$method_of_collection = "TBD"
metrics$method_of_collection[match(all_metadata$Run, metrics$SRA_ID)] = all_metadata$method

#add a column for assay type
metrics$assayType = "TBD"
metrics$assayType[match(all_metadata$Run, metrics$SRA_ID)] = all_metadata$Assay.Type

setwd("/Users/sana/Documents/AIM3/010-metrics_table/")
write.table(metrics, file = "metrics_df_forparsed_50.tsv", sep = "\t", row.names = FALSE)



