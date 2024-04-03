rm(list = ls())
setwd("/Users/sana/Documents/AIM3/021-full_pipeline/00-ALL_FILES")
library(stringr)

minimum_depth = 50
fixation_freq = 0.05
min_alt_dp = 50

## load vcf file names
vcfnames = read.delim("filenames", header = FALSE)
vcfnames = vcfnames[, 1]

#create vectors for coveage and freyja file names
samplenames = vector()
coveragenames = vector()
freyjanames = vector()
for(i in 1:length(vcfnames)){
  name = str_split(vcfnames[i], "\\.")
  samplenames = c(samplenames, unlist(name)[1])
  coveragenames = c(coveragenames, 
                    paste0(unlist(name)[1], ".sorted.filtered.primerTrim.bam.bed.coverage"))
  freyjanames = c(freyjanames, 
                  paste0(unlist(name)[1], ".sorted.filtered.primerTrim.bam_VOC.tsv"))
}

## function that filters the VCF file based on the thresholds
is_acceptable_snp = function(r, VCF, minimum_depth, fixation_freq, min_alt_dp){
  row = VCF[r, ]
  flag = 0
  #excludes indels from the file
  if(nchar(row$REF) > 1){
    return(0)
  }
  #excludes low depth calls
  if(row$DEPTH < minimum_depth){
    return(0)
  }
  #include if alt freq is within range
  if(row$VAF >= fixation_freq){
    if(row$VAF <= (1 - fixation_freq)){
      return(1)
    }
  }
  #include if alt depth is in range
  # if(row$ALT_COUNT >= min_alt_dp){
  #   return(1)
  # }
  #return 0 if non apply
  return(0)
}


metrics = data.frame(sample_id = vector(),
                     n_seg = vector(),
                     pi = vector(),
                     tajimasd = vector(),
                     min_50_site_count = vector(),
                     mean_depth_from_bed = vector(),
                     vcf_average_depth = vector(),
                     notes = vector())

for(i in 1:length(vcfnames)){
  print(samplenames[i])
  
  covfile = read.delim(coveragenames[i], header = FALSE)
  
  #if the vcf file is empty, then there were no snps called in the first place
  if(file.size(vcfnames[i]) == 0){
    metrics[nrow(metrics) + 1, ] = c(sample_id = samplenames[i],
                                     n_seg = 0,
                                     pi = 0,
                                     tajimasd = 0,
                                     min_50_site_count = length(which(covfile$V4 >= minimum_depth)),
                                     mean_depth_from_bed = mean(covfile$V4),
                                     vcf_average_depth = 0,
                                     notes = "no snps called")
    next
  }
  
  VCF = read.delim(vcfnames[i])
  tab = VCF[which(unlist(lapply(1:nrow(VCF), is_acceptable_snp, VCF, minimum_depth, fixation_freq, min_alt_dp)) == 1), ]
  
  if(nrow(tab) == 0){
    metrics[nrow(metrics) + 1, ] = c(sample_id = samplenames[i],
                                     n_seg = 0,
                                     pi = 0,
                                     tajimasd = 0,
                                     min_50_site_count = length(which(covfile$V4 >= minimum_depth)),
                                     mean_depth_from_bed = mean(covfile$V4),
                                     vcf_average_depth = 0,
                                     notes = "no eligible snps")
    next
  }
  
  
  n_seg = dim(tab)[1]
  tab$pi_count = tab$ALT_COUNT * (tab$DEPTH - tab$ALT_COUNT)
  tab$choose2 = tab$DEPTH * (tab$DEPTH - 1) * 0.5
  pi = (sum(tab$pi_count) / sum(tab$choose2)) * n_seg
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
  
  expected_sd = sqrt(e1 * n_seg + e2 * n_seg * (n_seg - 1))
  
  wattersons_theta = n_seg / a1
  
  if(expected_sd == 0){
    tajimas_d = float("NaN")
    metrics[nrow(metrics) + 1, ] = c(sample_id = samplenames[i],
                                     n_seg = n_seg,
                                     pi = pi,
                                     tajimasd = tajimas_d,
                                     min_50_site_count = length(which(covfile$V4 >= minimum_depth)),
                                     mean_depth_from_bed = mean(covfile$V4),
                                     vcf_average_depth = mean(tab$DEPTH),
                                     notes = "expected sd for tajima's d is 0")
  }
  if(expected_sd != 0){
    tajimas_d = (pi - wattersons_theta) / expected_sd
    metrics[nrow(metrics) + 1, ] = c(sample_id = samplenames[i],
                                     n_seg = n_seg,
                                     pi = pi,
                                     tajimasd = tajimas_d,
                                     min_50_site_count = length(which(covfile$V4 >= minimum_depth)),
                                     mean_depth_from_bed = mean(covfile$V4),
                                     vcf_average_depth = mean(tab$DEPTH),
                                     notes = "all good")
  }
  
}




#### add shannon and richness to the data

### the threshold for filtering abundances comes from the lab infected samples (script 012)
abundance_threshold = 0.025



#create a df for freyja information to be filled in
freyja_df = data.frame(SRA_id = vector(),
                       lineages = vector(),
                       abundances = vector(),
                       summary = vector())

for(i in 1:length(freyjanames)){
  freyja_file = read.delim(freyjanames[i])
  freyja_df[nrow(freyja_df) + 1, ] = c(samplenames[i],
                                       freyja_file[2, 2],
                                       freyja_file[3, 2],
                                       freyja_file[1, 2])
}


calculate_shannon_for_1_sample = function(rownum){
  sample = freyja_df[rownum, ]
  abundance_vector = as.numeric(unlist(str_split(sample$abundances, pattern = " ")))
  lineage_names = unlist(str_split(sample$lineages, pattern = " "))
  
  unfiltered_shannon = -sum(log(abundance_vector) * abundance_vector)
  unfiltered_richness = length(abundance_vector)
  
  if(length(which(abundance_vector < abundance_threshold)) > 0){
    ind = which(abundance_vector < abundance_threshold)
    abundance_vector = abundance_vector[-ind]
    lineage_names = lineage_names[-ind]
  }
  
  if(length(abundance_vector) == 0){
    filtered_shannon = "lineage split"
    filtered_richness = "lineage split"
    return(data.frame(unfiltered_shannon = unfiltered_shannon,
                      filtered_shannon = filtered_shannon,
                      unfiltered_richness = unfiltered_richness,
                      filtered_richness = filtered_richness,
                      linages_called = paste0(lineage_names, collapse = "/"),
                      abundances = paste0(abundance_vector, collapse = "/"),
                      freyja_summary = sample$summary))
    
  }
  
  filtered_shannon = -sum(log(abundance_vector) * abundance_vector)
  filtered_richness = length(abundance_vector)
  
  return(data.frame(unfiltered_shannon = unfiltered_shannon,
                    filtered_shannon = filtered_shannon,
                    unfiltered_richness = unfiltered_richness,
                    filtered_richness = filtered_richness,
                    linages_called = paste0(lineage_names, collapse = "/"),
                    abundances = paste0(abundance_vector, collapse = "/"),
                    freyja_summary = sample$summary))
  
}

shannon = as.data.frame((do.call(rbind, (lapply(1:nrow(freyja_df), calculate_shannon_for_1_sample)))))


metrics$filtered_shannon = unlist(shannon$filtered_shannon)
metrics$unfiltered_shannon = unlist(shannon$unfiltered_shannon)
metrics$unfiltered_richness= unlist(shannon$unfiltered_richness)
metrics$filtered_richness= unlist(shannon$filtered_richness)

metrics$freyja_summary = unlist(shannon$freyja_summary)
metrics$lineages_called = unlist(shannon$linages_called)
metrics$abundances = unlist(shannon$abundances)

metrics$Run = metrics$sample_id

metadata = read.delim("/users/sana/Documents/AIM3/021-full_pipeline/00-ALL_FILES/matched_metadata.tsv")

all = merge(metrics, metadata, by = "Run")


### filter all by depth and breadth of coverage:
all$min_50_site_count = as.numeric(all$min_50_site_count)
all$mean_depth_from_bed = as.numeric(all$mean_depth_from_bed)

all = all[-which(all$mean_depth_from_bed < 50), ]
all = all[-which(all$min_50_site_count < 15000), ]

write.table(all, file = 
              "/users/sana/Documents/AIM3/021-full_pipeline/1-metrics/freyja_lineage_details_strictcutoff_filtered.tsv",
            sep = "\t",
            row.names = FALSE)


## create sub table of the samples whose richness is > 1:

subtab = all[which(all$filtered_richness > 1), ]
write.table(subtab, file = 
              "/users/sana/Documents/AIM3/021-full_pipeline/1-metrics/high_richness_freyja_lineage_details_strictcutoff_filtered.tsv",
            sep = "\t",
            row.names = FALSE)
