## script to load every vcf sample one by one and check the frequency of GWAS hits in each one
rm(list = ls())
setwd("/Users/sana/Documents/AIM3/021-full_pipeline/5-gwas_comparison/")
library(stringr)

#threshold for minimum number of sites with enough coverage
thresh = 15000
minimum_depth = 50

metadata = read.delim("/users/sana/Documents/AIM3/021-full_pipeline/1-metrics/metrics.tsv")

mink = read.delim("mink.tsv")
mink$species = "mink"
print("ignore warning for loading mink table")
#ignore warning
deer = read.delim("deer.tsv")
deer$species = "deer"

hits = rbind(mink, deer)
head(hits)

####### load the dataframe with coverage report
coverage_report = read.delim("/users/sana/Documents/AIM3/019-GWAS_comparison/covreport.tsv")
## add file ame column to the coverage report
names = str_split(coverage_report$filename, "\\.")
coverage_report$SRA_ID = "tbd"
for(i in 1:length(names)){
  coverage_report$SRA_ID[i] = unlist(names[i])[1]
}
length(which(coverage_report$SRA_ID %in% metadata$Run))

coverage_report = coverage_report[which(coverage_report$SRA_ID %in% metadata$Run), ]
coverage_report = coverage_report[which(coverage_report$min_50_sites >= thresh), ]

metadata = metadata[which(metadata$Run %in% coverage_report$SRA_ID), ]

#load all file names
parsed_names = read.delim("/Users/sana/Documents/AIM3/021-full_pipeline/5-gwas_comparison/parsednames", header = FALSE)
parsed_names = parsed_names[, 1]
parsed_names_df = data.frame(filename = vector(),
                             samplename = vector())
for(i in 1:length(parsed_names)){
  parsed_names_df[i, ] = c(parsed_names[i],
                           unlist(str_split(parsed_names[i], "\\."))[1])
}

parsed_names = parsed_names[which(parsed_names_df$samplename %in% coverage_report$SRA_ID)]


dim(metadata)
create_list_of_hicov = function(species){
  sra_list = metadata$Run[which(metadata$Species == species)]
  coverage_files = paste0(sra_list, ".sorted.filtered.primerTrim.bam.bed.coverage")
  first = read.delim(coverage_files[1], header = FALSE)
  high_cov_sites = which(first$V4 >= minimum_depth)
  for(i in 2:length(coverage_files)){
    covfile = read.delim(coverage_files[i], header = FALSE)
    high_cov_sites = intersect(high_cov_sites, which(covfile$V4 > minimum_depth))
    #print(length(high_cov_sites))
  }
  return(high_cov_sites)
}

deer_hicov_list = create_list_of_hicov("deer")
mink_hicov_list = create_list_of_hicov("mink")




###make sure to check typeof hit positions for comparison
tab_test = read.delim(parsed_names[1])
if(!typeof(tab_test$POS) == typeof(hits$Pos.)){
  stop("data types not okay")
}

check_for_hit = function(hit_pos, filename){
  tab = read.delim(filename)
  if(length(which(tab$DEPTH < minimum_depth) > 0)){
    tab = tab[-which(tab$DEPTH < minimum_depth), ]
  }
  if(nrow(tab) == 0){
    return("not-found")
  }
  if(hit_pos %in% tab$POS){
    return(tab$INFO[which(tab$POS == hit_pos)])
  }
  if(!hit_pos %in% tab$POS){
    return("not-found")
  }
}


df = data.frame()
for(i in parsed_names){
  sample_name = unlist(str_split(i, "\\."))[1]
  df = rbind(df, c(sample_name, 
                   unlist(lapply(hits$Pos., check_for_hit, i))))
  colnames(df) = c("sample_name", as.character(hits$Pos.))
}

df$species = "tbd"
length(which(df$sample_name %in% metadata$Run))
df$species = metadata$Species[match(df$sample_name, metadata$Run)]
gwas_hits_per_sample = df[order(df$species), ]

write.table(df, file = "gwas_hits_per_sample_cutoff.tsv", sep = "\t", row.names = FALSE)


### summarize the counts in another table:


###load all data from script 18
mutations = colnames(gwas_hits_per_sample)[2:30]

df = data.frame(deer_notfound = vector(),
                deer_fixed = vector(),
                deer_segregating = vector(),
                mink_notfound = vector(),
                mink_fixed = vector(),
                mink_segregating = vector(),
                cat_notfound = vector(),
                cat_fixed = vector(),
                cat_segregating = vector(),
                dog_notfound = vector(),
                dog_fixed = vector(),
                dog_segregating = vector())

for(mut in mutations){
  vec = vector()
  for(species in c("deer", "mink", "cat", "dog")){
    col = gwas_hits_per_sample[which(gwas_hits_per_sample$species == species), which(colnames(gwas_hits_per_sample) == mut)]
    vec = c(vec, length(which(col == "not-found")),
            length(which(str_detect(col, "fixed"))),
            length(which(str_detect(col, "ambiguous"))))
  }
  vec = data.frame(matrix(data = vec, nrow = 1))
  colnames(vec) = colnames(df)
  df = rbind(df, vec)
}


df$mutations = str_remove_all(mutations, "X")

write.table(df, file = "counts_cutoff_new.tsv", sep = "\t", row.names = FALSE)








################# count the number of occurences for every site in every sample
library(phylotools)
ref = read.fasta("ref.MN908947.3.fasta")
sites = 1:nchar(ref$seq.text)
all_sites_matrix = data.frame(matrix(data = 0, 
                                     nrow = length(sites), ncol = nrow(gwas_hits_per_sample)))
colnames(all_sites_matrix) = gwas_hits_per_sample$sample_name

file_names = paste0(gwas_hits_per_sample$sample_name, ".", "freebayes_calling.consensus.vcf_table.tsv.unfiltered.parsed")


#### this secrion creates all_sites_matrix, nrow = all 29903 sites ncol = samples
### each element i,j says whether site i exists in sample j

check_for_hit_readonce_child = function(hit_pos, tab){
  if(nrow(tab) == 0){
    return("not-found")
  }
  if(hit_pos %in% tab$POS){
    return(tab$INFO[which(tab$POS == hit_pos)])
  }
  if(!hit_pos %in% tab$POS){
    return("not-found")
  }
}
check_for_hit_readonce_parent = function(filename){
  tab = read.delim(filename)
  if(length(which(tab$DEPTH < minimum_depth) > 0)){
    tab = tab[-which(tab$DEPTH < minimum_depth), ]
  }
  arr = unlist(lapply(sites, check_for_hit_readonce_child, tab))
  return(arr)
}

for(i in 1:length(file_names)){
  print(file_names[i])
  all_sites_matrix[, i] = check_for_hit_readonce_parent(file_names[i])
  
}


#### this matrix counts the number of each tag for every site
all_sites_counts = data.frame(deer_notfound = vector(),
                              deer_fixed = vector(),
                              deer_segregating = vector(),
                              mink_notfound = vector(),
                              mink_fixed = vector(),
                              mink_segregating = vector(),
                              cat_notfound = vector(),
                              cat_fixed = vector(),
                              cat_segregating = vector(),
                              dog_notfound = vector(),
                              dog_fixed = vector(),
                              dog_segregating = vector())


colnames_species = metadata$Species[match(colnames(all_sites_matrix), 
                                          metadata$Run)]

for(mut in sites){
  print(mut)
  vec = vector()
  for(species in c("deer", "mink", "cat", "dog")){
    col = all_sites_matrix[mut, which(colnames_species == species)]
    vec = c(vec, length(which(col == "not-found")),
            length(which(str_detect(col, "fixed"))),
            length(which(str_detect(col, "ambiguous"))))
  }
  vec = data.frame(matrix(data = vec, nrow = 1))
  colnames(vec) = colnames(all_sites_counts)
  all_sites_counts = rbind(all_sites_counts, vec)
}

all_sites_counts$gwas_hit = "non-GWAS-site"
all_sites_counts$gwas_hit[hits$Pos.[which(hits$species == "deer")]] = "deer-GWAS-hit"
all_sites_counts$gwas_hit[hits$Pos.[which(hits$species == "mink")]] = "mink-GWAS-hit"

#write.table(all_sites_counts, file = "all_sites_counts_50_cutoff.tsv", sep = "\t", row.names = FALSE)



#setwd("/users/sana/Documents/AIM3/019-GWAS_comparison/boxplots_for_all_sites/plots_4_filtered_sites/")
library(ggplot2)
s = 22
deer_segregating = ggplot(data = all_sites_counts, aes(x = gwas_hit, y = deer_segregating)) + 
  ggtitle("a.")+ ylab("Number of iSNVs in deer") + 
  theme_bw() + 
  geom_boxplot( fill = "#4A5899", alpha = 0.3, outlier.shape = NA) +
  geom_jitter(width = 0.25, height=0, size = 1, stroke = 1, color = "#35406E")+
  ylim(c(0,70))+ xlab("")+
  theme(axis.text.x = element_text(angle = 90,hjust = 0.8),
        text = element_text(size=s),
        plot.title = element_text(size = 15, face = "bold"))
ggsave(plot = deer_segregating, filename = "deer_segregating.pdf", device = "pdf")


deer_fixed = ggplot(data = all_sites_counts, aes(x = gwas_hit, y = deer_fixed)) + 
  ggtitle("b.") + ylab("Deer SNVs fixed between hosts") + 
  theme_bw() + geom_boxplot( fill = "#4A5899", alpha = 0.3, outlier.shape = NA)+
  geom_jitter(width = 0.25, height=0, size = 1, stroke = 1, color = "#35406E")+
  ylim(c(0,70)) + xlab("")+
  theme(axis.text.x = element_text(angle = 90,hjust = 0.8),
        text = element_text(size=s),
        plot.title = element_text(size = 15, face = "bold"))
ggsave(plot = deer_fixed, filename = "deer_fixed.pdf", device = "pdf")


mink_segregating = ggplot(data = all_sites_counts, aes(x = gwas_hit, y = mink_segregating)) + 
  ggtitle("c.") + ylab("Number of iSNVs in mink") + 
  theme_bw() + geom_boxplot( fill = "#4A5899", alpha = 0.3, outlier.shape = NA)+
  geom_jitter(width = 0.25, height=0, size = 1, stroke = 1, color = "#35406E")+
  ylim(c(0,70)) +xlab("")+
  theme(axis.text.x = element_text(angle = 90,hjust = 0.8),
        text = element_text(size=s),
        plot.title = element_text(size = 15, face = "bold"))
ggsave(plot = mink_segregating, filename = "mink_segregating.pdf", device = "pdf")

mink_fixed = ggplot(data = all_sites_counts, aes(x = gwas_hit, y = mink_fixed)) + 
  ggtitle("d.")+ ylab("Mink SNVs fixed between hosts") + 
  theme_bw() + geom_boxplot( fill = "#4A5899", alpha = 0.3, outlier.shape = NA)+
  ylim(c(0,70)) +xlab("")+
  geom_jitter(width = 0.25, height=0, size = 1, stroke = 1, color = "#35406E")+
  theme(axis.text.x = element_text(angle = 90,hjust = 0.8),
        text = element_text(size=s),
        plot.title = element_text(size = 15, face = "bold"))
ggsave(plot = mink_fixed, filename = "mink_fixed.pdf", device = "pdf")

library(ggpubr)
all = ggarrange(deer_segregating, deer_fixed,
                mink_segregating, mink_fixed, nrow = 1, ncol = 4)
ggsave(filename = "all_GWAS_counts.pdf", plot = all, device = "pdf", height = 10, width = 15)
ggsave(filename = "/users/sana/Documents/AIM3/021-full_pipeline/xx-Manuscript/figure4_final.pdf", plot = all, device = "pdf", height = 10, width = 15)



##### fishers exact test ############

#list of non-gwas sites on the genome
#non_gwas_sites = which(all_sites_counts$gwas_hit == "non-GWAS-hit")


##### compare deer segregating counts
non_gwas_sites = deer_hicov_list[-which(deer_hicov_list %in% hits$Pos.)]
compare_each_gwas_w_nongwas_deer_segregating = function(gwas_hit){
  print(gwas_hit)
  s = vector()
  for(nongwas in non_gwas_sites){
    fishers_matrix = matrix(data = c(
      all_sites_counts$deer_segregating[gwas_hit],
      all_sites_counts$deer_segregating[nongwas],
      all_sites_counts$deer_notfound[gwas_hit],
      all_sites_counts$deer_notfound[nongwas]
    ), nrow = 2, ncol = 2, byrow = TRUE)
    colnames(fishers_matrix) = c("GWAS", "non-GWAS")
    rownames(fishers_matrix) = c("segregating", "not-found")
    tstat = fisher.test(fishers_matrix)
    if(tstat$p.value < 0.05){
      if((fishers_matrix[1,1]/fishers_matrix[2,1]) > (fishers_matrix[1,2]/fishers_matrix[2,2])){
        s = c(s, nongwas)
      }
    }
  }
  return(length(s))
}


deer_counts_segregating = unlist(lapply(hits$Pos.[which(hits$Pos. %in% deer_hicov_list)], 
                                        compare_each_gwas_w_nongwas_deer_segregating))/length(non_gwas_sites)


##### compare mink segregating counts 
non_gwas_sites = mink_hicov_list[-which(mink_hicov_list %in% hits$Pos.)]
compare_each_gwas_w_nongwas_mink_segregating = function(gwas_hit){
  print(gwas_hit)
  s = vector()
  for(nongwas in non_gwas_sites){
    fishers_matrix = matrix(data = c(
      all_sites_counts$mink_segregating[gwas_hit],
      all_sites_counts$mink_segregating[nongwas],
      all_sites_counts$mink_notfound[gwas_hit],
      all_sites_counts$mink_notfound[nongwas]
    ), nrow = 2, ncol = 2, byrow = TRUE)
    colnames(fishers_matrix) = c("GWAS", "non-GWAS")
    rownames(fishers_matrix) = c("segregating", "not-found")
    tstat = fisher.test(fishers_matrix)
    if(tstat$p.value < 0.05){
      if((fishers_matrix[1,1]/fishers_matrix[2,1]) > (fishers_matrix[1,2]/fishers_matrix[2,2])){
        s = c(s, nongwas)
      }
    }
  }
  return(length(s))
}

mink_counts_segregating = unlist(lapply(hits$Pos.[which(hits$Pos. %in% mink_hicov_list)],
                                        compare_each_gwas_w_nongwas_mink_segregating))/length(non_gwas_sites)


##### compare deer fixed counts 
non_gwas_sites = deer_hicov_list[-which(deer_hicov_list %in% hits$Pos.)]
compare_each_gwas_w_nongwas_deer_fixed = function(gwas_hit){
  print(gwas_hit)
  s = vector()
  for(nongwas in non_gwas_sites){
    fishers_matrix = matrix(data = c(
      all_sites_counts$deer_fixed[gwas_hit],
      all_sites_counts$deer_fixed[nongwas],
      all_sites_counts$deer_notfound[gwas_hit],
      all_sites_counts$deer_notfound[nongwas]
    ), nrow = 2, ncol = 2, byrow = TRUE)
    colnames(fishers_matrix) = c("GWAS", "non-GWAS")
    rownames(fishers_matrix) = c("segregating", "not-found")
    tstat = fisher.test(fishers_matrix)
    if(tstat$p.value < 0.05){
      if((fishers_matrix[1,1]/fishers_matrix[2,1]) > (fishers_matrix[1,2]/fishers_matrix[2,2])){
        s = c(s, nongwas)
      }
    }
  }
  return(length(s))
}
deer_counts_fixed = unlist(lapply(hits$Pos.[which(hits$Pos. %in% deer_hicov_list)], 
                                  compare_each_gwas_w_nongwas_deer_fixed))/length(non_gwas_sites)

##### compare mink fixed counts 
non_gwas_sites = mink_hicov_list[-which(mink_hicov_list %in% hits$Pos.)]
compare_each_gwas_w_nongwas_mink_fixed = function(gwas_hit){
  print(gwas_hit)
  s = vector()
  for(nongwas in non_gwas_sites){
    fishers_matrix = matrix(data = c(
      all_sites_counts$mink_fixed[gwas_hit],
      all_sites_counts$mink_fixed[nongwas],
      all_sites_counts$mink_notfound[gwas_hit],
      all_sites_counts$mink_notfound[nongwas]
    ), nrow = 2, ncol = 2, byrow = TRUE)
    colnames(fishers_matrix) = c("GWAS", "non-GWAS")
    rownames(fishers_matrix) = c("segregating", "not-found")
    tstat = fisher.test(fishers_matrix)
    if(tstat$p.value < 0.05){
      if((fishers_matrix[1,1]/fishers_matrix[2,1]) > (fishers_matrix[1,2]/fishers_matrix[2,2])){
        s = c(s, nongwas)
      }
    }
  }
  return(length(s))
}

mink_counts_fixed = unlist(lapply(hits$Pos.[which(hits$Pos. %in% mink_hicov_list)],
                                  compare_each_gwas_w_nongwas_mink_fixed))/length(non_gwas_sites)



mink_counts_fixed
mink_counts_segregating
deer_counts_fixed
deer_counts_segregating

setwd("/users/sana/Documents/AIM3/021-full_pipeline/5-gwas_comparison/")
deer_df = data.frame(hits = hits$Pos.[which(hits$Pos. %in% deer_hicov_list)],
                     deer_fixed_percentage = deer_counts_fixed,
                     deer_segregating_percentage = deer_counts_segregating)

mink_df = data.frame(hits = hits$Pos.[which(hits$Pos. %in% mink_hicov_list)],
                     mink_fixed_percentage = mink_counts_fixed,
                     mink_segregating_percentage = mink_counts_segregating)

write.table(deer_df, file = "deer_sample_fisher_results.tsv", sep = "\t", row.names = FALSE)
write.table(mink_df, file = "mink_sample_fisher_results.tsv", sep = "\t", row.names = FALSE)

