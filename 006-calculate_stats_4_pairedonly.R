library(stringr)
filenames = read.delim("parsed_names_list.txt", header = FALSE)
filenames = filenames[, 1]
df = data.frame(SRA_ID = vector(),
                enough_avg_depth = vector(),
                average_depth = vector(),
                n_segregating_sites = vector(),
                pi = vector(),
                tajimas_D = vector())

for (n in filenames){
  print(n)
  ### load the table
  tab = read.delim(n)
  avg_dpth = mean(tab$DEPTH)
  
  #flag for average depth
  flag = "y"
  if(avg_dpth < 50){flag = "n"}
  id = unlist(str_split(n, "\\."))[1]
  
  #remove sites with low coverage
  if(length(which(tab$DEPTH < 50)) > 0){tab = tab[-which(tab$DEPTH < 50),]}
  
  if(length(which(tab$VAF > 0.95)) > 0){tab = tab[-which(tab$VAF > 0.95),]}
  
  
  if(dim(tab)[1] == 0){
    df[(nrow(df) + 1), ] = c(id, flag, avg_dpth, 0, "-", "-")
    next}
  print("here")
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
  
  df[(nrow(df) + 1), ] = c(id, flag, avg_dpth, n_seg_sites, pi, tajimas_d)
}



write.table(df, file = "dataframe.tsv", sep = "\t", row.names = FALSE)