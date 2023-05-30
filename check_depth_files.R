#load deptj names
dep = read.delim("/Users/sana/Documents/AIM3/007-all_runs_depth_summary_files/depth_file_list.tsv", header = FALSE)
dep = dep[, 1]

parsed = read.delim("/Users/sana/Documents/AIM3/006-all_runs_variant_files/1-pool/parsed_names.tsv", header = FALSE)
parsed = parsed[, 1]
 
noseg = read.delim("/Users/sana/Documents/AIM3/006-all_runs_variant_files/1-pool/noseg_names.tsv", header = FALSE)
noseg = noseg[, 1]

var = c(noseg, parsed)

var = str_split(var, "\\.")
dep = str_split(dep, "\\.")

dep_vec = vector()
var_vec = vector()
for(i in 1:length(var)){
  var_vec = c(var_vec, unlist(var[i])[1])
}
for(i in 1:length(dep)){
  dep_vec = c(dep_vec, unlist(dep[i])[1])
}



readset = read.delim("/Users/sana/Documents/AIM3/006-all_runs_variant_files/readset_for_paired_and_single.tsv")
r = readset[which(readset$Sample %in% var_vec[which(!var_vec %in% dep_vec)]), ]
write.table(r, file = "/Users/sana/Documents/AIM3/r.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
