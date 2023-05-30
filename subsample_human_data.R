### subsample human data from ncbi downloads
#I downloaded all the results with: SARS + REGION + Illumina + FATSQ + Paired (The last three I added in the sidebar)
library(stringr)
setwd("/users/sana/Documents/AIM3/015-human_matches/")

# add_month = function(data){
#   months = substr(data$Collection_Date, 1, 7)
#   data$month = months
#   return(data)
# }
# 
# #The netherlands
# netherlands = read.csv("Netherlands.txt")
# netherlands = add_month(netherlands)
# months = c("2020-03","2020-04", "2020-05", "2020-06", "2020-07", "2020-08")
# sub_netherlands = netherlands[which(netherlands$month %in% months), ]
# netherlands_out = sub_netherlands[sample(1:dim(sub_netherlands)[1])[1:10],]
# write.table(netherlands_out, file = "netherlands_subsample.tsv", sep = "\t", row.names = FALSE)
# 
# write_subsample = function(location, months, n){
#   data = read.csv(paste0(location, ".txt"))
#   data = add_month(data)
#   sub_data = data[which(data$month %in% months), ]
#   print("dim sub data:")
#   print(dim(sub_data))
#   out = sub_data[sample(1:dim(sub_data)[1])[1:n], ]
#   return(out)
# }


write_subsample = function(location){
  data = read.csv(paste0(location, ".txt"))
  write.table(data, file = paste0(location, "_table.tsv"), row.names = FALSE, sep = "\t")
}
write_subsample("Netherlands")
write_subsample("Colorado")
write_subsample("Italy")
write_subsample("Peru")
write_subsample("Madrid")
write_subsample("Kansas")
write_subsample("Quebec")
write_subsample("Ontario")
write_subsample("Ohio")
write_subsample("Iowa")
write_subsample("Greece")

