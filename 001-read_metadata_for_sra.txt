##
setwd("/Users/sana/Documents/AIM3/SRA_metadata")
tab = read.delim("/Users/sana/Documents/AIM3/SRA_metadata/SraRunTable-cat.txt", sep = ",")

acc = tab$Run
write.table(acc, file = "accession_cat.txt", sep = "\n", quote = FALSE, row.names = FALSE, col.names = FALSE)
