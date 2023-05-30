#### load SRA metadata files
cat = read.delim("/Users/sana/Documents/AIM3/001-SRA_metadata/SraRunTable-cat.txt", sep = ',')
dog = read.delim("/Users/sana/Documents/AIM3/001-SRA_metadata/SraRunTable-dog.txt", sep = ',')
mink = read.delim("/Users/sana/Documents/AIM3/001-SRA_metadata/SraRunTable-mink.txt", sep = ',')
deer = read.delim("/Users/sana/Documents/AIM3/001-SRA_metadata/SraRunTable-deer.txt", sep = ',')

cat = cat[which(cat$Run %in% coverage_df$id), ]
dog = dog[which(dog$Run %in% coverage_df$id), ]
deer = deer[which(deer$Run %in% coverage_df$id), ]
mink = mink[which(mink$Run %in% coverage_df$id), ]


write.table(cat, file = "cat.tsv", sep = "\t", row.names = FALSE)
write.table(dog, file = "dog.tsv", sep = "\t", row.names = FALSE)
write.table(deer, file = "deer.tsv", sep = "\t", row.names = FALSE)
write.table(mink, file = "mink.tsv", sep = "\t", row.names = FALSE)
