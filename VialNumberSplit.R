require(xlsx)
king2015 <- read.xlsx(file = "V:/WORK/Chinook/Westward/CSRI Westward Commercial Harvest 2014-2016/KMA Chinook Sample Selection 2015mbf.xlsx", sheetName = "Data", stringsAsFactors = FALSE)
str(king2015)

vials2besplit <- king2015$vial.number[which(king2015$XTR == "X")]
vials2besplit


unlist(strsplit(x = grep(pattern = "; ", x = vials2besplit, value = TRUE), split = "; "))


# Get individual vial numbers by splitting characters
test <- sapply(vials2besplit, function(vial) {
  if(";" %in% unlist(strsplit(x = vial, split = ""))){
    vial <- unlist(strsplit(x = vial, split = "; "))
  }
  if("-" %in% unlist(strsplit(x = vial, split = ""))) {
    vials <- sapply(vial, function(Vial) {as.numeric(unlist(strsplit(x = Vial, split = "-")))}, simplify = FALSE )
    lapply(vials, function(Vials) {seq(from = Vials[1], to = Vials[2], by = 1)})
  } else {
    as.numeric(vial)
  }
} )


# Get number of individuals per row
sampsize <- lapply(test, function(row) {length(unlist(row))})
str(sampsize, max.level = 0)
head(sampsize)

# check weird one which was a collection of two different vial ranges
test[352]

# Number of fish to extract for 2015 mixtures
length(unlist(x = test, recursive = TRUE, use.names = FALSE))

# Vector of vials to extract
sort(unlist(x = test, recursive = TRUE, use.names = FALSE))

# Write tables
king2015xtr <- king2015[which(king2015$XTR == "X"), ]
king2015xtr$SampleSize <- unlist(sampsize, use.names = FALSE)
write.xlsx(x = king2015xtr, file = "V:/WORK/Chinook/Westward/CSRI Westward Commercial Harvest 2014-2016/KMA Chinook Sample Selection 2015mbf.xlsx", sheetName = "NewData", row.names = FALSE, append = TRUE)

write.xlsx(x = sort(unlist(x = test, recursive = TRUE, use.names = FALSE)), file = "V:/WORK/Chinook/Westward/CSRI Westward Commercial Harvest 2014-2016/KMA Chinook Sample Selection 2015mbf.xlsx", sheetName = "ForLAB", row.names = FALSE, append = TRUE)
