rm(list = ls())
setwd("V:/Analysis/4_Westward/Chinook/CSRI Westward Commercial Harvest 2014-2016/Mixtures/Extraction")
load(file = "V:/Analysis/4_Westward/Chinook/CSRI Westward Commercial Harvest 2014-2016/Mixtures/Github-Westward-Chinook-Mixtures/KMA2016Extraction.RData")

# Had to remove an "X" in the "XTR" column that had "NO VIALS"

require(xlsx)
king2016 <- read.xlsx(file = "WWChinookSampleSelection2016 FINAL_KS.xls", sheetName = "Chinook", startRow = 3, stringsAsFactors = FALSE)
str(king2016)

king2016 <- king2016[, seq(which(dimnames(king2016)[[2]] == "comments"))]
str(king2016)

# What are possible values for the "XTR" column?
table(king2016$XTR)
sum(king2016$XTR == "X", na.rm = TRUE)
sum(king2016$XTR == "X ", na.rm = TRUE)

# Replace incorrect values in "XTR" column
king2016$XTR <- gsub(pattern = "X ", replacement = "X", x = king2016$XTR)

vials2besplit <- king2016$Vial.Number[which(king2016$XTR == "X")]
vials2besplit <- king2016$Vial.Number[1:443]
vials2besplit <- vials2besplit[!vials2besplit %in% "NO VIALS"]


strsplit(x = grep(pattern = ", ", x = vials2besplit, value = TRUE), split = ", ")


# Get individual vial numbers by splitting characters
test <- sapply(vials2besplit, function(vial) {
  if("," %in% unlist(strsplit(x = vial, split = ""))){
    vial <- unlist(strsplit(x = vial, split = ", "))
  }
  if("-" %in% unlist(strsplit(x = vial, split = ""))) {
    vials <- sapply(vial, function(Vial) {as.numeric(unlist(strsplit(x = Vial, split = "-")))}, simplify = FALSE)
    lapply(vials, function(Vials) {seq(from = min(Vials), to = max(Vials), by = 1)})
  } else {
    as.numeric(vial)
  }
} )

str(test, max.level = 0)
test[1:10]
str(test[181])



# Get number of individuals per row
sampsize <- lapply(test, function(row) {length(unlist(row))})
str(sampsize, max.level = 0)
head(sampsize)

# check weird one which was a collection of two different vial ranges
test[181]
sampsize[181]

# Number of fish to extract for 2015 mixtures
length(unlist(x = test, recursive = TRUE, use.names = FALSE))
sum(unlist(sampsize))

sum(unlist(sampsize)) / 95  # need to remove 4 fish


nfish <- as.numeric(readClipboard())
table(unlist(sampsize) == nfish)
which(!unlist(sampsize) == nfish)


# Vector of vials to extract
vials2extract <- unlist(x = test, recursive = TRUE, use.names = FALSE)

# Write tables
king2016xtr <- king2016[which(king2016$XTR == "X"), ]
king2016xtr$SampleSize <- unlist(sampsize, use.names = FALSE)
str(king2016xtr)


# Create data.frame with one row per fish
king2016xtrperfish <- apply(X = king2016xtr, 1, function(row) {matrix(data = rep(row, times = row["SampleSize"]), ncol = 18, byrow = TRUE)} )
king2016xtrperfish.mat <- Reduce(f = rbind, x = king2016xtrperfish)

str(king2016xtrperfish.mat)
king2016xtrperfish.df <- data.frame(king2016xtrperfish.mat, stringsAsFactors = FALSE)
head(king2016xtrperfish.df)
dimnames(king2016xtrperfish.df)[[2]] <- dimnames(king2016xtr)[[2]]
str(king2016xtrperfish.df)

king2016xtrperfish.df$Vial.Number <- vials2extract

# Pivot
table(king2016xtrperfish.df$Area, king2016xtrperfish.df$Strata)

# Remove 4 fish to get even plate numbers, one from each of the strata with 380 -> 379

vials2remove <- c(
  sample(x = subset(x = king2016xtrperfish.df, subset = Area == "Mainland" & Strata == "E")$Vial.Number, size = 1),
  sample(x = subset(x = king2016xtrperfish.df, subset = Area == "Mainland" & Strata == "L")$Vial.Number, size = 1),
  sample(x = subset(x = king2016xtrperfish.df, subset = Area == "Westside" & Strata == "E")$Vial.Number, size = 1),
  sample(x = subset(x = king2016xtrperfish.df, subset = Area == "Westside" & Strata == "L")$Vial.Number, size = 1)
)

length(sort(vials2extract[!vials2extract %in% vials2remove]))

king2016xtrperfish.df.final <- king2016xtrperfish.df[-which(king2016xtrperfish.df$Vial.Number %in% vials2remove), ]


# Pivot
table(king2016xtrperfish.df.final$Area, king2016xtrperfish.df.final$Strata)

dim(king2016xtrperfish.df.final)[1] / 29


# Convert to date
king2016xtrperfish.df.final$single.catch.date <- as.Date(king2016xtrperfish.df.final$single.catch.date)

str(king2016xtr)
str(king2016xtrperfish.df.final)
# Write tables

write.xlsx(x = king2016xtr, file = "WWChinookSampleSelection2016 FINAL_KS.xls", sheetName = "NewData", row.names = FALSE, append = TRUE)
write.xlsx(x = king2016xtrperfish.df.final, file = "WWChinookSampleSelection2016 FINAL_KS.xls", sheetName = "NewDataOneFishPerRow", row.names = FALSE, append = TRUE)
write.xlsx(x = sort(vials2extract[!vials2extract %in% vials2remove]), file = "WWChinookSampleSelection2016 FINAL_KS.xls", sheetName = "KKMAC16 ForLAB", row.names = FALSE, append = TRUE)

save.image(file = "KMA2016Extraction.RData")










require(xlsx)
king2016 <- read.xlsx(file = "WWChinookSampleSelection2016 FINAL_KS.xls", sheetName = "Chinook", startRow = 3, stringsAsFactors = FALSE)
str(king2016)

which(king2016$Vial.Number == "NO VIALS")
king2016 <- king2016[c(1:389, 391:443), seq(which(dimnames(king2016)[[2]] == "comments"))]
str(king2016)


vials2besplit <- king2016$Vial.Number



strsplit(x = grep(pattern = ", ", x = vials2besplit, value = TRUE), split = ", ")


# Get individual vial numbers by splitting characters
test <- sapply(vials2besplit, function(vial) {
  if("," %in% unlist(strsplit(x = vial, split = ""))){
    vial <- unlist(strsplit(x = vial, split = ", "))
  }
  if("-" %in% unlist(strsplit(x = vial, split = ""))) {
    vials <- sapply(vial, function(Vial) {as.numeric(unlist(strsplit(x = Vial, split = "-")))}, simplify = FALSE)
    lapply(vials, function(Vials) {seq(from = min(Vials), to = max(Vials), by = 1)})
  } else {
    as.numeric(vial)
  }
} )

str(test, max.level = 0)
test[1:10]
str(test[181])



# Get number of individuals per row
sampsize <- lapply(test, function(row) {length(unlist(row))})
str(sampsize, max.level = 0)
head(sampsize)

# check weird one which was a collection of two different vial ranges
test[181]
sampsize[181]

# Number of fish to extract for 2015 mixtures
length(unlist(x = test, recursive = TRUE, use.names = FALSE))
sum(unlist(sampsize))


# Write tables
king2016$SampleSize <- unlist(sampsize, use.names = FALSE)
str(king2016)


# Create data.frame with one row per fish
king2016perfish <- apply(X = king2016, 1, function(row) {matrix(data = rep(row, times = row["SampleSize"]), ncol = 18, byrow = TRUE)} )
king2016perfish.mat <- Reduce(f = rbind, x = king2016perfish)

str(king2016perfish.mat)
king2016perfish.df <- data.frame(king2016perfish.mat, stringsAsFactors = FALSE)
head(king2016perfish.df)
dimnames(king2016perfish.df)[[2]] <- dimnames(king2016)[[2]]
str(king2016perfish.df)

king2016perfish.df$VialNumber <- unlist(test, recursive = TRUE, use.names = FALSE)

king2016perfish.df.final <- king2016perfish.df[, c("VialNumber", "Area", "single.catch.date", "Catch.Date", "Card", "Strata", "Sampling.Port")]
str(king2016perfish.df.final)


write.xlsx(x = king2016perfish.df.final, file = "WWChinookSampleSelection2016 FINAL_KS.xls", sheetName = "AllVialsNewDataOneFishPerRow", row.names = FALSE, append = TRUE)
