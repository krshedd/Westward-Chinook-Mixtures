dates <- readClipboard()
dates

# How many different dates do we have?
unique(dates)

# What are the weird characters in here?
grep(pattern = "//", x = dates)
grep(pattern = ",", x = dates)

# Replace weird characters
dates <- gsub(pattern = "//", replacement = "/", x = dates)
dates <- gsub(pattern = ", ", replacement = "-", x = dates)

# Split based on -
twodates <- grep(pattern = "-", x = dates)

twodatesmat <- sapply(strsplit(x = dates[twodates], split = "-"), function(ind) {
  Month <- substr(x = ind[1], start = 1, stop = 1)
  Year <- "2014"
  c(paste(ind[1], Year, sep = "/"), paste(Month, ind[2], sep = "/"))
} )

newdates <- matrix(data = NA, nrow = length(dates), ncol = 2, dimnames = list(seq_along(dates), c("Begin", "End")))

newdates[twodates, ] <- t(twodatesmat)
newdates[-twodates, 1] <- dates[-twodates]
newdates[-twodates, 2] <- dates[-twodates]

newdates

table(newdates[, 1] == newdates[, 2])

require(xlsx)
write.xlsx(x = newdates, file = "V:/WORK/Chinook/Westward/CSRI Westward Commercial Harvest 2014-2016/KMA Chinook Genetics 2014 Extractions.xlsx", sheetName = "NewDates", col.names = TRUE, row.names = FALSE, append = TRUE)




## 2014 Extraction list
require(xlsx)
kin2014 <- read.xlsx(file = "V:/WORK/Chinook/Westward/CSRI Westward Commercial Harvest 2014-2016/KMA Chinook Genetics 2014 Extractions.xlsx", sheetName = "2014 KMA Chinook", stringsAsFactors = FALSE)
str(king2014)

# Birch's original selection
table(king2014$SAMPLING.AREA, king2014$Stratum, king2014$XTR)
OriginalVials2XTR <- sort(king2014$VIAL[which(king2014$XTR == "X")])


## Need to remove 5 fish from 2014 mixtures for plate math
## Briefly, there are 11 "extra" fish in 2014 and 2015 has "space" for 6, the other five will come from 2014 strata that have 380 fish (make them 379)

# Remove 1 fish randomly from Westside Early 2014
WestEarly2014Vials <- king2014$VIAL[which(king2014$XTR == "X" & king2014$Stratum == "Early" & king2014$SAMPLING.AREA == "Westside")]
length(WestEarly2014Vials)

WestEarly2014Vials2XTR <- sample(x = WestEarly2014Vials, size = 379, replace = FALSE)

WestEarly2014Vials2RMV <- WestEarly2014Vials[!WestEarly2014Vials %in% WestEarly2014Vials2XTR]

# Remove 1 fish randomly from Westside Late 2014
WestLate2014Vials <- king2014$VIAL[which(king2014$XTR == "X" & king2014$Stratum == "Late" & king2014$SAMPLING.AREA == "Westside")]
length(WestLate2014Vials)

WestLate2014Vials2XTR <- sample(x = WestLate2014Vials, size = 379, replace = FALSE)

WestLate2014Vials2RMV <- WestLate2014Vials[!WestLate2014Vials %in% WestLate2014Vials2XTR]

# Remove 1 fish randomly from Eastside Late 2014
EastLate2014Vials <- king2014$VIAL[which(king2014$XTR == "X" & king2014$Stratum == "Late" & king2014$SAMPLING.AREA == "Eastside")]
length(EastLate2014Vials)

EastLate2014Vials2XTR <- sample(x = EastLate2014Vials, size = 379, replace = FALSE)

EastLate2014Vials2RMV <- EastLate2014Vials[!EastLate2014Vials %in% EastLate2014Vials2XTR]

# Remove 1 fish randomly from Mainland Late 2014
MainLate2014Vials <- king2014$VIAL[which(king2014$XTR == "X" & king2014$Stratum == "Late" & king2014$SAMPLING.AREA == "Mainland")]
length(MainLate2014Vials)

MainLate2014Vials2XTR <- sample(x = MainLate2014Vials, size = 379, replace = FALSE)

MainLate2014Vials2RMV <- MainLate2014Vials[!MainLate2014Vials %in% MainLate2014Vials2XTR]

# Remove 1 fish randomly from SW Kodiak/Alitak Late 2014
SWAlitakLate2014Vials <- king2014$VIAL[which(king2014$XTR == "X" & king2014$Stratum == "Late" & king2014$SAMPLING.AREA == "SW Kodiak/Alitak")]
length(SWAlitakLate2014Vials)

SWAlitakLate2014Vials2XTR <- sample(x = SWAlitakLate2014Vials, size = 379, replace = FALSE)

SWAlitakLate2014Vials2RMV <- SWAlitakLate2014Vials[!SWAlitakLate2014Vials %in% SWAlitakLate2014Vials2XTR]


# Vials to remove
Vials2RMV <- sapply(objects(pattern = "RMV"), get)


FinalVials2XTR <- OriginalVials2XTR[!OriginalVials2XTR %in% Vials2RMV]

length(FinalVials2XTR)

length(OriginalVials2XTR) - length(FinalVials2XTR)

# New Selection
table(king2014$SAMPLING.AREA[king2014$VIAL %in% FinalVials2XTR], king2014$Stratum[king2014$VIAL %in% FinalVials2XTR], king2014$XTR[king2014$VIAL %in% FinalVials2XTR])


# Write extraction list ForLAB
write.xlsx(x = FinalVials2XTR, file = "V:/WORK/Chinook/Westward/CSRI Westward Commercial Harvest 2014-2016/KMA Chinook Genetics 2014 Extractions.xlsx", sheetName = "ForLAB", col.names = FALSE, row.names = FALSE, append = TRUE)




# Double check a table of what is being extracted
Extraction2014 <- as.numeric(readClipboard())
Extraction2015 <- as.numeric(readClipboard())

Extraction2014 <- Extraction2014[!is.na(Extraction2014)]
Extraction2015 <- Extraction2015[!is.na(Extraction2015)]

length(Extraction2014); length(Extraction2015)


BirchExtraction2014 <- read.xlsx(file = "V:/WORK/Chinook/Westward/CSRI Westward Commercial Harvest 2014-2016/KMA Chinook Genetics 2014 Extractions.xlsx", sheetName = "for R", stringsAsFactors = FALSE)
str(BirchExtraction2014)

length(BirchExtraction2014$VIAL)

table(table(BirchExtraction2014$VIAL))
which(table(BirchExtraction2014$VIAL) == 2)
table(table(Extraction2014))

# We are extracting these
sort(Extraction2014[!Extraction2014 %in% BirchExtraction2014$VIAL[which(BirchExtraction2014$XTR == "X")]])
sort(Extraction2014[!match(Extraction2014, BirchExtraction2014$VIAL[which(BirchExtraction2014$XTR == "X")])])
# We are not extracting these
sort(BirchExtraction2014$VIAL[!BirchExtraction2014$VIAL[which(BirchExtraction2014$XTR == "X")] %in% Extraction2014])
sort(BirchExtraction2014$VIAL[!match(BirchExtraction2014$VIAL[which(BirchExtraction2014$XTR == "X")], Extraction2014)])

sort(Extraction2014)



length(sort(Extraction2014))
length(sort(BirchExtraction2014$VIAL[which(BirchExtraction2014$XTR == "X")]))

writeClipboard(as.character(sort(Extraction2014)))
writeClipboard(as.character(sort(BirchExtraction2014$VIAL[which(BirchExtraction2014$XTR == "X")])))



vial.indexes <- match(Extraction2014, BirchExtraction2014$VIAL)
table(BirchExtraction2014$VIAL[vial.indexes] %in% Extraction2014)

table(BirchExtraction2014$SAMPLING.AREA[vial.indexes], BirchExtraction2014$Stratum[vial.indexes])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Add tissue data for LOKI ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("V:/Analysis/4_Westward/Chinook/CSRI Westward Commercial Harvest 2014-2016/Mixtures")
rm(list = ls())

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### KSPENC14 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in OceanAK data as data.table (lightning fast!)
require(data.table)
oceanak.dt <- fread(input = "OceanAK Tissue Dump/KSPENC14/GEN_SAMPLED_FISH_TISSUE.csv")  # amazing
str(oceanak.dt)
# Convert to data.frame
oceanak.df <- data.frame(oceanak.dt)
str(oceanak.df)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in collection raw datasheet from Birch
require(xlsx)
Spen_Chignik_Datasheet.df <- read.xlsx(file = "South Pen 2014/Extraction/WW Chinook Genetics 2014.xlsx", sheetName = "2014 Pen-Chig Chinook")
str(Spen_Chignik_Datasheet.df)



# Do the oceanak fish exist in datasheet?
sum(oceanak.df$FK_FISH_ID %in% Spen_Chignik_Datasheet.df$VIAL); dim(oceanak.df)[1]

# Do the datasheet fish exist in the oceanak fish?
table(Spen_Chignik_Datasheet.df$VIAL %in% oceanak.df$FK_FISH_ID); dim(oceanak.df)[1]

# Are there duplicates in the datasheet?
table(table(Spen_Chignik_Datasheet.df$VIAL))

names(which(table(Spen_Chignik_Datasheet.df$VIAL) == 2))



# Subset fish we want for date and area
Data_Area.df <- Spen_Chignik_Datasheet.df[Spen_Chignik_Datasheet.df$VIAL %in% oceanak.df$FK_FISH_ID, c("VIAL", "DATE.HARVESTED", "SAMPLING.AREA")]
str(Data_Area.df)

Data_Area.df$DATE.HARVESTED <- as.character(Data_Area.df$DATE.HARVESTED)
Data_Area.df$SAMPLING.AREA <- as.character(Data_Area.df$SAMPLING.AREA)

# Fix area
unique(Data_Area.df$SAMPLING.AREA)
Data_Area.df$SAMPLING.AREA <- gsub(pattern = "ern  ", replacement = "ern", x = Data_Area.df$SAMPLING.AREA)
unique(Data_Area.df$SAMPLING.AREA)

# Fix date
sum(is.na(Data_Area.df$DATE.HARVESTED))

x <- unique(Data_Area.df$DATE.HARVESTED)

as.Date(40336, origin = "1904-01-01", "%m/%d/%Y")
rep(format(x = as.Date(40336, origin = "1904-01-01"), "%m/%d/%Y"), 2)


rep(as.character(as.Date(40336, origin = "1904-01-01")), 2)


y <- grep(pattern = "-", x = x, value = TRUE)

unlist(strsplit(x = y[1], split = "-"))
unlist(strsplit(x = y[1], split = "-|/"))

sapply(y, function(i) {
  i.temp = unlist(strsplit(x = i, split = "-|/"))
  dates <- paste(i.temp[length(i.temp)], i.temp[1], 2014, sep = "-")
  format(c(as.Date(dates, "%b-%d-%Y"), as.Date(dates, "%b-%d-%Y")+1), "%m/%d/%Y")
  } )


as.Date("Jun-20", "%b-%d")


newdates <- t(sapply(Data_Area.df$DATE.HARVESTED, function(dat) {
  if(is.na(dat)) {
    rep(NA, 2)
  } else {
    if(grepl(pattern = "40", dat)) {
      rep(format(as.Date(as.numeric(dat), origin = "1904-01-01"), "%m/%d/%Y"), 2)
    } else {
      i.temp = unlist(strsplit(x = dat, split = "-|/"))
      dates <- paste(i.temp[length(i.temp)], i.temp[1], 2014, sep = "-")
      format(c(as.Date(dates, "%b-%d-%Y"), as.Date(dates, "%b-%d-%Y")+1), "%m/%d/%Y")
    }
  }
}))

str(newdates)


Data_Area.df$Start.Date <- newdates[, 1]
Data_Area.df$End.Date <- newdates[, 2]


oceanak.df$CAPTURE_DATE <- Data_Area.df$Start.Date[match(oceanak.df$FK_FISH_ID, Data_Area.df$VIAL)]
oceanak.df$END_CAPTURE_DATE <- Data_Area.df$End.Date[match(oceanak.df$FK_FISH_ID, Data_Area.df$VIAL)]
oceanak.df$CAPTURE_LOCATION <- Data_Area.df$SAMPLING.AREA[match(oceanak.df$FK_FISH_ID, Data_Area.df$VIAL)]


head(oceanak.df)

str(oceanak.df)

# Replace NAs
oceanak.df[is.na(oceanak.df)] <- ""
dimnames(oceanak.df)[[2]][1] <- "FK_COLLECTION_ID"

write.csv(x = oceanak.df, file = "OceanAK Tissue Dump/KSPENC14/GEN_SAMPLED_FISH_TISSUE Upload.csv", row.names = FALSE, quote = FALSE)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### KCHIGC14 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in OceanAK data as data.table (lightning fast!)
require(data.table)
oceanak.dt <- fread(input = "OceanAK Tissue Dump/KCHIGC14/GEN_SAMPLED_FISH_TISSUE.csv")  # amazing
str(oceanak.dt)
# Convert to data.frame
oceanak.df <- data.frame(oceanak.dt)
str(oceanak.df)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in collection raw datasheet from Birch
require(xlsx)
Spen_Chignik_Datasheet.df <- read.xlsx(file = "South Pen 2014/Extraction/WW Chinook Genetics 2014.xlsx", sheetName = "2014 Pen-Chig Chinook")
str(Spen_Chignik_Datasheet.df)



# Do the oceanak fish exist in datasheet?
sum(oceanak.df$FK_FISH_ID %in% Spen_Chignik_Datasheet.df$VIAL); dim(oceanak.df)[1]

# Do the datasheet fish exist in the oceanak fish?
table(Spen_Chignik_Datasheet.df$VIAL %in% oceanak.df$FK_FISH_ID); dim(oceanak.df)[1]

# Are there duplicates in the datasheet?
table(table(Spen_Chignik_Datasheet.df$VIAL))

names(which(table(Spen_Chignik_Datasheet.df$VIAL) == 2))



# Subset fish we want for date and area
Data_Area.df <- Spen_Chignik_Datasheet.df[Spen_Chignik_Datasheet.df$VIAL %in% oceanak.df$FK_FISH_ID, c("VIAL", "DATE.HARVESTED", "SAMPLING.AREA")]
str(Data_Area.df)

Data_Area.df$DATE.HARVESTED <- as.character(Data_Area.df$DATE.HARVESTED)
Data_Area.df$SAMPLING.AREA <- as.character(Data_Area.df$SAMPLING.AREA)

# Fix area
unique(Data_Area.df$SAMPLING.AREA)

# Fix date
sum(is.na(Data_Area.df$DATE.HARVESTED))

x <- unique(Data_Area.df$DATE.HARVESTED)
x

newdates <- t(sapply(Data_Area.df$DATE.HARVESTED, function(dat) {
  if(is.na(dat)) {
    rep(NA, 2)
  } else {
    if(grepl(pattern = "40", dat)) {
      rep(format(as.Date(as.numeric(dat), origin = "1904-01-01"), "%m/%d/%Y"), 2)
    } else {
      i.temp = unlist(strsplit(x = dat, split = "-|/"))
      dates <- paste(i.temp[length(i.temp)], i.temp[1], 2014, sep = "-")
      format(c(as.Date(dates, "%b-%d-%Y"), as.Date(dates, "%b-%d-%Y")+1), "%m/%d/%Y")
    }
  }
}))

str(newdates)


Data_Area.df$Start.Date <- newdates[, 1]
Data_Area.df$End.Date <- newdates[, 2]


oceanak.df$CAPTURE_DATE <- Data_Area.df$Start.Date[match(oceanak.df$FK_FISH_ID, Data_Area.df$VIAL)]
oceanak.df$END_CAPTURE_DATE <- Data_Area.df$End.Date[match(oceanak.df$FK_FISH_ID, Data_Area.df$VIAL)]
oceanak.df$CAPTURE_LOCATION <- Data_Area.df$SAMPLING.AREA[match(oceanak.df$FK_FISH_ID, Data_Area.df$VIAL)]


head(oceanak.df)

str(oceanak.df)

# Replace NAs
oceanak.df[is.na(oceanak.df)] <- ""
dimnames(oceanak.df)[[2]][1] <- "FK_COLLECTION_ID"

write.csv(x = oceanak.df, file = "OceanAK Tissue Dump/KCHIGC14/GEN_SAMPLED_FISH_TISSUE Upload.csv", row.names = FALSE, quote = FALSE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### KKODC14 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("V:/Analysis/4_Westward/Chinook/CSRI Westward Commercial Harvest 2014-2016/Mixtures")
rm(list = ls())
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in OceanAK data as data.table (lightning fast!)
require(data.table)
oceanak.dt <- fread(input = "OceanAK Tissue Dump/KKODC14/GEN_SAMPLED_FISH_TISSUE.csv")  # amazing
str(oceanak.dt)
# Convert to data.frame
oceanak.df <- data.frame(oceanak.dt)
str(oceanak.df)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in collection raw datasheet from Birch
require(xlsx)
Kodiak_Datasheet.df <- read.xlsx(file = "Extraction/KMA Chinook Genetics 2014 Extractions.xlsx", sheetName = "2014 KMA Chinook")
str(Kodiak_Datasheet.df)
Spen_Chignik_Datasheet.df <- read.xlsx(file = "South Pen 2014/Extraction/WW Chinook Genetics 2014.xlsx", sheetName = "2014 Pen-Chig Chinook")
str(Spen_Chignik_Datasheet.df)

oceanak.df$FK_FISH_ID[!oceanak.df$FK_FISH_ID %in% Kodiak_Spen_Chignik_Datasheet.df$VIAL] %in% Spen_Chignik_Datasheet.df$VIAL

# Combine SpenChig and Kodiak datasheets as some KKODC14 fish are from Chignik...
dimnames(Kodiak_Datasheet.df)[[2]] %in% dimnames(Spen_Chignik_Datasheet.df)[[2]]
dimnames(Spen_Chignik_Datasheet.df)[[2]] %in% dimnames(Kodiak_Datasheet.df)[[2]]

Kodiak_Spen_Chignik_Datasheet.df <- rbind(Spen_Chignik_Datasheet.df[, dimnames(Spen_Chignik_Datasheet.df)[[2]][-1]],
                                          Kodiak_Datasheet.df[, dimnames(Spen_Chignik_Datasheet.df)[[2]][-1]])
str(Kodiak_Spen_Chignik_Datasheet.df)


# Do the oceanak fish exist in datasheet?
sum(oceanak.df$FK_FISH_ID %in% Kodiak_Spen_Chignik_Datasheet.df$VIAL); dim(oceanak.df)[1]
oceanak.df$FK_FISH_ID[!oceanak.df$FK_FISH_ID %in% Kodiak_Spen_Chignik_Datasheet.df$VIAL]

# Do the datasheet fish exist in the oceanak fish?
table(Kodiak_Spen_Chignik_Datasheet.df$VIAL %in% oceanak.df$FK_FISH_ID); dim(oceanak.df)[1]

# Are there duplicates in the datasheet?
table(table(Kodiak_Spen_Chignik_Datasheet.df$VIAL))

names(which(table(Kodiak_Spen_Chignik_Datasheet.df$VIAL) == 2))






# Subset fish we want for date and area
Data_Area.df <- Kodiak_Spen_Chignik_Datasheet.df[Kodiak_Spen_Chignik_Datasheet.df$VIAL %in% oceanak.df$FK_FISH_ID, c("VIAL", "DATE.HARVESTED", "SAMPLING.AREA")]
str(Data_Area.df)

Data_Area.df$DATE.HARVESTED <- as.character(Data_Area.df$DATE.HARVESTED)
Data_Area.df$SAMPLING.AREA <- as.character(Data_Area.df$SAMPLING.AREA)

# Fix area
unique(Data_Area.df$SAMPLING.AREA)

# Fix date
sum(is.na(Data_Area.df$DATE.HARVESTED))

x <- unique(Data_Area.df$DATE.HARVESTED)
x

# y = x[3]
# i.temp = unlist(strsplit(x = y, split = "-|/|//|, "))
# dates <- paste(i.temp[length(i.temp)], i.temp[1], 2014, sep = "-")
# format(c(as.Date(dates, "%b-%d-%Y"), as.Date(dates, "%b-%d-%Y")+1), "%m/%d/%Y")

y = x[-grep(pattern = "40", x = x)]

sapply(y, function(i) {
  i.temp = unlist(strsplit(x = i, split = "-|/|//|, "))
  if(length(i.temp) == 3){
    dates <- rep(paste(2014, i.temp[1], i.temp[2], sep = "-"), 2)
  } else {
    dates <- c(paste(2014, i.temp[1], i.temp[2], sep = "-"),
               paste(2014, i.temp[1], i.temp[3], sep = "-"))
  }
  format(as.Date(dates), "%m/%d/%Y")
})


newdates <- t(sapply(Data_Area.df$DATE.HARVESTED, function(dat) {
  if(is.na(dat)) {
    rep(NA, 2)
  } else {
    if(grepl(pattern = "40", dat)) {
      rep(format(as.Date(as.numeric(dat), origin = "1904-01-01"), "%m/%d/%Y"), 2)
    } else {
      i.temp = unlist(strsplit(x = dat, split = "-|/|//|, "))
      if(length(i.temp) == 3){
        dates <- rep(paste(2014, i.temp[1], i.temp[2], sep = "-"), 2)
      } else {
        dates <- c(paste(2014, i.temp[1], i.temp[2], sep = "-"),
                   paste(2014, i.temp[1], i.temp[3], sep = "-"))
      }
      format(as.Date(dates), "%m/%d/%Y")
    }
  }
}))

str(newdates)


Data_Area.df$Start.Date <- newdates[, 1]
Data_Area.df$End.Date <- newdates[, 2]


oceanak.df$CAPTURE_DATE <- Data_Area.df$Start.Date[match(oceanak.df$FK_FISH_ID, Data_Area.df$VIAL)]
oceanak.df$END_CAPTURE_DATE <- Data_Area.df$End.Date[match(oceanak.df$FK_FISH_ID, Data_Area.df$VIAL)]
oceanak.df$CAPTURE_LOCATION <- Data_Area.df$SAMPLING.AREA[match(oceanak.df$FK_FISH_ID, Data_Area.df$VIAL)]

sum(is.na(oceanak.df$MESH_SIZE_COMMENT))

head(oceanak.df)
str(oceanak.df)

# Replace NAs
oceanak.df[is.na(oceanak.df)] <- ""
dimnames(oceanak.df)[[2]][1] <- "FK_COLLECTION_ID"
str(oceanak.df)

table(oceanak.df$MESH_SIZE_COMMENT)
oceanak.df$FK_FISH_ID[oceanak.df$MESH_SIZE_COMMENT == "#N/A"]


oceanak.df$MESH_SIZE_COMMENT[oceanak.df$MESH_SIZE_COMMENT == "#N/A"][-c(1:3)] <- ""

oceanak.df$FK_FISH_ID[oceanak.df$MESH_SIZE_COMMENT == "#N/A"]
oceanak.df$MESH_SIZE_COMMENT[oceanak.df$MESH_SIZE_COMMENT == "#N/A"] <- c("Early", "Early", "Late")

str(oceanak.df)


write.csv(x = oceanak.df, file = "OceanAK Tissue Dump/KKODC14/GEN_SAMPLED_FISH_TISSUE Upload.csv", row.names = FALSE, quote = FALSE)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### KALITC15 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("V:/Analysis/4_Westward/Chinook/CSRI Westward Commercial Harvest 2014-2016/Mixtures")
rm(list = ls())
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in OceanAK data as data.table (lightning fast!)
require(data.table)
oceanak.dt <- fread(input = "OceanAK Tissue Dump/KALITC15/GEN_SAMPLED_FISH_TISSUE.csv")  # amazing
str(oceanak.dt)
# Convert to data.frame
oceanak.df <- data.frame(oceanak.dt)
str(oceanak.df)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in collection raw datasheet from Birch
require(xlsx)
Kodiak_Datasheet.df <- read.xlsx(file = "Extraction/OLD/KMA Chinook Sample Selection 2015mbf.xlsx", sheetName = "Data", stringsAsFactors = FALSE)
str(Kodiak_Datasheet.df)



# First need one row per fish to begin with
vials2besplit <- Kodiak_Datasheet.df$vial.number
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
sampsize <- sapply(test, function(row) {length(unlist(row))})
str(sampsize, max.level = 0)
head(sampsize)

# check weird one which was a collection of two different vial ranges
test[352]

# Number of fish sampled for 2015 mixtures
length(unlist(x = test, recursive = TRUE, use.names = FALSE))

colnames(Kodiak_Datasheet.df)


# Is sample size equal to the number of fish thought to be sampled?
table(sampsize == apply(Kodiak_Datasheet.df[, 3:8], 1, function(samp.event) {sum(samp.event, na.rm = TRUE)}))

sampsize[!sampsize == apply(Kodiak_Datasheet.df[, 3:8], 1, function(samp.event) {sum(samp.event, na.rm = TRUE)})]
apply(Kodiak_Datasheet.df[, 3:8], 1, function(samp.event) {sum(samp.event, na.rm = TRUE)})[!apply(Kodiak_Datasheet.df[, 3:8], 1, function(samp.event) {sum(samp.event, na.rm = TRUE)}) == sampsize]

# Disregard the discrepancy with the Kodiak_Datasheet.df. There is a note on vials 212-213 that the 3rd sampled fish had no genetics sample taken

Kodiak_Datasheet.df$SampleSize <- sampsize

str(Kodiak_Datasheet.df)


# Create data.frame with one row per fish
king2015perfish <- apply(X = Kodiak_Datasheet.df, 1, function(row) {matrix(data = rep(row, times = row["SampleSize"]), ncol = 18, byrow = TRUE)} )
king2015perfish.mat <- Reduce(f = rbind, x = king2015perfish)

str(king2015perfish.mat)
king2015perfish.df <- data.frame(king2015perfish.mat, stringsAsFactors = FALSE)
head(king2015perfish.df)
dimnames(king2015perfish.df)[[2]] <- dimnames(Kodiak_Datasheet.df)[[2]]
str(king2015perfish.df)
length(unlist(test, recursive = TRUE, use.names = FALSE))

king2015perfish.df$VialNumber <- unlist(test, recursive = TRUE, use.names = FALSE)

king2015perfish.df.final <- king2015perfish.df[, c("VialNumber", "Sampling.Port", "Area", "single.day", "catch.dates", "Card..", "Strata")]
str(king2015perfish.df.final)

Kodiak_Datasheet_OneRowPerFish.df <- king2015perfish.df.final
str(Kodiak_Datasheet_OneRowPerFish.df)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Do the oceanak fish exist in datasheet?
sum(oceanak.df$FK_FISH_ID %in% Kodiak_Datasheet_OneRowPerFish.df$VialNumber); dim(oceanak.df)[1]
oceanak.df$FK_FISH_ID[!oceanak.df$FK_FISH_ID %in% Kodiak_Datasheet_OneRowPerFish.df$VialNumber]

# Do the datasheet fish exist in the oceanak fish?
table(Kodiak_Datasheet_OneRowPerFish.df$VialNumber %in% oceanak.df$FK_FISH_ID); dim(oceanak.df)[1]

# Are there duplicates in the datasheet?
table(table(Kodiak_Datasheet_OneRowPerFish.df$VialNumber))

names(which(table(Kodiak_Datasheet_OneRowPerFish.df$VialNumber) == 2))






# Subset fish we want for date and area
Data_Area.df <- Kodiak_Datasheet_OneRowPerFish.df[Kodiak_Datasheet_OneRowPerFish.df$VialNumber %in% oceanak.df$FK_FISH_ID, 
                                                  c("VialNumber", "catch.dates", "Area", "Strata")]
colnames(Data_Area.df) <- c("VIAL", "DATE.HARVESTED", "SAMPLING.AREA", "STRATA")
str(Data_Area.df)

Data_Area.df$DATE.HARVESTED <- as.character(Data_Area.df$DATE.HARVESTED)
Data_Area.df$SAMPLING.AREA <- as.character(Data_Area.df$SAMPLING.AREA)

# Fix area
unique(Data_Area.df$SAMPLING.AREA)

# Fix date
sum(is.na(Data_Area.df$DATE.HARVESTED))

x <- unique(Data_Area.df$DATE.HARVESTED)
x

z = x[grep(pattern = "42", x = x)]
sapply(z, function(i) {
  rep(format(as.Date(as.numeric(i), origin = "1899-12-30"), "%m/%d/%Y"), 2)
})

  
y = x[-grep(pattern = "42", x = x)]
y
sapply(y, function(i) {
  i.temp = unlist(strsplit(x = i, split = "-|/|//|, "))
  if(length(i.temp) == 3){
    dates <- rep(paste(2015, i.temp[1], i.temp[2], sep = "-"), 2)
  } else {
    dates <- c(paste(2015, i.temp[1], i.temp[2], sep = "-"),
               paste(2015, i.temp[3], i.temp[4], sep = "-"))
  }
  format(as.Date(dates), "%m/%d/%Y")
})


newdates <- t(sapply(Data_Area.df$DATE.HARVESTED, function(dat) {
  if(is.na(dat)) {
    rep(NA, 2)
  } else {
    if(grepl(pattern = "42", dat)) {
      rep(format(as.Date(as.numeric(dat), origin = "1899-12-30"), "%m/%d/%Y"), 2)
    } else {
      i.temp = unlist(strsplit(x = dat, split = "-|/|//|, "))
      if(length(i.temp) == 3){
        dates <- rep(paste(2015, i.temp[1], i.temp[2], sep = "-"), 2)
      } else {
        dates <- c(paste(2015, i.temp[1], i.temp[2], sep = "-"),
                   paste(2015, i.temp[3], i.temp[4], sep = "-"))
      }
      format(as.Date(dates), "%m/%d/%Y")
    }
  }
}))

str(newdates)


Data_Area.df$Start.Date <- newdates[, 1]
Data_Area.df$End.Date <- newdates[, 2]


oceanak.df$CAPTURE_DATE <- Data_Area.df$Start.Date[match(oceanak.df$FK_FISH_ID, Data_Area.df$VIAL)]
oceanak.df$END_CAPTURE_DATE <- Data_Area.df$End.Date[match(oceanak.df$FK_FISH_ID, Data_Area.df$VIAL)]
oceanak.df$CAPTURE_LOCATION <- Data_Area.df$SAMPLING.AREA[match(oceanak.df$FK_FISH_ID, Data_Area.df$VIAL)]
oceanak.df$MESH_SIZE_COMMENT <- Data_Area.df$STRATA[match(oceanak.df$FK_FISH_ID, Data_Area.df$VIAL)]

sum(is.na(oceanak.df$MESH_SIZE_COMMENT))

head(oceanak.df)
str(oceanak.df)

# Replace NAs
oceanak.df[is.na(oceanak.df)] <- ""
dimnames(oceanak.df)[[2]][1] <- "FK_COLLECTION_ID"
str(oceanak.df)

table(oceanak.df$MESH_SIZE_COMMENT, oceanak.df$CAPTURE_LOCATION)
str(oceanak.df)


write.csv(x = oceanak.df, file = "OceanAK Tissue Dump/KALITC15/GEN_SAMPLED_FISH_TISSUE Upload.csv", row.names = FALSE, quote = FALSE)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### KKODC15 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("V:/Analysis/4_Westward/Chinook/CSRI Westward Commercial Harvest 2014-2016/Mixtures")
rm(list = ls())
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in OceanAK data as data.table (lightning fast!)
require(data.table)
oceanak.dt <- fread(input = "OceanAK Tissue Dump/KKODC15/GEN_SAMPLED_FISH_TISSUE.csv")  # amazing
str(oceanak.dt)
# Convert to data.frame
oceanak.df <- data.frame(oceanak.dt)
str(oceanak.df)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in collection raw datasheet from Birch
require(xlsx)
Kodiak_Datasheet.df <- read.xlsx(file = "Extraction/OLD/KMA Chinook Sample Selection 2015mbf.xlsx", sheetName = "Data", stringsAsFactors = FALSE)
str(Kodiak_Datasheet.df)



# First need one row per fish to begin with
vials2besplit <- Kodiak_Datasheet.df$vial.number
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
sampsize <- sapply(test, function(row) {length(unlist(row))})
str(sampsize, max.level = 0)
head(sampsize)

# check weird one which was a collection of two different vial ranges
test[352]

# Number of fish sampled for 2015 mixtures
length(unlist(x = test, recursive = TRUE, use.names = FALSE))

colnames(Kodiak_Datasheet.df)


# Is sample size equal to the number of fish thought to be sampled?
table(sampsize == apply(Kodiak_Datasheet.df[, 3:8], 1, function(samp.event) {sum(samp.event, na.rm = TRUE)}))

sampsize[!sampsize == apply(Kodiak_Datasheet.df[, 3:8], 1, function(samp.event) {sum(samp.event, na.rm = TRUE)})]
apply(Kodiak_Datasheet.df[, 3:8], 1, function(samp.event) {sum(samp.event, na.rm = TRUE)})[!apply(Kodiak_Datasheet.df[, 3:8], 1, function(samp.event) {sum(samp.event, na.rm = TRUE)}) == sampsize]

# Disregard the discrepancy with the Kodiak_Datasheet.df. There is a note on vials 212-213 that the 3rd sampled fish had no genetics sample taken

Kodiak_Datasheet.df$SampleSize <- sampsize

str(Kodiak_Datasheet.df)


# Create data.frame with one row per fish
king2015perfish <- apply(X = Kodiak_Datasheet.df, 1, function(row) {matrix(data = rep(row, times = row["SampleSize"]), ncol = 18, byrow = TRUE)} )
king2015perfish.mat <- Reduce(f = rbind, x = king2015perfish)

str(king2015perfish.mat)
king2015perfish.df <- data.frame(king2015perfish.mat, stringsAsFactors = FALSE)
head(king2015perfish.df)
dimnames(king2015perfish.df)[[2]] <- dimnames(Kodiak_Datasheet.df)[[2]]
str(king2015perfish.df)
length(unlist(test, recursive = TRUE, use.names = FALSE))

king2015perfish.df$VialNumber <- unlist(test, recursive = TRUE, use.names = FALSE)

king2015perfish.df.final <- king2015perfish.df[, c("VialNumber", "Sampling.Port", "Area", "single.day", "catch.dates", "Card..", "Strata")]
str(king2015perfish.df.final)

Kodiak_Datasheet_OneRowPerFish.df <- king2015perfish.df.final
str(Kodiak_Datasheet_OneRowPerFish.df)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Do the oceanak fish exist in datasheet?
sum(oceanak.df$FK_FISH_ID %in% Kodiak_Datasheet_OneRowPerFish.df$VialNumber); dim(oceanak.df)[1]
oceanak.df$FK_FISH_ID[!oceanak.df$FK_FISH_ID %in% Kodiak_Datasheet_OneRowPerFish.df$VialNumber]

# Do the datasheet fish exist in the oceanak fish?
table(Kodiak_Datasheet_OneRowPerFish.df$VialNumber %in% oceanak.df$FK_FISH_ID); dim(oceanak.df)[1]

# Are there duplicates in the datasheet?
table(table(Kodiak_Datasheet_OneRowPerFish.df$VialNumber))

names(which(table(Kodiak_Datasheet_OneRowPerFish.df$VialNumber) == 2))


## Data sheet indicates that vials 1152-1157 were not used
str(Kodiak_Datasheet_OneRowPerFish.df)

new.rows.df <- data.frame(VialNumber = 1152:1157,
                          Sampling.Port = rep('', 6),
                          Area = rep('', 6),
                          single.day = rep('', 6),
                          catch.dates = rep('', 6),
                          Card.. = rep('', 6),
                          Strata = rep('', 6), stringsAsFactors = FALSE)
str(new.rows.df)
Kodiak_Datasheet_OneRowPerFish_new.df <- rbind(Kodiak_Datasheet_OneRowPerFish.df, new.rows.df)
sum(oceanak.df$FK_FISH_ID %in% Kodiak_Datasheet_OneRowPerFish_new.df$VialNumber); dim(oceanak.df)[1]
tail(Kodiak_Datasheet_OneRowPerFish_new.df, 10)



# Subset fish we want for date and area
Data_Area.df <- Kodiak_Datasheet_OneRowPerFish_new.df[Kodiak_Datasheet_OneRowPerFish_new.df$VialNumber %in% oceanak.df$FK_FISH_ID, 
                                                  c("VialNumber", "catch.dates", "Area", "Strata")]
colnames(Data_Area.df) <- c("VIAL", "DATE.HARVESTED", "SAMPLING.AREA", "STRATA")
str(Data_Area.df)

Data_Area.df$DATE.HARVESTED <- as.character(Data_Area.df$DATE.HARVESTED)
Data_Area.df$SAMPLING.AREA <- as.character(Data_Area.df$SAMPLING.AREA)

# Fix area
unique(Data_Area.df$SAMPLING.AREA)
table((Data_Area.df$SAMPLING.AREA), useNA = "always")


# Fix date
sum(is.na(Data_Area.df$DATE.HARVESTED))

x <- unique(Data_Area.df$DATE.HARVESTED)
x

z = x[grep(pattern = "42", x = x)]
sapply(z, function(i) {
  rep(format(as.Date(as.numeric(i), origin = "1899-12-30"), "%m/%d/%Y"), 2)
})


y = x[-grep(pattern = "42", x = x)]
y <- y[!y == '']
y
sapply(y, function(i) {
  i.temp = unlist(strsplit(x = i, split = "-|/|//|, "))
  if(length(i.temp) == 3){
    dates <- rep(paste(2015, i.temp[1], i.temp[2], sep = "-"), 2)
  } else {
    dates <- c(paste(2015, i.temp[1], i.temp[2], sep = "-"),
               paste(2015, i.temp[3], i.temp[4], sep = "-"))
  }
  format(as.Date(dates), "%m/%d/%Y")
})


newdates <- t(sapply(Data_Area.df$DATE.HARVESTED, function(dat) {
  if(dat == "") {
    rep(dat, 2)
  } else {
    if(grepl(pattern = "42", dat)) {
      rep(format(as.Date(as.numeric(dat), origin = "1899-12-30"), "%m/%d/%Y"), 2)
    } else {
      i.temp = unlist(strsplit(x = dat, split = "-|/|//|, "))
      if(length(i.temp) == 3){
        dates <- rep(paste(2015, i.temp[1], i.temp[2], sep = "-"), 2)
      } else {
        dates <- c(paste(2015, i.temp[1], i.temp[2], sep = "-"),
                   paste(2015, i.temp[3], i.temp[4], sep = "-"))
      }
      format(as.Date(dates), "%m/%d/%Y")
    }
  }
}))

str(newdates)
head(newdates, 10); tail(newdates, 10)


Data_Area.df$Start.Date <- newdates[, 1]
Data_Area.df$End.Date <- newdates[, 2]


oceanak.df$CAPTURE_DATE <- Data_Area.df$Start.Date[match(oceanak.df$FK_FISH_ID, Data_Area.df$VIAL)]
oceanak.df$END_CAPTURE_DATE <- Data_Area.df$End.Date[match(oceanak.df$FK_FISH_ID, Data_Area.df$VIAL)]
oceanak.df$CAPTURE_LOCATION <- Data_Area.df$SAMPLING.AREA[match(oceanak.df$FK_FISH_ID, Data_Area.df$VIAL)]
oceanak.df$MESH_SIZE_COMMENT <- Data_Area.df$STRATA[match(oceanak.df$FK_FISH_ID, Data_Area.df$VIAL)]

sum(is.na(oceanak.df$MESH_SIZE_COMMENT))

head(oceanak.df)
str(oceanak.df)

# Replace NAs
oceanak.df[is.na(oceanak.df)] <- ""
dimnames(oceanak.df)[[2]][1] <- "FK_COLLECTION_ID"
str(oceanak.df)

# Fix Mesh_size_comment
oceanak.df$MESH_SIZE_COMMENT <- gsub(pattern = "E ", replacement = "E", x = oceanak.df$MESH_SIZE_COMMENT)

table(oceanak.df$MESH_SIZE_COMMENT, oceanak.df$CAPTURE_LOCATION)
str(oceanak.df)


write.csv(x = oceanak.df, file = "OceanAK Tissue Dump/KKODC15/GEN_SAMPLED_FISH_TISSUE Upload.csv", row.names = FALSE, quote = FALSE)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### KLARSC15 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("V:/Analysis/4_Westward/Chinook/CSRI Westward Commercial Harvest 2014-2016/Mixtures")
rm(list = ls())
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in OceanAK data as data.table (lightning fast!)
require(data.table)
oceanak.dt <- fread(input = "OceanAK Tissue Dump/KLARSC15/GEN_SAMPLED_FISH_TISSUE.csv")  # amazing
str(oceanak.dt)
# Convert to data.frame
oceanak.df <- data.frame(oceanak.dt)
str(oceanak.df)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in collection raw datasheet from Birch
require(xlsx)
Kodiak_Datasheet.df <- read.xlsx(file = "Extraction/OLD/KMA Chinook Sample Selection 2015mbf.xlsx", sheetName = "Data", stringsAsFactors = FALSE)
str(Kodiak_Datasheet.df)



# First need one row per fish to begin with
vials2besplit <- Kodiak_Datasheet.df$vial.number
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
sampsize <- sapply(test, function(row) {length(unlist(row))})
str(sampsize, max.level = 0)
head(sampsize)

# check weird one which was a collection of two different vial ranges
test[352]

# Number of fish sampled for 2015 mixtures
length(unlist(x = test, recursive = TRUE, use.names = FALSE))

colnames(Kodiak_Datasheet.df)


# Is sample size equal to the number of fish thought to be sampled?
table(sampsize == apply(Kodiak_Datasheet.df[, 3:8], 1, function(samp.event) {sum(samp.event, na.rm = TRUE)}))

sampsize[!sampsize == apply(Kodiak_Datasheet.df[, 3:8], 1, function(samp.event) {sum(samp.event, na.rm = TRUE)})]
apply(Kodiak_Datasheet.df[, 3:8], 1, function(samp.event) {sum(samp.event, na.rm = TRUE)})[!apply(Kodiak_Datasheet.df[, 3:8], 1, function(samp.event) {sum(samp.event, na.rm = TRUE)}) == sampsize]

# Disregard the discrepancy with the Kodiak_Datasheet.df. There is a note on vials 212-213 that the 3rd sampled fish had no genetics sample taken

Kodiak_Datasheet.df$SampleSize <- sampsize

str(Kodiak_Datasheet.df)


# Create data.frame with one row per fish
king2015perfish <- apply(X = Kodiak_Datasheet.df, 1, function(row) {matrix(data = rep(row, times = row["SampleSize"]), ncol = 18, byrow = TRUE)} )
king2015perfish.mat <- Reduce(f = rbind, x = king2015perfish)

str(king2015perfish.mat)
king2015perfish.df <- data.frame(king2015perfish.mat, stringsAsFactors = FALSE)
head(king2015perfish.df)
dimnames(king2015perfish.df)[[2]] <- dimnames(Kodiak_Datasheet.df)[[2]]
str(king2015perfish.df)
length(unlist(test, recursive = TRUE, use.names = FALSE))

king2015perfish.df$VialNumber <- unlist(test, recursive = TRUE, use.names = FALSE)

king2015perfish.df.final <- king2015perfish.df[, c("VialNumber", "Sampling.Port", "Area", "single.day", "catch.dates", "Card..", "Strata")]
str(king2015perfish.df.final)

Kodiak_Datasheet_OneRowPerFish.df <- king2015perfish.df.final
str(Kodiak_Datasheet_OneRowPerFish.df)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Do the oceanak fish exist in datasheet?
sum(oceanak.df$FK_FISH_ID %in% Kodiak_Datasheet_OneRowPerFish.df$VialNumber); dim(oceanak.df)[1]
oceanak.df$FK_FISH_ID[!oceanak.df$FK_FISH_ID %in% Kodiak_Datasheet_OneRowPerFish.df$VialNumber]

# Do the datasheet fish exist in the oceanak fish?
table(Kodiak_Datasheet_OneRowPerFish.df$VialNumber %in% oceanak.df$FK_FISH_ID); dim(oceanak.df)[1]

# Are there duplicates in the datasheet?
table(table(Kodiak_Datasheet_OneRowPerFish.df$VialNumber))

names(which(table(Kodiak_Datasheet_OneRowPerFish.df$VialNumber) == 2))



# Subset fish we want for date and area
Data_Area.df <- Kodiak_Datasheet_OneRowPerFish.df[Kodiak_Datasheet_OneRowPerFish.df$VialNumber %in% oceanak.df$FK_FISH_ID, 
                                                      c("VialNumber", "catch.dates", "Area", "Strata")]
colnames(Data_Area.df) <- c("VIAL", "DATE.HARVESTED", "SAMPLING.AREA", "STRATA")
str(Data_Area.df)

Data_Area.df$DATE.HARVESTED <- as.character(Data_Area.df$DATE.HARVESTED)
Data_Area.df$SAMPLING.AREA <- as.character(Data_Area.df$SAMPLING.AREA)

# Fix area
unique(Data_Area.df$SAMPLING.AREA)
table((Data_Area.df$SAMPLING.AREA), useNA = "always")


# Fix date
sum(is.na(Data_Area.df$DATE.HARVESTED))

x <- unique(Data_Area.df$DATE.HARVESTED)
x

z = x[grep(pattern = "42", x = x)]
sapply(z, function(i) {
  rep(format(as.Date(as.numeric(i), origin = "1899-12-30"), "%m/%d/%Y"), 2)
})


y = x[-grep(pattern = "42", x = x)]
y
sapply(y, function(i) {
  i.temp = unlist(strsplit(x = i, split = "-|/|//|, "))
  if(length(i.temp) == 3){
    dates <- rep(paste(2015, i.temp[1], i.temp[2], sep = "-"), 2)
  } else {
    if(as.numeric(i.temp[1]) == 6) {
      dates <- c(paste(2015, i.temp[1], i.temp[2], sep = "-"),
                 paste(2015, i.temp[1], i.temp[3], sep = "-"))  
      
    } else {
      dates <- c(paste(2015, i.temp[1], i.temp[2], sep = "-"),
                 paste(2015, i.temp[3], i.temp[4], sep = "-"))  
    }
  }
  format(as.Date(dates), "%m/%d/%Y")
})


newdates <- t(sapply(Data_Area.df$DATE.HARVESTED, function(dat) {
  if(dat == "") {
    rep(dat, 2)
  } else {
    if(grepl(pattern = "42", dat)) {
      rep(format(as.Date(as.numeric(dat), origin = "1899-12-30"), "%m/%d/%Y"), 2)
    } else {
      i.temp = unlist(strsplit(x = dat, split = "-|/|//|, "))
      if(length(i.temp) == 3){
        dates <- rep(paste(2015, i.temp[1], i.temp[2], sep = "-"), 2)
      } else {
        if(as.numeric(i.temp[1]) == 6) {
          dates <- c(paste(2015, i.temp[1], i.temp[2], sep = "-"),
                     paste(2015, i.temp[1], i.temp[3], sep = "-"))  
          
        } else {
          dates <- c(paste(2015, i.temp[1], i.temp[2], sep = "-"),
                     paste(2015, i.temp[3], i.temp[4], sep = "-"))  
        }
      }
      format(as.Date(dates), "%m/%d/%Y")
    }
  }
}))

str(newdates)
head(newdates, 10); tail(newdates, 10)


Data_Area.df$Start.Date <- newdates[, 1]
Data_Area.df$End.Date <- newdates[, 2]


oceanak.df$CAPTURE_DATE <- Data_Area.df$Start.Date[match(oceanak.df$FK_FISH_ID, Data_Area.df$VIAL)]
oceanak.df$END_CAPTURE_DATE <- Data_Area.df$End.Date[match(oceanak.df$FK_FISH_ID, Data_Area.df$VIAL)]
oceanak.df$CAPTURE_LOCATION <- Data_Area.df$SAMPLING.AREA[match(oceanak.df$FK_FISH_ID, Data_Area.df$VIAL)]
oceanak.df$MESH_SIZE_COMMENT <- Data_Area.df$STRATA[match(oceanak.df$FK_FISH_ID, Data_Area.df$VIAL)]

sum(is.na(oceanak.df$MESH_SIZE_COMMENT))

head(oceanak.df)
str(oceanak.df)

# Replace NAs
oceanak.df[is.na(oceanak.df)] <- ""
dimnames(oceanak.df)[[2]][1] <- "FK_COLLECTION_ID"
str(oceanak.df)

# Fix Mesh_size_comment
table(oceanak.df$MESH_SIZE_COMMENT, oceanak.df$CAPTURE_LOCATION)
str(oceanak.df)


write.csv(x = oceanak.df, file = "OceanAK Tissue Dump/KLARSC15/GEN_SAMPLED_FISH_TISSUE Upload.csv", row.names = FALSE, quote = FALSE)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### KKMAC16 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("V:/Analysis/4_Westward/Chinook/CSRI Westward Commercial Harvest 2014-2016/Mixtures")
rm(list = ls())
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in OceanAK data as data.table (lightning fast!)
require(data.table)
oceanak.dt <- fread(input = "OceanAK Tissue Dump/KKMAC16/GEN_SAMPLED_FISH_TISSUE.csv")  # amazing
str(oceanak.dt)
# Convert to data.frame
oceanak.df <- data.frame(oceanak.dt)
str(oceanak.df)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in collection raw datasheet from Birch
require(xlsx)
Kodiak_Datasheet.df <- read.xlsx(file = "Extraction/2016 Chinook Samples for Genetics Lab.xlsx", sheetName = "Chinook Extraction by Fish", stringsAsFactors = FALSE)
str(Kodiak_Datasheet.df)


table(is.na(oceanak.df$FK_FISH_ID))
table(is.na(Kodiak_Datasheet.df$dna_vial))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Do the oceanak fish exist in datasheet?
sum(oceanak.df$FK_FISH_ID %in% Kodiak_Datasheet.df$dna_vial); dim(oceanak.df)[1]
oceanak.df$FK_FISH_ID[!oceanak.df$FK_FISH_ID %in% Kodiak_Datasheet.df$dna_vial]

# Do the datasheet fish exist in the oceanak fish?
table(Kodiak_Datasheet.df$dna_vial %in% oceanak.df$FK_FISH_ID); dim(oceanak.df)[1]

# Are there duplicates in the datasheet?
table(table(Kodiak_Datasheet.df$dna_vial))

names(which(table(Kodiak_Datasheet.df$dna_vial) == 2))




oceanak.df$CAPTURE_DATE <- Kodiak_Datasheet.df$Start.Date[match(oceanak.df$FK_FISH_ID, Kodiak_Datasheet.df$dna_vial)]
oceanak.df$END_CAPTURE_DATE <- Kodiak_Datasheet.df$End.Date[match(oceanak.df$FK_FISH_ID, Kodiak_Datasheet.df$dna_vial)]
oceanak.df$CAPTURE_LOCATION <- Kodiak_Datasheet.df$area_sampled[match(oceanak.df$FK_FISH_ID, Kodiak_Datasheet.df$dna_vial)]
oceanak.df$MESH_SIZE_COMMENT <- Kodiak_Datasheet.df$Stratra[match(oceanak.df$FK_FISH_ID, Kodiak_Datasheet.df$dna_vial)]
str(oceanak.df)


# Fix Capture Location
unique(oceanak.df$CAPTURE_LOCATION)
oceanak.df$CAPTURE_LOCATION[oceanak.df$CAPTURE_LOCATION == unique(oceanak.df$CAPTURE_LOCATION)[5]] <- "Mixed"

# Fix Capture Location
unique(oceanak.df$MESH_SIZE_COMMENT)
oceanak.df$MESH_SIZE_COMMENT[oceanak.df$MESH_SIZE_COMMENT == "Early"] <- "E"
oceanak.df$MESH_SIZE_COMMENT[oceanak.df$MESH_SIZE_COMMENT == "Late"] <- "L"



sum(is.na(oceanak.df$MESH_SIZE_COMMENT))

head(oceanak.df)
str(oceanak.df)

# Replace NAs
oceanak.df[is.na(oceanak.df)] <- ""
dimnames(oceanak.df)[[2]][1] <- "FK_COLLECTION_ID"
str(oceanak.df)

# Fix Mesh_size_comment
table(oceanak.df$MESH_SIZE_COMMENT, oceanak.df$CAPTURE_LOCATION)
str(oceanak.df)

# Reformat dates
oceanak.df$CAPTURE_DATE <- format(oceanak.df$CAPTURE_DATE, format = "%m/%d/%Y")
oceanak.df$END_CAPTURE_DATE <- format(oceanak.df$END_CAPTURE_DATE, format = "%m/%d/%Y")
head(oceanak.df)


write.csv(x = oceanak.df, file = "OceanAK Tissue Dump/KKMAC16/GEN_SAMPLED_FISH_TISSUE Upload.csv", row.names = FALSE, quote = FALSE)
# Done