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
### SPENC14
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
### SPENC14
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

