#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### KMA Chinook Mixtures 2014-2016 ####
# Kyle Shedd Mon Jul 11 15:39:50 2016
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
date()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Introduction ####
# The goal of this script is to analyze Chinook salmon mixtures from the KMA
# commercial harvest from 2014-2016 using a coastwide Chinook baseline
# containing 211 populations in 10 reporting groups characterized by 48 SNPs
# grouped into 46 loci (2 sets of linked SNPs are combined as haplotypes).
# All mixtures are to be analyzed with the program BAYES.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Specific Objectives ####
# This script will:
# 1) Import mixture data
# 2) Define spatio-temporal strata
# 3) Perform a data QC on mixtures
# 4) Prepare BAYES input files
# 5) Summarize BAYES results
# 6) Generate plots and tables of results


## Reporting groups
# 1. Russia
# 2. Coastal West Alaska/Yukon
# 3. Alaska Peninsula
# 4. Chignik
# 5. Kodiak
# 6. Cook Inlet
# 7. Copper
# 8. SE Alaska/NE Gulf of Alaska
# 9. British Columbia
# 10. West Coast U.S.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Initial Setup ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls(all = TRUE))
setwd("V:/Analysis/4_Westward/Chinook/CSRI Westward Commercial Harvest 2014-2016/Mixtures")
source("H:/R Source Scripts/Functions.GCL_KS.R")
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")
username <- "krshedd"
password <- "********"


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Pull genotypes from LOKI 2014/2015 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Get collection SILLYs
dir.create("Objects")

SPEN2014 <- c("KSPENC14", "KCHIGC14")
dput(SPEN2014, file = "Objects/SPEN2014.txt")

KMA2014 <- c("KKODC14")
dput(x = KMA2014, file = "Objects/KMA2014.txt")
KMA2015 <- c("KALITC15", "KKODC15", "KLARSC15")
dput(x = KMA2015, file = "Objects/KMA2015.txt")
# KMA2016 <- c("KALITC16", "KKODC16", "KLARSC16")
# dput(x = KMA2016, file = "Objects/KMA2016.txt")

## Create Locus Control
CreateLocusControl.GCL(markersuite = "Chinook_Kodiak_2016_48SNPs", username = username, password = password)

## Save original LocusControl
loci48 <- LocusControl$locusnames
mito.loci48 <- which(LocusControl$ploidy == 1)

dput(x = LocusControl, file = "Objects/OriginalLocusControl48.txt")
dput(x = loci48, file = "Objects/loci48.txt")
dput(x = mito.loci48, file = "Objects/mito.loci48.txt")

#~~~~~~~~~~~~~~~~~~
## Pull all data for each silly code and create .gcl objects for each
LOKI2R.GCL(sillyvec = c(KMA2014, KMA2015, SPEN2014), username = username, password = password)

rm(username, password)
objects(pattern = "\\.gcl")

## Save unaltered .gcl's as back-up:
dir.create("Raw genotypes")
dir.create("Raw genotypes/OriginalCollections")
invisible(sapply(c(KMA2014, KMA2015, SPEN2014), function(silly) {dput(x = get(paste(silly, ".gcl", sep = '')), file = paste("Raw genotypes/OriginalCollections/" , silly, ".txt", sep = ''))} )); beep(8)

## Original sample sizes by SILLY
collection.size.original <- sapply(c(KMA2014, KMA2015, SPEN2014), function(silly) get(paste(silly, ".gcl", sep = ""))$n)
# matrix(data = collection.size.original, ncol = 2, dimnames = list(c("Alitak", "Ayakulik", "Karluk", "Uganik", "Uyak"), 2014:2015))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Clean workspace; dget .gcl objects and Locus Control ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls(all = TRUE))
setwd("V:/Analysis/4_Westward/Chinook/CSRI Westward Commercial Harvest 2014-2016/Mixtures")
# This sources all of the new GCL functions to this workspace
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")
source("H:/R Source Scripts/Functions.GCL_KS.R")

## Get objects
LocusControl <- dget(file = "Objects/OriginalLocusControl48.txt")

KMAobjects <- list.files(path = "Objects", recursive = FALSE)
KMAobjects <- KMAobjects[!KMAobjects %in% c("OriginalLocusControl48.txt")]
KMAobjects

invisible(sapply(KMAobjects, function(objct) {assign(x = unlist(strsplit(x = objct, split = ".txt")), value = dget(file = paste(getwd(), "Objects", objct, sep = "/")), pos = 1) })); beep(2)


## Get un-altered mixtures
invisible(sapply(c(SPEN2014, KMA2014, KMA2015), function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/OriginalCollections/", silly, ".txt", sep = "")), pos = 1)} )); beep(2)
objects(pattern = "\\.gcl")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Define strata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## All fish have a capture date?
sapply(c(SPEN2014, KMA2014, KMA2015), function(silly) {sum(is.na(get(paste(silly, ".gcl", sep = ''))$attributes$CAPTURE_DATE))} )  # Zeros are good

## Confirming samples sizes by date
sapply(c(SPEN2014, KMA2014, KMA2015), function(silly) {table(get(paste(silly, ".gcl", sep = ''))$attributes$CAPTURE_DATE)} )

str(KCHIGC14.gcl)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Pooling roadmap
# 1) Pool all 2014 and all 2015 samples into two sillys by year
# 2) Pool geographically within year by "CAPTURE_LOCATION"
# 3) Pool temporally within geographic by year by "CAPTURE_DATE"

#~~~~~~~~~~~~~~~~~~
## South Pen + Chignik 2014
KSPENCHIG14.gcl$n; KCHIGC14.gcl$n
table(KKODC14.gcl$attributes$CAPTURE_LOCATION)

# Grab Chignik fish out of KKODC14
KKODC14_KCHIG14IDs <- AttributesToIDs.GCL(silly = "KKODC14", attribute = "CAPTURE_LOCATION", matching = "Chignik Outside")
KKODC14_KCHIG14IDs <- list(na.omit(KKODC14_KCHIG14IDs))
names(KKODC14_KCHIG14IDs) <- "KKODC14"

# Pool Chignik and South Pen fish
PoolCollections.GCL(collections = "KKODC14", loci = loci48, IDs = KKODC14_KCHIG14IDs, newname = "KKODC14_KCHIG14")
str(KKODC14_KCHIG14.gcl)

PoolCollections.GCL(collections = c(SPEN2014, "KKODC14_KCHIG14"), loci = loci48, IDs = NULL, newname = "KSPENCHIG14")
str(KSPENCHIG14.gcl)
table(KSPENCHIG14.gcl$attributes$CAPTURE_LOCATION)


#~~~~~~~~~~~~~~~~~~
## Kodiak 2014, remove Chignik fish
str(KKODC14.gcl)
table(KKODC14.gcl$attributes$CAPTURE_LOCATION)
unique(KKODC14.gcl$attributes$CAPTURE_LOCATION)[1:4]

# Grab non-Chignik fish out of KKODC14
KKODC14IDs <- AttributesToIDs.GCL(silly = "KKODC14", attribute = "CAPTURE_LOCATION", matching = unique(KKODC14.gcl$attributes$CAPTURE_LOCATION)[1:4])
KKODC14IDs <- list(na.omit(KKODC14IDs))
names(KKODC14IDs) <- "KKODC14"

# Pool Chignik and South Pen fish
PoolCollections.GCL(collections = "KKODC14", loci = loci48, IDs = KKODC14IDs, newname = "KMA2014")
str(KMA2014.gcl)

table(KMA2014.gcl$attributes$CAPTURE_LOCATION, KMA2014.gcl$attributes$MESH_SIZE_COMMENT)

#~~~~~~~~~~~~~~~~~~
## Kodiak 2015
KMA2015
str(KKODC15.gcl)
sapply(KMA2015, function(silly) {table(get(paste(silly, ".gcl", sep = ''))$attributes$CAPTURE_LOCATION)})

# Pool Chignik and South Pen fish
PoolCollections.GCL(collections = KMA2015, loci = loci48, IDs = NULL, newname = "KMA2015")
str(KMA2015.gcl)

table(KMA2015.gcl$attributes$CAPTURE_LOCATION, KMA2015.gcl$attributes$MESH_SIZE_COMMENT)


KMA2015.gcl$attributes$MESH_SIZE_COMMENT <- gsub(pattern = "E", replacement = "Early", x = KMA2015.gcl$attributes$MESH_SIZE_COMMENT)
KMA2015.gcl$attributes$MESH_SIZE_COMMENT <- gsub(pattern = "L", replacement = "Late", x = KMA2015.gcl$attributes$MESH_SIZE_COMMENT)


#~~~~~~~~~~~~~~~~~~
# Dput these sillys for now
dir.create("Raw genotypes/PoolingByYear")
invisible(sapply(c("KMA2014", "KMA2015", "KSPENCHIG14"), function(silly) {dput(x = get(paste(silly, ".gcl", sep = '')), file = paste("Raw genotypes/PoolingByYear/" , silly, ".txt", sep = ''))} )); beep(8)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Spatio temporal strata for 2014 Kodiak samples

silly = "KMA2014"
loci = loci48
geostrata = cbind(object = unique(KMA2014.gcl$attributes$CAPTURE_LOCATION), nm = c("KALITC14", "KEASTC14", "KMAINC14", "KWESTC14"))
tempstrata = cbind(object = unique(KMA2014.gcl$attributes$MESH_SIZE_COMMENT), nm = c("1_Early", "2_Late"))
geostrata; tempstrata

table(KMA2014.gcl$attributes$CAPTURE_LOCATION, KMA2014.gcl$attributes$MESH_SIZE_COMMENT)

my.gcl <- get(paste(silly, ".gcl", sep = ''))
for(i in seq(nrow(geostrata))) {
  geo <- geostrata[i, ]
  for(j in seq(nrow(tempstrata))) {
    temp <- tempstrata[j, ]
    
    if(geo["object"] == "Mainland" & temp["object"] == "Early") {next}
    
    geoIDs <- AttributesToIDs.GCL(silly = silly, attribute = "CAPTURE_LOCATION", matching = geo["object"])
    tempIDs <- AttributesToIDs.GCL(silly = silly, attribute = "MESH_SIZE_COMMENT", matching = temp["object"])
    IDs <- geoIDs[geoIDs %in% tempIDs]
    IDs <- list(na.omit(IDs))
    names(IDs) <- silly
    PoolCollections.GCL(collections = silly, loci = loci, IDs = IDs, newname = paste(geo["nm"], temp["nm"], sep = "_"))
    print(get(paste(geo["nm"], "_", temp["nm"], ".gcl", sep = ''))$n)
  }
}

grep(pattern = ".gcl", x = objects(pattern = "14_"), value = TRUE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Spatio temporal strata for 2015 Kodiak samples
silly = "KMA2015"
loci = loci48
geostrata = cbind(object = unique(KMA2015.gcl$attributes$CAPTURE_LOCATION), nm = c("KALITC15", "KEASTC15", "KWESTC15", "KMAINC15"))
tempstrata = cbind(object = unique(KMA2015.gcl$attributes$MESH_SIZE_COMMENT), nm = c("2_Late", "1_Early"))
geostrata; tempstrata

table(KMA2015.gcl$attributes$CAPTURE_LOCATION, KMA2015.gcl$attributes$MESH_SIZE_COMMENT)

my.gcl <- get(paste(silly, ".gcl", sep = ''))
for(i in seq(nrow(geostrata))) {
  geo <- geostrata[i, ]
  for(j in seq(nrow(tempstrata))) {
    temp <- tempstrata[j, ]
    
    if(geo["object"] == "Mainland" & temp["object"] == "Early") {next}
    
    geoIDs <- AttributesToIDs.GCL(silly = silly, attribute = "CAPTURE_LOCATION", matching = geo["object"])
    tempIDs <- AttributesToIDs.GCL(silly = silly, attribute = "MESH_SIZE_COMMENT", matching = temp["object"])
    IDs <- geoIDs[geoIDs %in% tempIDs]
    IDs <- list(na.omit(IDs))
    names(IDs) <- silly
    PoolCollections.GCL(collections = silly, loci = loci, IDs = IDs, newname = paste(geo["nm"], temp["nm"], sep = "_"))
    print(list("First Last Fish" = range(as.numeric(unlist(IDs))), "n" = get(paste(geo["nm"], "_", temp["nm"], ".gcl", sep = ''))$n))
  }
}

grep(pattern = ".gcl", x = objects(pattern = "15_"), value = TRUE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Mixture sillyvec
KMA2014Strata <- unlist(strsplit(x = grep(pattern = "14_", x = objects(pattern = "\\.gcl"), value = TRUE), split = "\\.gcl"))
KMA2014Strata <- KMA2014Strata[-5]
KMA2015Strata <- unlist(strsplit(x = grep(pattern = "15_", x = objects(pattern = "\\.gcl"), value = TRUE), split = "\\.gcl"))
KMA2014_2015Strata <- c("KSPENCHIG14", KMA2014Strata, KMA2015Strata)

dput(x = KMA2014Strata, file = "Objects/KMA2014Strata.txt")
dput(x = KMA2015Strata, file = "Objects/KMA2015Strata.txt")
dput(x = KMA2014_2015Strata, file = "Objects/KMA2014_2015Strata.txt")

# Confirm sample sizes
sapply(KMA2014_2015Strata, function(silly) get(paste(silly, ".gcl", sep = ""))$n)

# View as tables by year
require(reshape)
samp.df.2014 <- data.frame(t(sapply(KMA2014Strata, function(strata) {
  location <- unlist(strsplit(x = strata, split = "_"))[1]
  temporal.strata <- as.numeric(unlist(strsplit(x = strata, split = "_"))[2])
  n <- get(paste(strata, ".gcl", sep = ""))$n
  c(location = location, temporal.strata = temporal.strata, n = n)
})))
cast(data = samp.df.2014, location ~ temporal.strata)

samp.df.2015 <- data.frame(t(sapply(KMA2015Strata, function(strata) {
  location <- unlist(strsplit(x = strata, split = "_"))[1]
  temporal.strata <- as.numeric(unlist(strsplit(x = strata, split = "_"))[2])
  n <- get(paste(strata, ".gcl", sep = ""))$n
  c(location = location, temporal.strata = temporal.strata, n = n)
})))
cast(data = samp.df.2015, location ~ temporal.strata)



# dput mixture sillys
dir.create("Raw genotypes/OriginalCollections_Strata")
invisible(sapply(KMA2014_2015Strata, function(silly) {dput(x = get(paste(silly, ".gcl", sep = '')), file = paste("Raw genotypes/OriginalCollections_Strata/" , silly, ".txt", sep = ''))} )); beep(8)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Clean workspace; dget .gcl objects and Locus Control ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls(all = TRUE))
setwd("V:/Analysis/4_Westward/Chinook/CSRI Westward Commercial Harvest 2014-2016/Mixtures")
# This sources all of the new GCL functions to this workspace
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")
source("H:/R Source Scripts/Functions.GCL_KS.R")

## Get objects
LocusControl <- dget(file = "Objects/OriginalLocusControl48.txt")

KMAobjects <- list.files(path = "Objects", recursive = FALSE)
KMAobjects <- KMAobjects[!KMAobjects %in% c("OriginalLocusControl48.txt")]
KMAobjects

invisible(sapply(KMAobjects, function(objct) {assign(x = unlist(strsplit(x = objct, split = ".txt")), value = dget(file = paste(getwd(), "Objects", objct, sep = "/")), pos = 1) })); beep(2)


## Get un-altered mixtures
invisible(sapply(KMA2014_2015Strata, function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/OriginalCollections_Strata/", silly, ".txt", sep = "")), pos = 1)} )); beep(2)
objects(pattern = "\\.gcl")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Data QC/Massage ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Confirm appropriate strata
sapply(KMA2014_2015Strata, function(mix) {table(get(paste(mix, ".gcl", sep = ''))$attributes$CAPTURE_LOCATION,
                                                get(paste(mix, ".gcl", sep = ''))$attributes$MESH_SIZE_COMMENT)})
str(KSPENCHIG14.gcl)

# Begin QC
require(xlsx)

KMA2014_2015Strata

KMA2014_2015Strata_SampleSizes <- matrix(data = NA, nrow = length(KMA2014_2015Strata), ncol = 4, 
                                         dimnames = list(KMA2014_2015Strata, c("Genotyped", "Missing", "Duplicate", "Final")))

#### Check loci
## Get sample size by locus
Original_KMA2014_2015Strata_SampleSizebyLocus <- SampSizeByLocus.GCL(sillyvec = KMA2014_2015Strata, loci = loci48)
min(Original_KMA2014_2015Strata_SampleSizebyLocus)  ## 264
apply(Original_KMA2014_2015Strata_SampleSizebyLocus, 1, min) / apply(Original_KMA2014_2015Strata_SampleSizebyLocus, 1, max)  # Good, all > 0.9

Original_KMA2014_2015Strata_PercentbyLocus <- apply(Original_KMA2014_2015Strata_SampleSizebyLocus, 1, function(row) {row / max(row)} )
which(apply(Original_KMA2014_2015Strata_PercentbyLocus, 2, min) < 0.8)  # no re-runs!

require(lattice)
new.colors <- colorRampPalette(c("black", "white"))
levelplot(t(Original_KMA2014_2015Strata_PercentbyLocus), 
          col.regions = new.colors, 
          at = seq(from = 0, to = 1, length.out = 100), 
          main = "% Genotyped", xlab = "SILLY", ylab = "Locus", 
          scales = list(x = list(rot = 90)), 
          aspect = "fill")  # aspect = "iso" will make squares


#### Check individuals
## View Histogram of Failure Rate by Strata
invisible(sapply(KMA2014_2015Strata, function(mix) {
  my.gcl <- get(paste(mix, ".gcl", sep = ''))
  failure <- apply(my.gcl$scores[, , 1], 1, function(ind) {sum(ind == "0") / length(ind)} )
  hist(x = failure, main = mix, xlab = "Failure Rate", col = 8, xlim = c(0, 1), ylim = c(0, 20), breaks = seq(from = 0, to = 1, by = 0.02))
  abline(v = 0.2, lwd = 3)
}))

## Experimental method to look for a way to remove "missing" loci fish based on a strata specific failure rate (rather than 0.8 globally, as done normally)
# Experimented with loess, smooth.split, polynomial regressions, outlier (Franz), and Cooks Distance methods to model outliers.
# Thought of a clever idea to look at the "shoulder" of the cumulative success rate, but this concept is a bit weird.
# Planning to keep it simple!
invisible(sapply(KMA2014_2015Strata, function(mix) {
  my.gcl <- get(paste(mix, ".gcl", sep = ''))
  success <- apply(my.gcl$scores[, , 1], 1, function(ind) {sum(!ind == "0") / length(ind)} )
  plot(x = sort(success, decreasing = TRUE), main = mix, xlab = "Sorted Individuatls", col = 8, type = "p", pch = 16, ylim = c(0, 1), ylab = "Success Rate")
  
  cutoff.floor = 0.9
  cutoff <- min(sort(success, decreasing = TRUE)[sort(success, decreasing = TRUE) > (seq(success) / length(success))])
  cutoff <- pmin(cutoff, cutoff.floor)
  points(y = sort(success, decreasing = TRUE)[sort(success, decreasing = TRUE) < cutoff],
         x = seq(success)[sort(success, decreasing = TRUE) < cutoff], pch = 16, col = "black")
  
  abline(h = 0.8, col = "red", lwd = 3)
  abline(h = cutoff.floor, col = "red", lwd = 3, lty = 3)
  segments(x0 = 0, x1 = length(success), y0 = 0, y1 = 1, lwd = 3)
  
  # points(smooth.spline(sort(success, decreasing = TRUE), spar = 0.2), type = "l", col = 1, lwd = 3)
  # 
  # y.loess <- loess(sort(success, decreasing = TRUE) ~ seq(success), span = 0.75, data.frame(x = seq(success), y = sort(success, decreasing = TRUE)))
  # y.predict <- predict(y.loess, data.frame(x = seq(success)))
  # lines(x = seq(success), y = y.predict, col = "black", lwd = 3)
  # 
  # y.loess <- loess(sort(success, decreasing = TRUE) ~ seq(success), span = 0.5, data.frame(x = seq(success), y = sort(success, decreasing = TRUE)))
  # y.predict <- predict(y.loess, data.frame(x = seq(success)))
  # lines(x = seq(success), y = y.predict, col = "blue", lwd = 3)
  # 
  # y.loess <- loess(sort(success, decreasing = TRUE) ~ seq(success), span = 0.3, data.frame(x = seq(success), y = sort(success, decreasing = TRUE)))
  # y.predict <- predict(y.loess, data.frame(x = seq(success)))
  # lines(x = seq(success), y = y.predict, col = "green", lwd = 3)
  # 
  # fit1 <- lm(sort(success, decreasing = TRUE) ~ poly(seq(success) + 100, 1, raw = TRUE))
  # fit2 <- lm(sort(success, decreasing = TRUE) ~ poly(seq(success) + 100, 2, raw = TRUE))
  # fit3 <- lm(sort(success, decreasing = TRUE) ~ poly(seq(success) + 100, 3, raw = TRUE))
  # fit4 <- lm(sort(success, decreasing = TRUE) ~ poly(seq(success) + 100, 4, raw = TRUE))
  # 
  # points(y = predict(fit1, data.frame(x = seq(success) + 100)), x = seq(success), type = "l", col = 1, lwd = 3)
  # points(y = predict(fit2, data.frame(x = seq(success) + 100)), x = seq(success), type = "l", col = 2, lwd = 3)
  # points(y = predict(fit3, data.frame(x = seq(success) + 100)), x = seq(success), type = "l", col = 3, lwd = 3)
  # points(y = predict(fit4, data.frame(x = seq(success) + 100)), x = seq(success), type = "l", col = 4, lwd = 3)
  
  text(x = 0, y = 0.5, labels = paste("cutoff", sum(success < cutoff), sep = "_"), cex = 1.2, pos = 4)
  text(x = 0, y = 0.4, labels = paste("cutoff.floor", sum(success < cutoff.floor), sep = "_"), pos = 4)
  text(x = 0, y = 0.3, labels = paste("0.8", sum(success < 0.8), sep = "_"), pos = 4)
  
}))


### Initial
## Get number of individuals per silly before removing missing loci individuals
Original_KMA2014_2015Strata_ColSize <- sapply(paste(KMA2014_2015Strata, ".gcl", sep = ''), function(x) get(x)$n)
KMA2014_2015Strata_SampleSizes[, "Genotyped"] <- Original_KMA2014_2015Strata_ColSize


### Missing
## Remove individuals with >20% missing data
KMA2014_2015Strata_MissLoci <- RemoveIndMissLoci.GCL(sillyvec = KMA2014_2015Strata, proportion = 0.8)

## Get number of individuals per silly after removing missing loci individuals
ColSize_KMA2014_2015Strata_PostMissLoci <- sapply(paste(KMA2014_2015Strata, ".gcl", sep = ''), function(x) get(x)$n)
KMA2014_2015Strata_SampleSizes[, "Missing"] <- Original_KMA2014_2015Strata_ColSize-ColSize_KMA2014_2015Strata_PostMissLoci


### Duplicate
## Check within collections for duplicate individuals.
KMA2014_2015Strata_DuplicateCheck95MinProportion <- CheckDupWithinSilly.GCL(sillyvec = KMA2014_2015Strata, loci = loci48, quantile = NULL, minproportion = 0.95)
KMA2014_2015Strata_DuplicateCheckReportSummary <- sapply(KMA2014_2015Strata, function(x) KMA2014_2015Strata_DuplicateCheck95MinProportion[[x]]$report)
KMA2014_2015Strata_DuplicateCheckReportSummary

## Remove duplicate individuals
KMA2014_2015Strata_RemovedDups <- RemoveDups.GCL(KMA2014_2015Strata_DuplicateCheck95MinProportion)

## Get number of individuals per silly after removing duplicate individuals
ColSize_KMA2014_2015Strata_PostDuplicate <- sapply(paste(KMA2014_2015Strata, ".gcl", sep = ''), function(x) get(x)$n)
KMA2014_2015Strata_SampleSizes[, "Duplicate"] <- ColSize_KMA2014_2015Strata_PostMissLoci-ColSize_KMA2014_2015Strata_PostDuplicate


### Final
KMA2014_2015Strata_SampleSizes[, "Final"] <- ColSize_KMA2014_2015Strata_PostDuplicate
KMA2014_2015Strata_SampleSizes

dir.create("Output")
write.xlsx(KMA2014_2015Strata_SampleSizes, file = "Output/KMA2014_2015Strata_SampleSizes.xlsx")
dput(x = KMA2014_2015Strata_SampleSizes, file = "Objects/KMA2014_2015Strata_SampleSizes.txt")



# dput postQC mixture sillys
dir.create("Raw genotypes/OriginalCollections_Strata_PostQC")
invisible(sapply(KMA2014_2015Strata, function(silly) {dput(x = get(paste(silly, ".gcl", sep = '')), file = paste("Raw genotypes/OriginalCollections_Strata_PostQC/" , silly, ".txt", sep = ''))} )); beep(8)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Clean workspace; dget .gcl objects and Locus Control ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls(all = TRUE))
setwd("V:/Analysis/4_Westward/Chinook/CSRI Westward Commercial Harvest 2014-2016/Mixtures")
# This sources all of the new GCL functions to this workspace
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")
source("H:/R Source Scripts/Functions.GCL_KS.R")

## Get objects
LocusControl <- dget(file = "Objects/OriginalLocusControl48.txt")

KMAobjects <- list.files(path = "Objects", recursive = FALSE)
KMAobjects <- KMAobjects[!KMAobjects %in% c("OriginalLocusControl48.txt")]
KMAobjects

invisible(sapply(KMAobjects, function(objct) {assign(x = unlist(strsplit(x = objct, split = ".txt")), value = dget(file = paste(getwd(), "Objects", objct, sep = "/")), pos = 1) })); beep(2)


## Get un-altered mixtures
invisible(sapply(KMA2014_2015Strata, function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/OriginalCollections_Strata_PostQC/", silly, ".txt", sep = "")), pos = 1)} )); beep(2)
objects(pattern = "\\.gcl")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Get/Create MSA Objects ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("V:/Analysis/4_Westward/Chinook/CSRI Westward Commercial Harvest 2014-2016/Baseline/Objects")

colors10 <- dget(file = "colors10.txt")
groups10 <- dget(file = "groups10.txt")
groups10short <- dget(file = "groups10short.txt")
groupvec10 <- dget(file = "groupvec10.txt")
KMA211Pops <- dget(file = "KMA211Pops.txt")
KMA211Pops42Loci.baseline <- dget(file = "KMA211Pops42Loci.baseline.txt")
loci42 <- dget(file = "loci42.txt")
popnames211 <- dget(file = "popnames211.txt")

setwd("V:/Analysis/4_Westward/Chinook/CSRI Westward Commercial Harvest 2014-2016/Mixtures/Objects")

dput(x = colors10, file = "colors10.txt")
dput(x = groups10, file = "groups10.txt")
dput(x = groups10short, file = "groups10short.txt")
dput(x = groupvec10, file = "groupvec10.txt")
dput(x = KMA211Pops, file = "KMA211Pops.txt")
dput(x = KMA211Pops42Loci.baseline, file = "KMA211Pops42Loci.baseline.txt")
dput(x = loci42, file = "loci42.txt")
dput(x = popnames211, file = "popnames211.txt")

# Starting with a Regionally flat prior, then rolling prior afterwards
KMA211Pops10FlatPrior <- Prior.GCL(groupvec = groupvec10, groupweights = rep(1 / 10, 10), minval = 0.01)
dput(x = KMA211Pops10FlatPrior, file = "KMA211Pops10FlatPrior.txt")

# Initial values
KMA211PopsInits <- MultiChainInits.GCL(npops = length(KMA211Pops), nchains = 5, prop = 0.9)
dput(x = KMA211PopsInits, file = "KMA211PopsInits.txt")

# Seeds
KMA211PopsChinookSeeds <- matrix(sample(seq(10000), 3 * 5), nrow = 3)
dput(x = KMA211PopsChinookSeeds, file = "KMA211PopsChinookSeeds.txt")

# Groups 2 Row
groups10tworows <- c("Russia", "CWAK\nYukon", "North\nPen", "Chignik", "Kodiak", "Cook\nInlet", "Copper", "SEAK", "British\nColumbia", "West\nCoast US")
dput(x = groups10tworows, file = "groups10tworows.txt")

setwd("V:/Analysis/4_Westward/Chinook/CSRI Westward Commercial Harvest 2014-2016/Mixtures")


dir.create("BAYES")
sapply(c("Baseline", "Control", "Mixture", "Output"), function(fldr) {dir.create(paste("BAYES", fldr, sep = "/"))})
file.copy(from = "V:/Analysis/4_Westward/Chinook/CSRI Westward Commercial Harvest 2014-2016/Baseline/BAYES/Baseline/KMA211Pops42Loci.bse", 
          to = "BAYES/Baseline/KMA211Pops42Loci.bse")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 1 MSA files for BAYES 2014 Early Strata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Indetify Strata to Run
KMA2014Strata_1_Early <- grep(pattern = "1_Early", x = KMA2014Strata, value = TRUE)
Round1Mixtures_2014 <- c(KMA2014Strata_1_Early, "KMAINC14_2_Late", "KSPENCHIG14")
dput(x = Round1Mixtures_2014, file = "Objects/Round1Mixtures_2014.txt")

## Dumping Mixture files
KMA21142MixtureFormat <- CreateMixture.GCL(sillys = KMA2014Strata_1_Early[1], loci = loci42, IDs = NULL, mixname = KMA2014Strata_1_Early[1],
                                           dir = "BAYES/Mixture", type = "BAYES", PT = FALSE)
dput(KMA21142MixtureFormat, file = "Objects/KMA21142MixtureFormat.txt")

sapply(Round1Mixtures_2014, function(Mix) {CreateMixture.GCL(sillys = Mix, loci = loci42, IDs = NULL, mixname = Mix, dir = "BAYES/Mixture", type = "BAYES", PT = FALSE)} )

## Dumping Control files
sapply(Round1Mixtures_2014, function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA211Pops, loci = loci42, mixname = Mix, basename = "KMA211Pops42Loci", suffix = "", nreps = 40000, nchains = 5,
                        groupvec = groupvec10, priorvec = KMA211Pops10FlatPrior, initmat = KMA211PopsInits, dir = "BAYES/Control",
                        seeds = KMA211PopsChinookSeeds, thin = c(1, 1, 100), mixfortran = KMA21142MixtureFormat, basefortran = KMA211Pops42Loci.baseline, switches = "F T F T T T F")
})

## Create output directory
sapply(Round1Mixtures_2014, function(Mix) {dir.create(paste("BAYES/Output/", Mix, sep = ""))})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Go run BAYES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Round 1 Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round1Mixtures_2014_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = seq(groups10), groupnames = groups10, 
                                                              maindir = "BAYES/Output", mixvec = Round1Mixtures_2014, prior = "",  
                                                              ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)

# Dput 1) estimates stats + posterior output & 2) estimates stats
dir.create("Estimates objects")
dput(Round1Mixtures_2014_Estimates, file = "Estimates objects/Round1Mixtures_2014_Estimates.txt")
dput(Round1Mixtures_2014_Estimates$Stats, file = "Estimates objects/Round1Mixtures_2014_EstimatesStats.txt")

Round1Mixtures_2014_Estimates <- dget(file = "Estimates objects/Round1Mixtures_2014_Estimates.txt")
Round1Mixtures_2014_EstimatesStats <- dget(file = "Estimates objects/Round1Mixtures_2014_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2
sapply(Round1Mixtures_2014_Estimates$Stats, function(Mix) {Mix[, "GR"]})
sapply(Round1Mixtures_2014_Estimates$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})
require(gplots)
sapply(Round1Mixtures_2014, function(Mix) {
  BarPlot <- barplot2(Round1Mixtures_2014_EstimatesStats[[Mix]][, "GR"], col = "blue", ylim = c(1, pmax(1.5, max(Round1Mixtures_2014_EstimatesStats[[Mix]][, "GR"]))), ylab = "Gelman-Rubin", type = "h", xpd = FALSE, main = Mix, names.arg = '')
  abline(h = 1.2, lwd = 3, xpd = FALSE)
  text(x = BarPlot, y = 1, labels = groups10tworows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
})

# Quick look at raw posterior output
str(Round1Mixtures_2014_Estimates$Output)
Round1Mixtures_2014_Header <- setNames(object = c("Southwest Kodiak / Alitak Early June 1-July 5, 2014",
                                                  "Eastside Kodiak Early June 1-July 5, 2014",
                                                  "Westwide Kodiak Early June 1-July 5, 2014",
                                                  "Mainland Kodiak Late July 6-August 5, 2014",
                                                  "South Peninsula / Chignik June 1-August 5, 2014"), 
                                       nm = Round1Mixtures_2014)
dput(x = Round1Mixtures_2014_Header, file = "Objects/Round1Mixtures_2014_Header.txt")


PlotPosterior <- function(mixvec = NULL, output, header = NULL, groups, colors = NULL, set.mfrow, thin = 10, chains = 5){
  if(is.null(colors)) {colors <- rep("black", length(groups))}
  if(is.null(mixvec)) {mixvec <- names(output)}
  
  par(mfrow = set.mfrow, mar = c(2.1, 2.1, 1.1, 1.1), oma = c(4.1, 4.1, 3.1, 1.1))
  
  invisible(sapply(mixvec, function(Mix) {
    invisible(sapply(seq(groups), function(i) {
      RG <- output[[Mix]][, i]
      RG <- RG[seq(1, length(RG), thin)]
      plot(RG, type = "l", ylim = c(0,1), xlab = "", ylab = "")
      abline(v = seq(0, length(RG), length(RG)/chains), xpd = FALSE)
      text(x = length(RG)/2, y = 0.96, labels = groups[i], col = colors[i], cex = 1.2, font = 2)} ))
    if(prod(set.mfrow) - length(groups) == 1) {plot.new()}
    if(prod(set.mfrow) - length(groups) == 2) {plot.new(); plot.new()}
    if(!is.null(header)) {
      mtext(text = header[Mix], side = 3, outer = TRUE, cex = 1.5)
    } else {
      mtext(text = Mix, side = 3, outer = TRUE, cex = 1.5)
    }
    mtext(text = paste("Iteration (", chains, " chain[s] each)", sep = ''), side = 1, outer = TRUE, cex = 1.5, line = 1.5)
    mtext(text = "Posterior", side = 2, outer = TRUE, cex = 1.5, line = 1.5)
  }))
  
  par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1), oma = c(0, 0, 0, 0))
}
dput(x = PlotPosterior, file = "Objects/PlotPosterior.txt")

PlotPosterior(mixvec = Round1Mixtures_2014, output = Round1Mixtures_2014_Estimates$Output, 
              groups = groups10, colors = colors10, 
              header = Round1Mixtures_2014_Header, set.mfrow = c(5, 2), thin = 10)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Round 1 Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Barplots
QuickBarplot <- function(mixvec, estimatesstats, groups, groups2rows = NULL, header) {
  while(!require(gplots, quietly = TRUE)){install.packages("gplots")}
  
  # if(is.null(groups2rows)) {groups2rows <- groups}
  if(length(mixvec) != length(header)) {stop("Header is not the same length as mixvec!")}
  for(i in seq(mixvec)) {
    header[i] <- paste(header[i], "\nn = ", get(paste(mixvec[i], ".gcl", sep = ''))$n, sep = '')
  }
  if("Stats" %in% names(estimatesstats)) {estimatesstats <- estimatesstats$Stats}
  
  par(mfrow = c(1, 1), mar = c(3.1, 5.1, 4.1, 2.1), oma = rep(0, 4))
  sapply(mixvec, function(Mix) {
    Barplot <- barplot2(height = estimatesstats[[Mix]][, "median"] * 100,
                        beside = TRUE, plot.ci = TRUE, ci.lwd = 1,
                        ci.l = estimatesstats[[Mix]][, "5%"] * 100,
                        ci.u = estimatesstats[[Mix]][, "95%"] * 100,
                        ylim = c(0, 100), col = "blue", yaxt = 'n', xaxt = 'n',
                        main = header[Mix], ylab = "Percent of Mixture", cex.lab = 2, cex.main = 2)
    axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = 1.5)
    abline(h = 0, xpd = FALSE)
    
    if(is.null(groups2rows)) {
      text(x = Barplot[, 1], y = -1, labels = groups, srt = 90, adj =  1, xpd = TRUE, cex = 0.5)
    } else {
      mtext(text = groups2rows, side = 1, line = 1, at = Barplot[, 1], adj = 0.5, cex = 0.6)
    }
  })
  par(mar = c(5.1, 4.1, 4.1, 2.1))
}

dput(x = QuickBarplot, file = "Objects/QuickBarplot.txt")

QuickBarplot(mixvec = Round1Mixtures_2014, estimatesstats = Round1Mixtures_2014_EstimatesStats, groups = groups10, groups2rows = groups10tworows, header = Round1Mixtures_2014_Header)


## Make violin plots of posteriors with RGs sorted

ViolinPlot <- function(mixvec = NULL, estimates, groups, colors, header, wex = 1, thin = 10) {
  while(!require(vioplot, quietly = TRUE)){install.packages("vioplot")}
  
  if(is.null(mixvec)) {mixvec <- names(estimates$Stats)}
  
  par(mar = c(5.6, 4.6, 3.6, 1.1))
  sapply(mixvec, function(Mix) {
    plot(estimates$Stats[[Mix]][, "median"], cex = 3, pch = 16, col = colors, ylab = "Proportion of Mixture", ylim = c(0, 1), xlab = "", axes = FALSE, main = header[[Mix]], cex.main = 2, cex.lab = 1.5)
    sapply(seq(groups), function(i) {vioplot(estimates$Output[[Mix]][seq(from = 1, to = nrow(estimates$Output[[Mix]]), by = thin), i], at = i, horizontal = FALSE, col = colors[i], border = TRUE, drawRect = FALSE, rectCol = colors[i], add = TRUE, wex = wex, lwd = 2)})
    arrows(x0 = seq(groups), y0 = estimates$Stats[[Mix]][, "5%"], x1 = seq(groups), y1 = estimates$Stats[[Mix]][, "95%"], angle = 90, code = 3, length = 0.2, lwd = 2)
    points(estimates$Stats[[Mix]][, "median"], cex = 2, pch = 21, col = "white", bg = colors, lwd = 3)
    axis(side = 2, lwd = 3, cex.axis = 1.5)
    text(x = (seq(groups)) - 0.35, y = 0, labels = groups, srt = 60, pos = 1, offset = 2.5, xpd = TRUE)
    axis(side = 1, labels = NA, at = seq(groups), pos = 0, lwd = 2, tick = FALSE)
    abline(h = 0, lwd = 3, xpd = FALSE)
  } )
  
}

dput(x = ViolinPlot, file = "Objects/ViolinPlot.txt")

ViolinPlot(estimates = Round1Mixtures_2014_Estimates, groups = groups10tworows, colors = colors10, header = Round1Mixtures_2014_Header)
rm(Round1Mixtures_2014_Estimates)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 2 MSA files for BAYES 2014 Late Strata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Indetify Strata to Run
KMA2014Strata_2_Late <- grep(pattern = "2_Late", x = KMA2014Strata, value = TRUE)
Round2Mixtures_2014 <- KMA2014Strata_2_Late[-which(KMA2014Strata_2_Late == "KMAINC14_2_Late")]  # remove mainland late, as already done
dput(x = Round2Mixtures_2014, file = "Objects/Round2Mixtures_2014.txt")

# Create rolling prior based on Round 1 estimates
Round2Mixtures_2014_Prior <- sapply(Round1Mixtures_2014_EstimatesStats[1:3], function(Mix) {Prior.GCL(groupvec = groupvec10, groupweights = Mix[, "mean"], minval = 0.01)}, simplify = FALSE)
names(Round2Mixtures_2014_Prior) <- gsub(pattern = "1_Early", replacement = "2_Late", x = names(Round2Mixtures_2014_Prior))  # This changes the names
dput(x = Round2Mixtures_2014_Prior, file = "Objects/Round2Mixtures_2014_Prior.txt")
str(Round2Mixtures_2014_Prior)


## Dumping Mixture files
sapply(Round2Mixtures_2014, function(Mix) {CreateMixture.GCL(sillys = Mix, loci = loci42, IDs = NULL, mixname = Mix, dir = "BAYES/Mixture", type = "BAYES", PT = FALSE)} )

## Dumping Control files
sapply(Round2Mixtures_2014, function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA211Pops, loci = loci42, mixname = Mix, basename = "KMA211Pops42Loci", suffix = "", nreps = 40000, nchains = 5,
                        groupvec = groupvec10, priorvec = Round2Mixtures_2014_Prior[[Mix]], initmat = KMA211PopsInits, dir = "BAYES/Control",
                        seeds = KMA211PopsChinookSeeds, thin = c(1, 1, 100), mixfortran = KMA21142MixtureFormat, basefortran = KMA211Pops42Loci.baseline, switches = "F T F T T T F")
})

## Create output directory
sapply(Round2Mixtures_2014, function(Mix) {dir.create(paste("BAYES/Output/", Mix, sep = ""))})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Go run BAYES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Round 2 Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round2Mixtures_2014_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = seq(groups10), groupnames = groups10, 
                                                              maindir = "BAYES/Output", mixvec = Round2Mixtures_2014, prior = "",  
                                                              ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)

# Dput 1) estimates stats + posterior output & 2) estimates stats
dput(Round2Mixtures_2014_Estimates, file = "Estimates objects/Round2Mixtures_2014_Estimates.txt")
dput(Round2Mixtures_2014_Estimates$Stats, file = "Estimates objects/Round2Mixtures_2014_EstimatesStats.txt")

Round2Mixtures_2014_Estimates <- dget(file = "Estimates objects/Round2Mixtures_2014_Estimates.txt")
Round2Mixtures_2014_EstimatesStats <- dget(file = "Estimates objects/Round2Mixtures_2014_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2
sapply(Round2Mixtures_2014_Estimates$Stats, function(Mix) {Mix[, "GR"]})
sapply(Round2Mixtures_2014_Estimates$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})
require(gplots)
sapply(Round2Mixtures_2014, function(Mix) {
  BarPlot <- barplot2(Round2Mixtures_2014_EstimatesStats[[Mix]][, "GR"], col = "blue", ylim = c(1, pmax(1.5, max(Round2Mixtures_2014_EstimatesStats[[Mix]][, "GR"]))), ylab = "Gelman-Rubin", type = "h", xpd = FALSE, main = Mix, names.arg = '')
  abline(h = 1.2, lwd = 3, xpd = FALSE)
  text(x = BarPlot, y = 1, labels = groups10tworows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
})

# Quick look at raw posterior output
str(Round2Mixtures_2014_Estimates$Output)
Round2Mixtures_2014_Header <- setNames(object = c("Southwest Kodiak / Alitak Late July 6-August 5, 2014",
                                                  "Eastside Kodiak Late July 6-August 5, 2014",
                                                  "Westwide Kodiak Late July 6-August 5, 2014"), 
                                       nm = Round2Mixtures_2014)
dput(x = Round2Mixtures_2014_Header, file = "Objects/Round2Mixtures_2014_Header.txt")


PlotPosterior(mixvec = Round2Mixtures_2014, output = Round2Mixtures_2014_Estimates$Output, 
              groups = groups10, colors = colors10, 
              header = Round2Mixtures_2014_Header, set.mfrow = c(5, 2), thin = 10)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Round 2 Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Barplots
QuickBarplot(mixvec = Round2Mixtures_2014, estimatesstats = Round2Mixtures_2014_EstimatesStats, groups = groups10, groups2rows = groups10tworows, header = Round2Mixtures_2014_Header)

## Make violin plots of posteriors with RGs sorted
ViolinPlot(estimates = Round2Mixtures_2014_Estimates, groups = groups10tworows, colors = colors10, header = Round2Mixtures_2014_Header)
rm(Round2Mixtures_2014_Estimates)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Create Final Estimates Objects for Each Temporal Strata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Round1Mixtures_2014_EstimatesStats <- dget("Estimates objects/Round1Mixtures_2014_EstimatesStats.txt")
Round2Mixtures_2014_EstimatesStats <- dget("Estimates objects/Round2Mixtures_2014_EstimatesStats.txt")

# View as tables by year
require(reshape)
samp.df.2014 <- data.frame(t(sapply(KMA2014Strata, function(strata) {
  location <- unlist(strsplit(x = strata, split = "_"))[1]
  temporal.strata <- as.numeric(unlist(strsplit(x = strata, split = "_"))[2])
  n <- get(paste(strata, ".gcl", sep = ""))$n
  c(location = location, temporal.strata = temporal.strata, n = n)
})))
cast(data = samp.df.2014, location ~ temporal.strata)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1_Early
str(Round1Mixtures_2014_EstimatesStats)
KMA2014Strata_1_Early_EstimatesStats <- Round1Mixtures_2014_EstimatesStats[1:3]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2_Late
str(Round2Mixtures_2014_EstimatesStats)
KMA2014Strata_2_Late_EstimatesStats <- c(Round2Mixtures_2014_EstimatesStats,
                                         Round1Mixtures_2014_EstimatesStats[4])
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SouthPen/Chignik
KSPENCHIG2014_EstimatesStats <- Round1Mixtures_2014_EstimatesStats[5]


str(KMA2014Strata_1_Early_EstimatesStats)
str(KMA2014Strata_2_Late_EstimatesStats)
str(KSPENCHIG2014_EstimatesStats)


dir.create("Estimates objects/Final")
dput(x = KMA2014Strata_1_Early_EstimatesStats, file = "Estimates objects/Final/KMA2014Strata_1_Early_EstimatesStats.txt")
dput(x = KMA2014Strata_2_Late_EstimatesStats, file = "Estimates objects/Final/KMA2014Strata_2_Late_EstimatesStats.txt")
dput(x = KSPENCHIG2014_EstimatesStats, file = "Estimates objects/Final/KSPENCHIG2014_EstimatesStats.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot stock composition results 2014 KMA Mixtures ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
KMA2014Strata_1_Early_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2014Strata_1_Early_EstimatesStats.txt")
KMA2014Strata_2_Late_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2014Strata_2_Late_EstimatesStats.txt")
KSPENCHIG2014_EstimatesStats <- dget(file = "Estimates objects/Final/KSPENCHIG2014_EstimatesStats.txt")

KMA2014Strata_EstimatesStats <- c(KMA2014Strata_1_Early_EstimatesStats, 
                                  KMA2014Strata_2_Late_EstimatesStats, 
                                  KSPENCHIG2014_EstimatesStats)
dput(x = KMA2014Strata_EstimatesStats, file = "Estimates objects/Final/KMA2014Strata_EstimatesStats.txt")
KMA2014Strata_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2014Strata_EstimatesStats.txt")

str(KMA2014Strata_EstimatesStats)


KMA2014 <- c("KALITC14", "KEASTC14", "KWESTC14", "KMAINC14")
dput(KMA2014, file = "Objects/KMA2014.txt")
TempMix14 <- sapply(KMA2014, function(geo) {grep(pattern = geo, x = names(KMA2014Strata_EstimatesStats), value = TRUE)} )

Legend14 <- setNames(object = c("June 1-July 5", "July 6-August 5"), 
                     nm = c("1_Early", "2_Late"))
TempLegend14 <- sapply(KMA2014, function(geo) {
  Legend14[sapply(TempMix14[[geo]], function(strata) {unlist(strsplit(x = strata, split = paste(geo, "_", sep = '')))[2]} )]
} )

GeoHeader <- setNames(object = c(paste0("SW Kodiak/Alitak 255, 256, 257"),
                                 paste0("Eastside Kodiak/Afognak 252, 258, 259"),
                                 paste0("Westside Kodiak/Afognak 251, 253, 254"),
                                 paste0("Mainland 262")),
                      nm = unlist(strsplit(x = KMA2014, split = "14")))

ProportionColors <- colorpanel(n = 2, low = "blue", high = "white")
TempProportionColors14 <- sapply(KMA2014, function(geo) {
  ProportionColors[sapply(TempMix14[[geo]], function(strata) {as.numeric(unlist(strsplit(x = strata, split = "_"))[2])} )]
} )

Estimates <- KMA2014Strata_EstimatesStats
Groups <- groups10
Groups2Rows <- groups10tworows
cex.lab <- 1.5
cex.yaxis <- 1.3
cex.xaxis <- 0.6
cex.main <- 1.7
cex.leg <- 1.3
ci.lwd <- 2.5

dir.create("Figures")
dir.create("Figures/2014")
require(devEMF)

sapply(names(TempMix14), function(geomix) {
  emf(file = paste("Figures/2014/", geomix, ".emf", sep = ''), width = 8.5, height = 6.5, family = "sans", bg = "white")
  par(mar = c(2.1, 4.1, 2.6, 0.6))
  
  Barplot1 <- barplot2(height = t(sapply(TempMix14[[geomix]], function(tempmix) {Estimates[[tempmix]][, "median"]})) * 100, 
                       beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                       ci.l = t(sapply(TempMix14[[geomix]], function(tempmix) {Estimates[[tempmix]][, "5%"]})) * 100, 
                       ci.u = t(sapply(TempMix14[[geomix]], function(tempmix) {Estimates[[tempmix]][, "95%"]})) * 100, 
                       ylim = c(0, 100), col = TempProportionColors14[[geomix]], cex.axis = cex.yaxis, yaxt = "n", xaxt = 'n')
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
  legend(legend = TempLegend14[[geomix]], x = "topleft", fill = TempProportionColors14[[geomix]], border = "black", bty = "n", cex = cex.leg, title="2014")
  abline(h = 0, xpd = FALSE)
  
  mtext(text = "Percentage of Catch", side = 2, cex = cex.yaxis, line = 3)
  mtext(text = Groups2Rows, side = 1, line = 1, at = apply(Barplot1, 2, mean), adj = 0.5, cex = cex.xaxis)
  mtext(text = GeoHeader[unlist(strsplit(geomix, split = "14"))], side = 3, cex = cex.main, line = 1)
  dev.off()
})





Estimates <- KMA2014Strata_EstimatesStats
geomix = "KSPENCHIG14"

emf(file = paste("Figures/2014/", geomix, ".emf", sep = ''), width = 8.5, height = 6.5, family = "sans", bg = "white")
par(mar = c(2.1, 4.1, 2.6, 0.6))

Barplot1 <- barplot2(height = t(Estimates[[geomix]][, "median"]) * 100, 
                     beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                     ci.l = t(Estimates[[geomix]][, "5%"]) * 100, 
                     ci.u = t(Estimates[[geomix]][, "95%"]) * 100, 
                     ylim = c(0, 100), col = "blue", cex.axis = cex.yaxis, yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
legend(legend = "June 1-August 5", x = "topleft", fill = "blue", border = "black", bty = "n", cex = cex.leg, title="2014")
abline(h = 0, xpd = FALSE)

mtext(text = "Percentage of Catch", side = 2, cex = cex.yaxis, line = 3)
mtext(text = Groups2Rows, side = 1, line = 1, at = apply(Barplot1, 2, mean), adj = 0.5, cex = cex.xaxis)
mtext(text = paste0("South Peninsula/Chignik Outside 282", "\u2013", "285; 272, 273, 275"), side = 3, cex = cex.main, line = 1)
dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 2014 Harvest Data ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
harvest14 <- read.csv(file = "Harvest/Salmon Catch by Day and Stat Area 2014 Date.csv", as.is = TRUE)
str(harvest14)

# Convert Date
harvest14$Date.Landed <- as.Date(harvest14$Date.Landed, format = "%Y-%m-%d")
harvest14$Date.Fishing.Began <- as.Date(harvest14$Date.Fishing.Began, format = "%Y-%m-%d")
harvest14$Date.Fishing.Ended <- as.Date(harvest14$Date.Fishing.Ended, format = "%Y-%m-%d")

# Create Temporal Strata
harvest14$Strata <- ifelse(harvest14$Date.Fishing.Began >= as.Date("2014-06-01") & harvest14$Date.Fishing.Began <= as.Date("2014-07-05"),
                           "Early",
                           ifelse(harvest14$Date.Fishing.Began >= as.Date("2014-07-06") & harvest14$Date.Fishing.Began <= as.Date("2014-08-05"),
                                  "Late",
                                  NA))
table(harvest14$Strata)

# Create Geographic Strata
unique(harvest14$Stat.Area)

Stat.Area.Eastside <- unique(harvest14$Stat.Area)[c(
  which(unique(harvest14$Stat.Area) >= 25800 & unique(harvest14$Stat.Area) <=25999),
  which(unique(harvest14$Stat.Area) >= 25200 & unique(harvest14$Stat.Area) <=25299))]

Stat.Area.Westside <- unique(harvest14$Stat.Area)[c(
  which(unique(harvest14$Stat.Area) >= 25100 & unique(harvest14$Stat.Area) <=25199),
  which(unique(harvest14$Stat.Area) >= 25300 & unique(harvest14$Stat.Area) <=25499))]

Stat.Area.SWAlitak <- unique(harvest14$Stat.Area)[c(
  which(unique(harvest14$Stat.Area) >= 25500 & unique(harvest14$Stat.Area) <=25799))]

Stat.Area.Mainland <- unique(harvest14$Stat.Area)[c(
  which(unique(harvest14$Stat.Area) >= 26200))]


length(c(Stat.Area.Eastside, Stat.Area.Westside, Stat.Area.SWAlitak, Stat.Area.Mainland)); length(unique(harvest14$Stat.Area))
unique(harvest14$Stat.Area)[!unique(harvest14$Stat.Area) %in% c(Stat.Area.Eastside, Stat.Area.Westside, Stat.Area.SWAlitak, Stat.Area.Mainland)]


harvest14$Geo <- ifelse(harvest14$Stat.Area %in% Stat.Area.Mainland, "Mainland",
                        ifelse(harvest14$Stat.Area %in% Stat.Area.SWAlitak, "SWAlitak",
                               ifelse(harvest14$Stat.Area %in% Stat.Area.Eastside, "Eastside",
                                      ifelse(harvest14$Stat.Area %in% Stat.Area.Westside, "Westside", NA
                                      ))))

table(harvest14$Strata, harvest14$Geo)

require(plyr)
ddply(.data = harvest14, ~Geo+Strata, summarise, harvest = sum(Number))
daply(.data = harvest14, ~Geo+Strata, summarise, harvest = sum(Number))

aggregate(x = harvest14$Number, by = list(harvest14$Strata, harvest14$Geo), FUN = sum, simplify = FALSE)
cast(aggregate(Number ~ Geo + Strata, data = harvest14, sum), Geo ~ Strata, value = "Number")

md <- melt(data = harvest14, id.vars = c("Strata", "Geo"), measure.vars = "Number", na.rm = TRUE)
cast(md, Geo~Strata, sum)


HarvestByStrata2014 <- data.matrix(cast(aggregate(Number ~ Geo + Strata, data = harvest14, sum), Geo ~ Strata, value = "Number")[c(3, 1, 4, 2), 2:3])
dimnames(HarvestByStrata2014) <- list(KMA2014, c("1_Early", "2_Late"))
HarvestByStrata2014
dput(x = HarvestByStrata2014, file = "Objects/HarvestByStrata2014.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot stock specific harvest results 2014 KMA Mixtures ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Multiply by harvest
HarvestByStrata2014


KMA2014Strata_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2014Strata_EstimatesStats.txt")
str(KMA2014Strata_EstimatesStats)
names(KMA2014Strata_EstimatesStats)
dimnames(KMA2014Strata_EstimatesStats[[1]])


KMA2014Strata_HarvestEstimatesStats <- sapply(names(KMA2014Strata_EstimatesStats)[-8], function(strata) {
  strata.split <- unlist(strsplit(x = strata, split = "_"))
  strata.split <- c(strata.split[1], paste(c(strata.split[2], strata.split[3]), collapse = "_"))
  
  cbind(KMA2014Strata_EstimatesStats[[strata]][, c("mean", "sd", "median", "5%", "95%")] * HarvestByStrata2014[strata.split[1], strata.split[2]],
        KMA2014Strata_EstimatesStats[[strata]][, c("P=0", "GR")])
}, simplify = FALSE )

dput(x = KMA2014Strata_HarvestEstimatesStats, file = "Estimates objects/Final/KMA2014Strata_HarvestEstimatesStats.txt")
KMA2014Strata_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2014Strata_HarvestEstimatesStats.txt")
str(KMA2014Strata_HarvestEstimatesStats)

# What should ymax be?
max(sapply(KMA2014Strata_HarvestEstimatesStats, function(strata) strata[, "95%"]))


TempMix14 <- sapply(KMA2014, function(geo) {grep(pattern = geo, x = names(KMA2014Strata_EstimatesStats), value = TRUE)} )

Legend14 <- setNames(object = c("June 1-July 5", "July 6-August 5"), 
                     nm = c("1_Early", "2_Late"))
TempLegend14 <- sapply(KMA2014, function(geo) {
  Legend14[sapply(TempMix14[[geo]], function(strata) {unlist(strsplit(x = strata, split = paste(geo, "_", sep = '')))[2]} )]
} )

GeoHeader <- setNames(object = c(paste0("SW Kodiak/Alitak 255, 256, 257"),
                                 paste0("Eastside Kodiak/Afognak 252, 258, 259"),
                                 paste0("Westside Kodiak/Afognak 251, 253, 254"),
                                 paste0("Mainland 262")),
                      nm = unlist(strsplit(x = KMA2014, split = "14")))


HarvestColors <- colorpanel(n = 2, low = "green", high = "white")
TempHarvestColors14 <- sapply(KMA2014, function(geo) {
  HarvestColors[sapply(TempMix14[[geo]], function(strata) {as.numeric(unlist(strsplit(x = strata, split = "_"))[2])} )]
} )

Estimates <- KMA2014Strata_HarvestEstimatesStats
Groups <- groups10
Groups2Rows <- groups10tworows
cex.lab <- 1.5
cex.yaxis <- 1.3
cex.xaxis <- 0.6
cex.main <- 1.7
cex.leg <- 1.3
ci.lwd <- 2.5
ymax <- 1500

sapply(names(TempMix14), function(geomix) {
  emf(file = paste("Figures/2014/", geomix, "Harvest.emf", sep = ''), width = 8.5, height = 6.5, family = "sans", bg = "white")
  par(mar = c(2.1, 4.1, 2.6, 0.6))
  
  Barplot1 <- barplot2(height = t(sapply(TempMix14[[geomix]], function(tempmix) {Estimates[[tempmix]][, "median"]})), 
                       beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                       ci.l = t(sapply(TempMix14[[geomix]], function(tempmix) {Estimates[[tempmix]][, "5%"]})), 
                       ci.u = t(sapply(TempMix14[[geomix]], function(tempmix) {Estimates[[tempmix]][, "95%"]})), 
                       ylim = c(0, ymax), col = TempHarvestColors14[[geomix]], cex.axis = cex.yaxis, yaxt = "n", xaxt = 'n')
  axis(side = 2, at = seq(0, ymax, 500), labels = formatC(x = seq(0, ymax, 500), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
  legend(legend = TempLegend14[[geomix]], x = "topleft", fill = TempHarvestColors14[[geomix]], border = "black", bty = "n", cex = cex.leg, title="2014")
  abline(h = 0, xpd = FALSE)
  
  mtext(text = "Number of Fish Harvested", side = 2, cex = cex.yaxis, line = 3)
  mtext(text = Groups2Rows, side = 1, line = 1, at = apply(Barplot1, 2, mean), adj = 0.5, cex = cex.xaxis)
  mtext(text = GeoHeader[unlist(strsplit(geomix, split = "14"))], side = 3, cex = cex.main, line = 1)
  dev.off()
})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Stratified Regional Roll-Ups 2014 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
HarvestByStrata2014


# For this to work properly, the output mixture directories all need to be in the same place (can be moved back afterwards)
sapply(KMA2014[1:3], function(geomix) {
  assign(x = paste(geomix, "_Annual_Stratified", sep = ''), 
         value = StratifiedEstimator.GCL(
           groupvec = seq(groups10), groupnames = groups10, maindir = "BAYES/Output", 
           mixvec = paste(geomix, colnames(HarvestByStrata2014)[!is.na(HarvestByStrata2014[geomix,])], sep = "_"), 
           catchvec = HarvestByStrata2014[geomix, !is.na(HarvestByStrata2014[geomix,])], 
           newname = paste(geomix, "_Annual_Stratified", sep = ''), priorname = '',
           ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1), 
         pos = 1)
  dput(x = get(paste(geomix, "_Annual_Stratified", sep = '')), file = paste("Estimates objects/", geomix, "_Annual_Stratified.txt", sep = ''))
} ); beep(5)

str(KALITC14_Annual_Stratified)

Round1Mixtures_2014_Estimates <- dget(file = "Estimates objects/Round1Mixtures_2014_Estimates.txt")
KMAINC14_Annual_Stratified <- list(Stats = Round1Mixtures_2014_Estimates$Stats$KMAINC14_2_Late,
                                   Output = Round1Mixtures_2014_Estimates$Output$KMAINC14_2_Late)
dput(x = KMAINC14_Annual_Stratified, file = "Estimates objects/KMAINC14_Annual_Stratified.txt")

# Create a list object with all Stratified Annual Rollups for the "Estimates objects/Final" folder
KMA2014_Annual_EstimatesStats <- sapply(KMA2014, function(geomix) {
  Stats <- get(paste(geomix, "_Annual_Stratified", sep = ''))$Stats
  Stats
}, simplify = FALSE)
str(KMA2014_Annual_EstimatesStats)
dput(x = KMA2014_Annual_EstimatesStats, file = "Estimates objects/Final/KMA2014_Annual_EstimatesStats.txt")




# Create a list object with all Stratified Annual Harvest
KMA2014_Annual_HarvestEstimatesStats <- sapply(KMA2014, function(strata) {
  cbind(KMA2014_Annual_EstimatesStats[[strata]][, c("mean", "sd", "median", "5%", "95%")] * sum(HarvestByStrata2014[strata, ], na.rm = TRUE),
        KMA2014_Annual_EstimatesStats[[strata]][, c("P=0", "GR")])
}, simplify = FALSE )

dput(x = KMA2014_Annual_HarvestEstimatesStats, file = "Estimates objects/Final/KMA2014_Annual_HarvestEstimatesStats.txt")
str(KMA2014_Annual_HarvestEstimatesStats)






# Create a matrix of annual means
Annual2014_Stratified_Estimates <- sapply(KMA2014, function(geomix) {
  KMA2014_Annual_EstimatesStats[[geomix]][, "mean"]
})

# Create a matrix of early strata means
KMA2014Strata_1_Early_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2014Strata_1_Early_EstimatesStats.txt")
EarlyStrata2014_Estimates <- sapply(KMA2014Strata_1_Early_EstimatesStats, function(geomix) {
  geomix[, "mean"]
})
colnames(EarlyStrata2014_Estimates) <- KMA2014[-4]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot Annual vs. Early Means
par(mar = c(4.1, 5.1, 3.1, 1.1))
sapply(KMA2014[-4], function(geomix) {
  Barplot <- barplot2(height = rbind(Annual2014_Stratified_Estimates[, geomix], EarlyStrata2014_Estimates[, geomix]) * 100, 
                      beside = TRUE, col = c("blue", "white"), yaxt = 'n', xaxt = 'n', main = geomix, ylab = "Precent of Mixture",
                      cex.lab = 2, cex.main = 2, ylim = c(0, 100))
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = ",", digits = 0, format = "f"), cex.axis = 1.5)
  abline(h = 0, xpd = FALSE)
  text(x = colMeans(Barplot), y = -1, labels = groups10tworows, srt = 90, adj = 1, xpd = TRUE, cex = 0.6)
  legend("topleft", legend = c("Annual Means", "Early Strata Means"), fill = c("blue", "white"), bty = 'n', cex = 1.5)
})




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 1 MSA files for BAYES 2015 Early Strata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require(reshape)
samp.df.2015 <- data.frame(t(sapply(KMA2015Strata, function(strata) {
  location <- unlist(strsplit(x = strata, split = "_"))[1]
  temporal.strata <- as.numeric(unlist(strsplit(x = strata, split = "_"))[2])
  n <- get(paste(strata, ".gcl", sep = ""))$n
  c(location = location, temporal.strata = temporal.strata, n = n)
})))
cast(data = samp.df.2015, location ~ temporal.strata)

# Indetify Strata to Run
KMA2015Strata_1_Early <- grep(pattern = "1_Early", x = KMA2015Strata, value = TRUE)
Round1Mixtures_2015 <- c(KMA2015Strata_1_Early, "KMAINC15_2_Late")
dput(x = Round1Mixtures_2015, file = "Objects/Round1Mixtures_2015.txt")

# Create rolling prior based on 2014 Round 1 estimates
Round1Mixtures_2014_EstimatesStats <- dget(file = "Estimates objects/Round1Mixtures_2014_EstimatesStats.txt")

Round1Mixtures_2015_Prior <- sapply(Round1Mixtures_2014_EstimatesStats[1:4], function(Mix) {
  Prior.GCL(groupvec = groupvec10, groupweights = Mix[, "mean"], minval = 0.01)}, simplify = FALSE)
names(Round1Mixtures_2015_Prior) <- gsub(pattern = "C14", replacement = "C15", x = names(Round1Mixtures_2015_Prior))  # This changes the names
dput(x = Round1Mixtures_2015_Prior, file = "Objects/Round1Mixtures_2015_Prior.txt")
str(Round1Mixtures_2015_Prior)

# Verify
sapply(Round1Mixtures_2015, function(geomix) {plot(as.vector(Round1Mixtures_2015_Prior[[geomix]]), type = "h", main = geomix)})

## Dumping Mixture files
sapply(Round1Mixtures_2015, function(Mix) {CreateMixture.GCL(sillys = Mix, loci = loci42, IDs = NULL, mixname = Mix, dir = "BAYES/Mixture", type = "BAYES", PT = FALSE)} )

## Dumping Control files
sapply(Round1Mixtures_2015, function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA211Pops, loci = loci42, mixname = Mix, basename = "KMA211Pops42Loci", suffix = "", nreps = 40000, nchains = 5,
                        groupvec = groupvec10, priorvec = Round1Mixtures_2015_Prior[[Mix]], initmat = KMA211PopsInits, dir = "BAYES/Control",
                        seeds = KMA211PopsChinookSeeds, thin = c(1, 1, 100), mixfortran = KMA21142MixtureFormat, basefortran = KMA211Pops42Loci.baseline, switches = "F T F T T T F")
})

## Create output directory
sapply(Round1Mixtures_2015, function(Mix) {dir.create(paste("BAYES/Output/", Mix, sep = ""))})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Go run BAYES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Round 1 2015 Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round1Mixtures_2015_Estimates <- CustomCombineBAYESOutput.GCL(
  groupvec = seq(groups10), groupnames = groups10 ,
  maindir = "BAYES/Output", 
  mixvec = Round1Mixtures_2015, prior = "",  
  ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)

# Dput 1) estimates stats + posterior output & 2) estimates stats
dput(Round1Mixtures_2015_Estimates, file = "Estimates objects/Round1Mixtures_2015_Estimates.txt")
dput(Round1Mixtures_2015_Estimates$Stats, file = "Estimates objects/Round1Mixtures_2015_EstimatesStats.txt")

Round1Mixtures_2015_Estimates <- dget(file = "Estimates objects/Round1Mixtures_2015_Estimates.txt")
Round1Mixtures_2015_EstimatesStats <- dget(file = "Estimates objects/Round1Mixtures_2015_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2
sapply(Round1Mixtures_2015_Estimates$Stats, function(Mix) {Mix[, "GR"]})
sapply(Round1Mixtures_2015_Estimates$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})
require(gplots)
par(mfrow = c(1, 1), mar = c(3.1, 4.6, 3.1, 1.1), oma = rep(0, 4))
sapply(Round1Mixtures_2015, function(Mix) {
  BarPlot <- barplot2(Round1Mixtures_2015_EstimatesStats[[Mix]][, "GR"], col = "blue", 
                      ylim = c(1, pmax(1.5, max(Round1Mixtures_2015_EstimatesStats[[Mix]][, "GR"]))), ylab = "Gelman-Rubin", 
                      xpd = FALSE, main = Mix, names.arg = '', cex.lab = 1.5, cex.main = 2)
  abline(h = 1.2, lwd = 3, xpd = FALSE)
  text(x = BarPlot, y = 1, labels = groups10tworows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
})

# Quick look at raw posterior output
str(Round1Mixtures_2015_Estimates$Output)
Round1Mixtures_2015_Header <- setNames(object = c("Southwest Kodiak / Alitak Early June 1-July 5, 2015",
                                                  "Eastside Kodiak Early June 1-July 5, 2015",
                                                  "Westwide Kodiak Early June 1-July 5, 2015",
                                                  "Mainland Kodiak Late July 6-August 5, 2015"), 
                                       nm = Round1Mixtures_2015)
dput(x = Round1Mixtures_2015_Header, file = "Objects/Round1Mixtures_2015_Header.txt")

PlotPosterior(mixvec = Round1Mixtures_2015, output = Round1Mixtures_2015_Estimates$Output, 
              groups = groups10tworows, colors = colors10, 
              header = Round1Mixtures_2015_Header, set.mfrow = c(5, 2), thin = 10)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Round 1 2015 Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Quick barplot
QuickBarplot(mixvec = Round1Mixtures_2015, estimatesstats = Round1Mixtures_2015_EstimatesStats, groups = groups10, groups2rows = groups10tworows, header = Round1Mixtures_2015_Header)

## Make violin plots of posteriors with RGs sorted
ViolinPlot(mixvec = Round1Mixtures_2015, estimates = Round1Mixtures_2015_Estimates, groups = groups10tworows, colors = colors10, header = Round1Mixtures_2015_Header)
rm(Round1Mixtures_2015_Estimates)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 2 MSA files for BAYES 2015 Late Strata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Indetify Strata to Run
KMA2015Strata_2_Late <- grep(pattern = "2_Late", x = KMA2015Strata, value = TRUE)
Round2Mixtures_2015 <- KMA2015Strata_2_Late[-which(KMA2015Strata_2_Late == "KMAINC15_2_Late")]  # remove mainland late, as already done
dput(x = Round2Mixtures_2015, file = "Objects/Round2Mixtures_2015.txt")

# Create rolling prior based on Round 1 estimates
Round2Mixtures_2015_Prior <- sapply(Round1Mixtures_2015_EstimatesStats[1:3], function(Mix) {Prior.GCL(groupvec = groupvec10, groupweights = Mix[, "mean"], minval = 0.01)}, simplify = FALSE)
names(Round2Mixtures_2015_Prior) <- gsub(pattern = "1_Early", replacement = "2_Late", x = names(Round2Mixtures_2015_Prior))  # This changes the names
dput(x = Round2Mixtures_2015_Prior, file = "Objects/Round2Mixtures_2015_Prior.txt")
str(Round2Mixtures_2015_Prior)

## Dumping Mixture files
sapply(Round2Mixtures_2015, function(Mix) {CreateMixture.GCL(sillys = Mix, loci = loci42, IDs = NULL, mixname = Mix, dir = "BAYES/Mixture", type = "BAYES", PT = FALSE)} )

## Dumping Control files
sapply(Round2Mixtures_2015, function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA211Pops, loci = loci42, mixname = Mix, basename = "KMA211Pops42Loci", suffix = "", nreps = 40000, nchains = 5,
                        groupvec = groupvec10, priorvec = Round2Mixtures_2015_Prior[[Mix]], initmat = KMA211PopsInits, dir = "BAYES/Control",
                        seeds = KMA211PopsChinookSeeds, thin = c(1, 1, 100), mixfortran = KMA21142MixtureFormat, basefortran = KMA211Pops42Loci.baseline, switches = "F T F T T T F")
})

## Create output directory
sapply(Round2Mixtures_2015, function(Mix) {dir.create(paste("BAYES/Output/", Mix, sep = ""))})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Go run BAYES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Round 2 Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round2Mixtures_2015_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = seq(groups10), groupnames = groups10, 
                                                              maindir = "BAYES/Output", mixvec = Round2Mixtures_2015, prior = "",  
                                                              ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)

# Dput 1) estimates stats + posterior output & 2) estimates stats
dput(Round2Mixtures_2015_Estimates, file = "Estimates objects/Round2Mixtures_2015_Estimates.txt")
dput(Round2Mixtures_2015_Estimates$Stats, file = "Estimates objects/Round2Mixtures_2015_EstimatesStats.txt")

Round2Mixtures_2015_Estimates <- dget(file = "Estimates objects/Round2Mixtures_2015_Estimates.txt")
Round2Mixtures_2015_EstimatesStats <- dget(file = "Estimates objects/Round2Mixtures_2015_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2
sapply(Round2Mixtures_2015_Estimates$Stats, function(Mix) {Mix[, "GR"]})
sapply(Round2Mixtures_2015_Estimates$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})
require(gplots)
sapply(Round2Mixtures_2015, function(Mix) {
  BarPlot <- barplot2(Round2Mixtures_2015_EstimatesStats[[Mix]][, "GR"], col = "blue", ylim = c(1, pmax(1.5, max(Round2Mixtures_2015_EstimatesStats[[Mix]][, "GR"]))), ylab = "Gelman-Rubin", type = "h", xpd = FALSE, main = Mix, names.arg = '')
  abline(h = 1.2, lwd = 3, xpd = FALSE)
  text(x = BarPlot, y = 1, labels = groups10tworows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
})

# Quick look at raw posterior output
str(Round2Mixtures_2015_Estimates$Output)
Round2Mixtures_2015_Header <- setNames(object = c("Southwest Kodiak / Alitak Late July 6-August 5, 2015",
                                                  "Eastside Kodiak Late July 6-August 5, 2015",
                                                  "Westwide Kodiak Late July 6-August 5, 2015"), 
                                       nm = Round2Mixtures_2015)
dput(x = Round2Mixtures_2015_Header, file = "Objects/Round2Mixtures_2015_Header.txt")


PlotPosterior(mixvec = Round2Mixtures_2015, output = Round2Mixtures_2015_Estimates$Output, 
              groups = groups10, colors = colors10, 
              header = Round2Mixtures_2015_Header, set.mfrow = c(5, 2), thin = 10)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Round 2 Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Barplots
QuickBarplot(mixvec = Round2Mixtures_2015, estimatesstats = Round2Mixtures_2015_EstimatesStats, groups = groups10, groups2rows = groups10tworows, header = Round2Mixtures_2015_Header)

## Make violin plots of posteriors with RGs sorted
ViolinPlot(estimates = Round2Mixtures_2015_Estimates, groups = groups10tworows, colors = colors10, header = Round2Mixtures_2015_Header)
rm(Round2Mixtures_2015_Estimates)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Create Final Estimates Objects for Each Temporal Strata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Round1Mixtures_2015_EstimatesStats <- dget("Estimates objects/Round1Mixtures_2015_EstimatesStats.txt")
Round2Mixtures_2015_EstimatesStats <- dget("Estimates objects/Round2Mixtures_2015_EstimatesStats.txt")

# View as tables by year
require(reshape)
samp.df.2015 <- data.frame(t(sapply(KMA2015Strata, function(strata) {
  location <- unlist(strsplit(x = strata, split = "_"))[1]
  temporal.strata <- as.numeric(unlist(strsplit(x = strata, split = "_"))[2])
  n <- get(paste(strata, ".gcl", sep = ""))$n
  c(location = location, temporal.strata = temporal.strata, n = n)
})))
cast(data = samp.df.2015, location ~ temporal.strata)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1_Early
str(Round1Mixtures_2015_EstimatesStats)
KMA2015Strata_1_Early_EstimatesStats <- Round1Mixtures_2015_EstimatesStats[1:3]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2_Late
str(Round2Mixtures_2015_EstimatesStats)
KMA2015Strata_2_Late_EstimatesStats <- c(Round2Mixtures_2015_EstimatesStats,
                                         Round1Mixtures_2015_EstimatesStats[4])

str(KMA2015Strata_1_Early_EstimatesStats)
str(KMA2015Strata_2_Late_EstimatesStats)


dput(x = KMA2015Strata_1_Early_EstimatesStats, file = "Estimates objects/Final/KMA2015Strata_1_Early_EstimatesStats.txt")
dput(x = KMA2015Strata_2_Late_EstimatesStats, file = "Estimates objects/Final/KMA2015Strata_2_Late_EstimatesStats.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot stock composition results 2015 KMA Mixtures ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
KMA2015Strata_1_Early_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2015Strata_1_Early_EstimatesStats.txt")
KMA2015Strata_2_Late_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2015Strata_2_Late_EstimatesStats.txt")

KMA2015Strata_EstimatesStats <- c(KMA2015Strata_1_Early_EstimatesStats, 
                                  KMA2015Strata_2_Late_EstimatesStats)
dput(x = KMA2015Strata_EstimatesStats, file = "Estimates objects/Final/KMA2015Strata_EstimatesStats.txt")
KMA2015Strata_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2015Strata_EstimatesStats.txt")

str(KMA2015Strata_EstimatesStats)


KMA2015 <- c("KALITC15", "KEASTC15", "KWESTC15", "KMAINC15")
dput(KMA2015, file = "Objects/KMA2015.txt")
TempMix15 <- sapply(KMA2015, function(geo) {grep(pattern = geo, x = names(KMA2015Strata_EstimatesStats), value = TRUE)} )

Legend15 <- setNames(object = c("June 1-July 5", "July 6-August 5"), 
                     nm = c("1_Early", "2_Late"))
TempLegend15 <- sapply(KMA2015, function(geo) {
  Legend15[sapply(TempMix15[[geo]], function(strata) {unlist(strsplit(x = strata, split = paste(geo, "_", sep = '')))[2]} )]
} )

GeoHeader <- setNames(object = c(paste0("SW Kodiak/Alitak 255, 256, 257"),
                                 paste0("Eastside Kodiak/Afognak 252, 258, 259"),
                                 paste0("Westside Kodiak/Afognak 251, 253, 254"),
                                 paste0("Mainland 262")),
                      nm = unlist(strsplit(x = KMA2015, split = "15")))

ProportionColors <- colorpanel(n = 2, low = "blue", high = "white")
TempProportionColors15 <- sapply(KMA2015, function(geo) {
  ProportionColors[sapply(TempMix15[[geo]], function(strata) {as.numeric(unlist(strsplit(x = strata, split = "_"))[2])} )]
} )

Estimates <- KMA2015Strata_EstimatesStats
Groups <- groups10
Groups2Rows <- groups10tworows
cex.lab <- 1.5
cex.yaxis <- 1.3
cex.xaxis <- 0.6
cex.main <- 1.7
cex.leg <- 1.3
ci.lwd <- 2.5

dir.create("Figures/2015")
require(devEMF)

sapply(names(TempMix15), function(geomix) {
  emf(file = paste("Figures/2015/", geomix, ".emf", sep = ''), width = 8.5, height = 6.5, family = "sans", bg = "white")
  par(mar = c(2.1, 4.1, 2.6, 0.6))
  
  Barplot1 <- barplot2(height = t(sapply(TempMix15[[geomix]], function(tempmix) {Estimates[[tempmix]][, "median"]})) * 100, 
                       beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                       ci.l = t(sapply(TempMix15[[geomix]], function(tempmix) {Estimates[[tempmix]][, "5%"]})) * 100, 
                       ci.u = t(sapply(TempMix15[[geomix]], function(tempmix) {Estimates[[tempmix]][, "95%"]})) * 100, 
                       ylim = c(0, 100), col = TempProportionColors15[[geomix]], cex.axis = cex.yaxis, yaxt = "n", xaxt = 'n')
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
  legend(legend = TempLegend15[[geomix]], x = "topleft", fill = TempProportionColors15[[geomix]], border = "black", bty = "n", cex = cex.leg, title="2015")
  abline(h = 0, xpd = FALSE)
  
  mtext(text = "Percentage of Catch", side = 2, cex = cex.yaxis, line = 3)
  mtext(text = Groups2Rows, side = 1, line = 1, at = apply(Barplot1, 2, mean), adj = 0.5, cex = cex.xaxis)
  mtext(text = GeoHeader[unlist(strsplit(geomix, split = "15"))], side = 3, cex = cex.main, line = 1)
  dev.off()
})



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 2015 Harvest Data ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
harvest15 <- read.csv(file = "Harvest/Salmon Catch by Day and Stat Area 2015 Date.csv", as.is = TRUE)
str(harvest15)

# Convert Date
harvest15$Date.Landed <- as.Date(harvest15$Date.Landed, format = "%Y-%m-%d")
harvest15$Date.Fishing.Began <- as.Date(harvest15$Date.Fishing.Began, format = "%Y-%m-%d")
harvest15$Date.Fishing.Ended <- as.Date(harvest15$Date.Fishing.Ended, format = "%Y-%m-%d")

# Create Temporal Strata
harvest15$Strata <- ifelse(harvest15$Date.Fishing.Began >= as.Date("2015-06-01") & harvest15$Date.Fishing.Began <= as.Date("2015-07-05"),
                           "Early",
                           ifelse(harvest15$Date.Fishing.Began >= as.Date("2015-07-06") & harvest15$Date.Fishing.Began <= as.Date("2015-08-05"),
                                  "Late",
                                  NA))
table(harvest15$Strata)

# Create Geographic Strata
unique(harvest15$Stat.Area)

Stat.Area.Eastside <- unique(harvest15$Stat.Area)[c(
  which(unique(harvest15$Stat.Area) >= 25800 & unique(harvest15$Stat.Area) <=25999),
  which(unique(harvest15$Stat.Area) >= 25200 & unique(harvest15$Stat.Area) <=25299))]

Stat.Area.Westside <- unique(harvest15$Stat.Area)[c(
  which(unique(harvest15$Stat.Area) >= 25100 & unique(harvest15$Stat.Area) <=25199),
  which(unique(harvest15$Stat.Area) >= 25300 & unique(harvest15$Stat.Area) <=25499))]

Stat.Area.SWAlitak <- unique(harvest15$Stat.Area)[c(
  which(unique(harvest15$Stat.Area) >= 25500 & unique(harvest15$Stat.Area) <=25799))]

Stat.Area.Mainland <- unique(harvest15$Stat.Area)[c(
  which(unique(harvest15$Stat.Area) >= 26200))]


length(c(Stat.Area.Eastside, Stat.Area.Westside, Stat.Area.SWAlitak, Stat.Area.Mainland)); length(unique(harvest15$Stat.Area))
unique(harvest15$Stat.Area)[!unique(harvest15$Stat.Area) %in% c(Stat.Area.Eastside, Stat.Area.Westside, Stat.Area.SWAlitak, Stat.Area.Mainland)]


harvest15$Geo <- ifelse(harvest15$Stat.Area %in% Stat.Area.Mainland, "Mainland",
                        ifelse(harvest15$Stat.Area %in% Stat.Area.SWAlitak, "SWAlitak",
                               ifelse(harvest15$Stat.Area %in% Stat.Area.Eastside, "Eastside",
                                      ifelse(harvest15$Stat.Area %in% Stat.Area.Westside, "Westside", NA
                                      ))))

table(harvest15$Strata, harvest15$Geo)

require(plyr)
ddply(.data = harvest15, ~Geo+Strata, summarise, harvest = sum(Number))
daply(.data = harvest15, ~Geo+Strata, summarise, harvest = sum(Number))

aggregate(harvest15$Number, by = list(harvest15$Strata, harvest15$Geo), FUN = sum)
aggregate(Number ~ Geo + Strata, data = harvest15, sum)
cast(aggregate(Number ~ Geo + Strata, data = harvest15, sum), Geo ~ Strata, value = "Number")

md <- melt(data = harvest15, id.vars = c("Strata", "Geo"), measure.vars = "Number", na.rm = TRUE)
cast(md, Geo~Strata, sum)


HarvestByStrata2015 <- data.matrix(cast(aggregate(Number ~ Geo + Strata, data = harvest15, sum), Geo ~ Strata, value = "Number")[c(3, 1, 4, 2), 2:3])
dimnames(HarvestByStrata2015) <- list(KMA2015, c("1_Early", "2_Late"))
HarvestByStrata2015
dput(x = HarvestByStrata2015, file = "Objects/HarvestByStrata2015.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot stock specific harvest results 2015 KMA Mixtures ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Multiply by harvest
HarvestByStrata2015


KMA2015Strata_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2015Strata_EstimatesStats.txt")
str(KMA2015Strata_EstimatesStats)
names(KMA2015Strata_EstimatesStats)
dimnames(KMA2015Strata_EstimatesStats[[1]])


KMA2015Strata_HarvestEstimatesStats <- sapply(names(KMA2015Strata_EstimatesStats)[-8], function(strata) {
  strata.split <- unlist(strsplit(x = strata, split = "_"))
  strata.split <- c(strata.split[1], paste(c(strata.split[2], strata.split[3]), collapse = "_"))
  
  cbind(KMA2015Strata_EstimatesStats[[strata]][, c("mean", "sd", "median", "5%", "95%")] * HarvestByStrata2015[strata.split[1], strata.split[2]],
        KMA2015Strata_EstimatesStats[[strata]][, c("P=0", "GR")])
}, simplify = FALSE )

dput(x = KMA2015Strata_HarvestEstimatesStats, file = "Estimates objects/Final/KMA2015Strata_HarvestEstimatesStats.txt")
KMA2015Strata_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2015Strata_HarvestEstimatesStats.txt")
str(KMA2015Strata_HarvestEstimatesStats)

# What should ymax be?
max(sapply(KMA2015Strata_HarvestEstimatesStats, function(strata) strata[, "95%"]))


TempMix15 <- sapply(KMA2015, function(geo) {grep(pattern = geo, x = names(KMA2015Strata_EstimatesStats), value = TRUE)} )

Legend15 <- setNames(object = c("June 1-July 5", "July 6-August 5"), 
                     nm = c("1_Early", "2_Late"))
TempLegend15 <- sapply(KMA2015, function(geo) {
  Legend15[sapply(TempMix15[[geo]], function(strata) {unlist(strsplit(x = strata, split = paste(geo, "_", sep = '')))[2]} )]
} )

GeoHeader <- setNames(object = c(paste0("SW Kodiak/Alitak 255, 256, 257"),
                                 paste0("Eastside Kodiak/Afognak 252, 258, 259"),
                                 paste0("Westside Kodiak/Afognak 251, 253, 254"),
                                 paste0("Mainland 262")),
                      nm = unlist(strsplit(x = KMA2015, split = "15")))


HarvestColors <- colorpanel(n = 2, low = "green", high = "white")
TempHarvestColors15 <- sapply(KMA2015, function(geo) {
  HarvestColors[sapply(TempMix15[[geo]], function(strata) {as.numeric(unlist(strsplit(x = strata, split = "_"))[2])} )]
} )

Estimates <- KMA2015Strata_HarvestEstimatesStats
Groups <- groups10
Groups2Rows <- groups10tworows
cex.lab <- 1.5
cex.yaxis <- 1.3
cex.xaxis <- 0.6
cex.main <- 1.7
cex.leg <- 1.3
ci.lwd <- 2.5
ymax <- 1500

sapply(names(TempMix15), function(geomix) {
  emf(file = paste("Figures/2015/", geomix, "Harvest.emf", sep = ''), width = 8.5, height = 6.5, family = "sans", bg = "white")
  par(mar = c(2.1, 4.1, 2.6, 0.6))
  
  Barplot1 <- barplot2(height = t(sapply(TempMix15[[geomix]], function(tempmix) {Estimates[[tempmix]][, "median"]})), 
                       beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                       ci.l = t(sapply(TempMix15[[geomix]], function(tempmix) {Estimates[[tempmix]][, "5%"]})), 
                       ci.u = t(sapply(TempMix15[[geomix]], function(tempmix) {Estimates[[tempmix]][, "95%"]})), 
                       ylim = c(0, ymax), col = TempHarvestColors15[[geomix]], cex.axis = cex.yaxis, yaxt = "n", xaxt = 'n')
  axis(side = 2, at = seq(0, ymax, 500), labels = formatC(x = seq(0, ymax, 500), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
  legend(legend = TempLegend15[[geomix]], x = "topleft", fill = TempHarvestColors15[[geomix]], border = "black", bty = "n", cex = cex.leg, title="2015")
  abline(h = 0, xpd = FALSE)
  
  mtext(text = "Number of Fish Harvested", side = 2, cex = cex.yaxis, line = 3)
  mtext(text = Groups2Rows, side = 1, line = 1, at = apply(Barplot1, 2, mean), adj = 0.5, cex = cex.xaxis)
  mtext(text = GeoHeader[unlist(strsplit(geomix, split = "15"))], side = 3, cex = cex.main, line = 1)
  dev.off()
})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Stratified Regional Roll-Ups 2015 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
HarvestByStrata2015


# For this to work properly, the output mixture directories all need to be in the same place (can be moved back afterwards)
sapply(KMA2015[1:3], function(geomix) {
  assign(x = paste(geomix, "_Annual_Stratified", sep = ''), 
         value = StratifiedEstimator.GCL(
           groupvec = seq(groups10), groupnames = groups10, maindir = "BAYES/Output", 
           mixvec = paste(geomix, colnames(HarvestByStrata2015)[!is.na(HarvestByStrata2015[geomix,])], sep = "_"), 
           catchvec = HarvestByStrata2015[geomix, !is.na(HarvestByStrata2015[geomix,])], 
           newname = paste(geomix, "_Annual_Stratified", sep = ''), priorname = '',
           ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1), 
         pos = 1)
  dput(x = get(paste(geomix, "_Annual_Stratified", sep = '')), file = paste("Estimates objects/", geomix, "_Annual_Stratified.txt", sep = ''))
} ); beep(5)

str(KALITC15_Annual_Stratified)

Round1Mixtures_2015_Estimates <- dget(file = "Estimates objects/Round1Mixtures_2015_Estimates.txt")
KMAINC15_Annual_Stratified <- list(Stats = Round1Mixtures_2015_Estimates$Stats$KMAINC15_2_Late,
                                   Output = Round1Mixtures_2015_Estimates$Output$KMAINC15_2_Late)
dput(x = KMAINC15_Annual_Stratified, file = "Estimates objects/KMAINC15_Annual_Stratified.txt")

# Create a list object with all Stratified Annual Rollups for the "Estimates objects/Final" folder
KMA2015_Annual_EstimatesStats <- sapply(KMA2015, function(geomix) {
  Stats <- get(paste(geomix, "_Annual_Stratified", sep = ''))$Stats
  Stats
}, simplify = FALSE)
str(KMA2015_Annual_EstimatesStats)
dput(x = KMA2015_Annual_EstimatesStats, file = "Estimates objects/Final/KMA2015_Annual_EstimatesStats.txt")




# Create a list object with all Stratified Annual Harvest
KMA2015_Annual_HarvestEstimatesStats <- sapply(KMA2015, function(strata) {
  cbind(KMA2015_Annual_EstimatesStats[[strata]][, c("mean", "sd", "median", "5%", "95%")] * sum(HarvestByStrata2015[strata, ], na.rm = TRUE),
        KMA2015_Annual_EstimatesStats[[strata]][, c("P=0", "GR")])
}, simplify = FALSE )

dput(x = KMA2015_Annual_HarvestEstimatesStats, file = "Estimates objects/Final/KMA2015_Annual_HarvestEstimatesStats.txt")
str(KMA2015_Annual_HarvestEstimatesStats)




# Create a matrix of annual means
Annual2015_Stratified_Estimates <- sapply(KMA2015, function(geomix) {
  KMA2015_Annual_EstimatesStats[[geomix]][, "mean"]
})

# Create a matrix of early strata means
KMA2015Strata_1_Early_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2015Strata_1_Early_EstimatesStats.txt")
EarlyStrata2015_Estimates <- sapply(KMA2015Strata_1_Early_EstimatesStats, function(geomix) {
  geomix[, "mean"]
})
colnames(EarlyStrata2015_Estimates) <- KMA2015[-4]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot Annual vs. Early Means
par(mar = c(4.1, 5.1, 3.1, 1.1))
sapply(KMA2015[-4], function(geomix) {
  Barplot <- barplot2(height = rbind(Annual2015_Stratified_Estimates[, geomix], EarlyStrata2015_Estimates[, geomix]) * 100, 
                      beside = TRUE, col = c("blue", "white"), yaxt = 'n', xaxt = 'n', main = geomix, ylab = "Precent of Mixture",
                      cex.lab = 2, cex.main = 2, ylim = c(0, 100))
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = ",", digits = 0, format = "f"), cex.axis = 1.5)
  abline(h = 0, xpd = FALSE)
  text(x = colMeans(Barplot), y = -1, labels = groups10tworows, srt = 90, adj = 1, xpd = TRUE, cex = 0.6)
  legend("topleft", legend = c("Annual Means", "Early Strata Means"), fill = c("blue", "white"), bty = 'n', cex = 1.5)
})



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 2016 Harvest Data ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
harvest16 <- read.csv(file = "Harvest/Salmon Catch by Day and Stat Area 2016 Date.csv", as.is = TRUE)
str(harvest16)

# Convert Date
harvest16$Date.Landed <- as.Date(harvest16$Date.Landed, format = "%Y-%m-%d")
harvest16$Date.Fishing.Began <- as.Date(harvest16$Date.Fishing.Began, format = "%Y-%m-%d")
harvest16$Date.Fishing.Ended <- as.Date(harvest16$Date.Fishing.Ended, format = "%Y-%m-%d")

# Create Temporal Strata
harvest16$Strata <- ifelse(harvest16$Date.Landed >= as.Date("2016-06-01") & harvest16$Date.Landed <= as.Date("2016-07-05"),
                           "Early",
                           ifelse(harvest16$Date.Landed >= as.Date("2016-07-06") & harvest16$Date.Landed <= as.Date("2016-08-05"),
                                  "Late",
                                  NA))
table(harvest16$Strata)

# Create Geographic Strata
unique(harvest16$Stat.Area)

Stat.Area.Eastside <- unique(harvest06$Stat.Area)[c(
  which(unique(harvest06$Stat.Area) >= 25800 & unique(harvest06$Stat.Area) <=25999),
  which(unique(harvest06$Stat.Area) >= 25200 & unique(harvest06$Stat.Area) <=25299))]

Stat.Area.Westside <- unique(harvest06$Stat.Area)[c(
  which(unique(harvest06$Stat.Area) >= 25100 & unique(harvest06$Stat.Area) <=25199),
  which(unique(harvest06$Stat.Area) >= 25300 & unique(harvest06$Stat.Area) <=25499))]

Stat.Area.SWAlitak <- unique(harvest06$Stat.Area)[c(
  which(unique(harvest06$Stat.Area) >= 25500 & unique(harvest06$Stat.Area) <=25799))]

Stat.Area.Mainland <- unique(harvest06$Stat.Area)[c(
  which(unique(harvest06$Stat.Area) >= 26200))]


length(c(Stat.Area.Eastside, Stat.Area.Westside, Stat.Area.SWAlitak, Stat.Area.Mainland)); length(unique(harvest16$Stat.Area))
unique(harvest16$Stat.Area)[!unique(harvest16$Stat.Area) %in% c(Stat.Area.Eastside, Stat.Area.Westside, Stat.Area.SWAlitak, Stat.Area.Mainland)]


harvest16$Geo <- ifelse(harvest16$Stat.Area %in% Stat.Area.Mainland, "Mainland",
                        ifelse(harvest16$Stat.Area %in% Stat.Area.SWAlitak, "SWAlitak",
                               ifelse(harvest16$Stat.Area %in% Stat.Area.Eastside, "Eastside",
                                      ifelse(harvest16$Stat.Area %in% Stat.Area.Westside, "Westside", NA
                                      ))))

table(harvest16$Strata, harvest16$Geo)

require(plyr)
ddply(.data = harvest16, ~Geo+Strata, summarise, harvest = sum(Number))
daply(.data = harvest16, ~Geo+Strata, summarise, harvest = sum(Number))

aggregate(harvest16$Number, by = list(harvest16$Strata, harvest16$Geo), FUN = sum)


