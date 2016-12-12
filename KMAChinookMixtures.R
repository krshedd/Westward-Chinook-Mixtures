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
KMAobjects <- KMAobjects[!KMAobjects %in% c("OriginalLocusControl48.txt", "OLD")]
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


invisible(sapply(c("KMA2014", "KMA2015"), function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/PoolingByYear/", silly, ".txt", sep = "")), pos = 1)} )); beep(2)


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
KMAobjects <- KMAobjects[!KMAobjects %in% c("OriginalLocusControl48.txt", "OLD")]
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
groups10tworows <- c("Russia\n", "Eastern\nBering", "North\nPen", "Chignik\n", "Kodiak\n", "Cook\nInlet", "Copper\n", "SEAK\n", "British\nColumbia", "West\nCoast US")
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
  emf(file = paste("Figures/2014/", geomix, ".emf", sep = ''), width = 8.5, height = 6.5, family = "serif", bg = "white")
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




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Size Parameters
Groups <- groups10
Groups2Rows <- groups10tworows
cex.lab <- 1.5
cex.xaxis <- 0.5
cex.yaxis <- 1.3
cex.leg <- 1.1
ci.lwd <- 2.5
geomix = "KSPENCHIG14"


emf(file = paste("Figures/2014/", geomix, ".emf", sep = ''), width = 6, height = 5.75, family = "serif", bg = "white")

# Three barplot layout
layoutmat <- matrix(data=c(  1, 2,
                             3, 4), nrow = 2, ncol = 2, byrow = TRUE)

layout(mat = layoutmat, widths = c(0.075, 1, 1), heights = c(1, 0.1))
par(mar = rep(0, 4))
par(family = "times")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Y-axis label
plot.new()
text(x = 0.25, y = 0.5, labels = "Percentage of Catch", srt = 90, cex = cex.lab)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2014 Barplot
Estimates <- KMA2014Strata_EstimatesStats
geomix = "KSPENCHIG14"
par(mar = c(1, 1, 1, 1))
Barplot14 <- barplot2(height = t(Estimates[[geomix]][, "median"]) * 100, 
                      beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                      ci.l = t(Estimates[[geomix]][, "5%"]) * 100, 
                      ci.u = t(Estimates[[geomix]][, "95%"]) * 100, 
                      ylim = c(0, 100), col = "blue", yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
legend(legend = "June 1-August 5", x = "topleft", fill = "blue", border = "black", bty = "n", cex = cex.leg, title="2014")
abline(h = 0, xpd = FALSE)

mtext(text = Groups2Rows, side = 1, line = 1, at = apply(Barplot14, 2, mean), adj = 0.5, cex = cex.xaxis)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Blank Corner
par(mar = rep(0, 4))
plot.new()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## x-axis label
par(mar = rep(0, 4))
plot.new()
text(x = 0.5, y = 0.25, labels = "Reporting Group", cex = cex.lab)


dev.off()






Estimates <- KMA2014Strata_EstimatesStats
Groups <- groups10
Groups2Rows <- groups10tworows
cex.lab <- 1.5
cex.yaxis <- 1.3
cex.xaxis <- 0.6
cex.main <- 1.7
cex.leg <- 1.3
ci.lwd <- 2.5

geomix = "KSPENCHIG14"

emf(file = paste("Figures/2014/", geomix, ".emf", sep = ''), width = 6, height = 5.75, family = "times", bg = "white")
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
mtext(text = paste0("South Peninsula / Chignik Outside 282", "\u2013", "285; 272, 273, 275"), side = 3, cex = cex.main, line = 1)
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
  emf(file = paste("Figures/2014/", geomix, "Harvest.emf", sep = ''), width = 8.5, height = 6.5, family = "serif", bg = "white")
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
  emf(file = paste("Figures/2015/", geomix, ".emf", sep = ''), width = 8.5, height = 6.5, family = "serif", bg = "white")
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
require(reshape)
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
  emf(file = paste("Figures/2015/", geomix, "Harvest.emf", sep = ''), width = 8.5, height = 6.5, family = "serif", bg = "white")
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
KMAobjects <- KMAobjects[!KMAobjects %in% c("OriginalLocusControl48.txt", "OLD")]
KMAobjects

invisible(sapply(KMAobjects, function(objct) {assign(x = unlist(strsplit(x = objct, split = ".txt")), value = dget(file = paste(getwd(), "Objects", objct, sep = "/")), pos = 1) })); beep(2)


## Get un-altered 2014-2015 mixtures
invisible(sapply(KMA2014_2015Strata, function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/OriginalCollections_Strata_PostQC/", silly, ".txt", sep = "")), pos = 1)} )); beep(2)
objects(pattern = "\\.gcl")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Pull genotypes from LOKI 2016 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Pull all data for each silly code and create .gcl objects for each
KMA2016 <- c("KKMAC16")
dput(x = KMA2016, file = "Objects/KMA2016.txt")

LOKI2R.GCL(sillyvec = c(KMA2016), username = username, password = password)

rm(username, password)
objects(pattern = "\\.gcl")

## Save unaltered .gcl's as back-up:
invisible(sapply(c(KMA2016), function(silly) {dput(x = get(paste(silly, ".gcl", sep = '')), file = paste("Raw genotypes/OriginalCollections/" , silly, ".txt", sep = ''))} )); beep(8)

## Original sample sizes by SILLY
collection.size.original <- sapply(c(KMA2016), function(silly) get(paste(silly, ".gcl", sep = ""))$n)
# matrix(data = collection.size.original, ncol = 2, dimnames = list(c("Alitak", "Ayakulik", "Karluk", "Uganik", "Uyak"), 2014:2015))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Define strata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## All fish have a capture date?
sapply(c(KMA2016), function(silly) {sum(is.na(get(paste(silly, ".gcl", sep = ''))$attributes$CAPTURE_DATE))} )  # Zeros are good

## Confirming samples sizes by date
sapply(c(KMA2016), function(silly) {table(get(paste(silly, ".gcl", sep = ''))$attributes$CAPTURE_DATE)} )

str(KKMAC16.gcl)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Pooling roadmap
# 1) Pool geographically within year by "CAPTURE_LOCATION"
# 2) Pool temporally within geographic by year by "CAPTURE_DATE"

#~~~~~~~~~~~~~~~~~~
# Dput these sillys for now
dput(x = KKMAC16.gcl, file = "Raw genotypes/PoolingByYear/KMA2016.txt")
KMA2016.gcl <- KKMAC16.gcl

KMA2016.gcl <- dget(file = "Raw genotypes/PoolingByYear/KMA2016.txt")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Spatio temporal strata for 2016 Kodiak samples

silly = "KMA2016"
loci = loci48
geostrata = cbind(object = unique(KMA2016.gcl$attributes$CAPTURE_LOCATION), nm = c("KEASTC16", "KMAINC16", "KWESTC16", "KALITC16"))
tempstrata = cbind(object = unique(KMA2016.gcl$attributes$MESH_SIZE_COMMENT), nm = c("1_Early", "2_Late"))
geostrata; tempstrata

table(KMA2016.gcl$attributes$CAPTURE_LOCATION, KMA2016.gcl$attributes$MESH_SIZE_COMMENT)

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

grep(pattern = ".gcl", x = objects(pattern = "16_"), value = TRUE)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Mixture sillyvec
KMA2016Strata <- unlist(strsplit(x = grep(pattern = "16_", x = objects(pattern = "\\.gcl"), value = TRUE), split = "\\.gcl"))
KMA2014_2016Strata <- c("KSPENCHIG14", KMA2014Strata, KMA2015Strata, KMA2016Strata)

dput(x = KMA2016Strata, file = "Objects/KMA2016Strata.txt")
dput(x = KMA2014_2016Strata, file = "Objects/KMA2014_2016Strata.txt")

# Confirm sample sizes
sapply(KMA2016Strata, function(silly) get(paste(silly, ".gcl", sep = ""))$n)

# View as tables by year
require(reshape)
samp.df.2016 <- data.frame(t(sapply(KMA2016Strata, function(strata) {
  location <- unlist(strsplit(x = strata, split = "_"))[1]
  temporal.strata <- as.numeric(unlist(strsplit(x = strata, split = "_"))[2])
  n <- get(paste(strata, ".gcl", sep = ""))$n
  c(location = location, temporal.strata = temporal.strata, n = n)
})))
cast(data = samp.df.2016, location ~ temporal.strata)



# dput mixture sillys
invisible(sapply(KMA2016Strata, function(silly) {dput(x = get(paste(silly, ".gcl", sep = '')), file = paste("Raw genotypes/OriginalCollections_Strata/" , silly, ".txt", sep = ''))} )); beep(8)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Data QC/Massage ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Confirm appropriate strata
sapply(KMA2016Strata, function(mix) {table(get(paste(mix, ".gcl", sep = ''))$attributes$CAPTURE_LOCATION,
                                           get(paste(mix, ".gcl", sep = ''))$attributes$MESH_SIZE_COMMENT)})

# Begin QC
require(xlsx)

KMA2016Strata

KMA2016Strata_SampleSizes <- matrix(data = NA, nrow = length(KMA2016Strata), ncol = 4, 
                                    dimnames = list(KMA2016Strata, c("Genotyped", "Missing", "Duplicate", "Final")))

#### Check loci
## Get sample size by locus
Original_KMA2016Strata_SampleSizebyLocus <- SampSizeByLocus.GCL(sillyvec = KMA2016Strata, loci = loci48)
min(Original_KMA2016Strata_SampleSizebyLocus)  ## 268
apply(Original_KMA2016Strata_SampleSizebyLocus, 1, min) / apply(Original_KMA2016Strata_SampleSizebyLocus, 1, max)  # Good, all > 0.9

Original_KMA2016Strata_PercentbyLocus <- apply(Original_KMA2016Strata_SampleSizebyLocus, 1, function(row) {row / max(row)} )
which(apply(Original_KMA2016Strata_PercentbyLocus, 2, min) < 0.8)  # no re-runs!

require(lattice)
new.colors <- colorRampPalette(c("black", "white"))
levelplot(t(Original_KMA2016Strata_PercentbyLocus), 
          col.regions = new.colors, 
          at = seq(from = 0, to = 1, length.out = 100), 
          main = "% Genotyped", xlab = "SILLY", ylab = "Locus", 
          scales = list(x = list(rot = 90)), 
          aspect = "fill")  # aspect = "iso" will make squares


#### Check individuals
## View Histogram of Failure Rate by Strata
invisible(sapply(KMA2016Strata, function(mix) {
  my.gcl <- get(paste(mix, ".gcl", sep = ''))
  failure <- apply(my.gcl$scores[, , 1], 1, function(ind) {sum(ind == "0") / length(ind)} )
  hist(x = failure, main = mix, xlab = "Failure Rate", col = 8, xlim = c(0, 1), ylim = c(0, 20), breaks = seq(from = 0, to = 1, by = 0.02))
  abline(v = 0.2, lwd = 3)
}))

## Experimental method to look for a way to remove "missing" loci fish based on a strata specific failure rate (rather than 0.8 globally, as done normally)
# Experimented with loess, smooth.split, polynomial regressions, outlier (Franz), and Cooks Distance methods to model outliers.
# Thought of a clever idea to look at the "shoulder" of the cumulative success rate, but this concept is a bit weird.
# Planning to keep it simple!
invisible(sapply(KMA2016Strata, function(mix) {
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
Original_KMA2016Strata_ColSize <- sapply(paste(KMA2016Strata, ".gcl", sep = ''), function(x) get(x)$n)
KMA2016Strata_SampleSizes[, "Genotyped"] <- Original_KMA2016Strata_ColSize


### Missing
## Remove individuals with >20% missing data
KMA2016Strata_MissLoci <- RemoveIndMissLoci.GCL(sillyvec = KMA2016Strata, proportion = 0.8)

## Get number of individuals per silly after removing missing loci individuals
ColSize_KMA2016Strata_PostMissLoci <- sapply(paste(KMA2016Strata, ".gcl", sep = ''), function(x) get(x)$n)
KMA2016Strata_SampleSizes[, "Missing"] <- Original_KMA2016Strata_ColSize-ColSize_KMA2016Strata_PostMissLoci


### Duplicate
## Check within collections for duplicate individuals.
KMA2016Strata_DuplicateCheck95MinProportion <- CheckDupWithinSilly.GCL(sillyvec = KMA2016Strata, loci = loci48, quantile = NULL, minproportion = 0.95)
KMA2016Strata_DuplicateCheckReportSummary <- sapply(KMA2016Strata, function(x) KMA2016Strata_DuplicateCheck95MinProportion[[x]]$report)
KMA2016Strata_DuplicateCheckReportSummary

## Remove duplicate individuals
KMA2016Strata_RemovedDups <- RemoveDups.GCL(KMA2016Strata_DuplicateCheck95MinProportion)

## Get number of individuals per silly after removing duplicate individuals
ColSize_KMA2016Strata_PostDuplicate <- sapply(paste(KMA2016Strata, ".gcl", sep = ''), function(x) get(x)$n)
KMA2016Strata_SampleSizes[, "Duplicate"] <- ColSize_KMA2016Strata_PostMissLoci-ColSize_KMA2016Strata_PostDuplicate


### Final
KMA2016Strata_SampleSizes[, "Final"] <- ColSize_KMA2016Strata_PostDuplicate
KMA2016Strata_SampleSizes


write.xlsx(KMA2016Strata_SampleSizes, file = "Output/KMA2016Strata_SampleSizes.xlsx")
dput(x = KMA2016Strata_SampleSizes, file = "Objects/KMA2016Strata_SampleSizes.txt")


### Create a final sample size matrix for 2014-2016
KMA2014_2016Strata_SampleSizes <- rbind(KMA2014_2015Strata_SampleSizes, KMA2016Strata_SampleSizes)
dput(x = KMA2014_2016Strata_SampleSizes, file = "Objects/KMA2014_2016Strata_SampleSizes.txt")
write.xlsx(KMA2014_2016Strata_SampleSizes, file = "Output/KMA2014_2016Strata_SampleSizes.xlsx")



# dput postQC mixture sillys
invisible(sapply(KMA2016Strata, function(silly) {dput(x = get(paste(silly, ".gcl", sep = '')), file = paste("Raw genotypes/OriginalCollections_Strata_PostQC/" , silly, ".txt", sep = ''))} )); beep(8)



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
KMAobjects <- KMAobjects[!KMAobjects %in% c("OriginalLocusControl48.txt", "OLD")]
KMAobjects

invisible(sapply(KMAobjects, function(objct) {assign(x = unlist(strsplit(x = objct, split = ".txt")), value = dget(file = paste(getwd(), "Objects", objct, sep = "/")), pos = 1) })); beep(2)


## Get un-altered mixtures
invisible(sapply(KMA2014_2016Strata, function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/OriginalCollections_Strata_PostQC/", silly, ".txt", sep = "")), pos = 1)} )); beep(2)
objects(pattern = "\\.gcl")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 1 MSA files for BAYES 2016 Early Strata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require(reshape)
samp.df.2016 <- data.frame(t(sapply(KMA2016Strata, function(strata) {
  location <- unlist(strsplit(x = strata, split = "_"))[1]
  temporal.strata <- as.numeric(unlist(strsplit(x = strata, split = "_"))[2])
  n <- get(paste(strata, ".gcl", sep = ""))$n
  c(location = location, temporal.strata = temporal.strata, n = n)
})))
cast(data = samp.df.2016, location ~ temporal.strata)

# Indetify Strata to Run
KMA2016Strata_1_Early <- grep(pattern = "1_Early", x = KMA2016Strata, value = TRUE)
Round1Mixtures_2016 <- KMA2016Strata_1_Early
dput(x = Round1Mixtures_2016, file = "Objects/Round1Mixtures_2016.txt")

# Create rolling prior based on 2015 Round 1 estimates
Round1Mixtures_2015_EstimatesStats <- dget(file = "Estimates objects/Round1Mixtures_2015_EstimatesStats.txt")
str(Round1Mixtures_2015_EstimatesStats)

Round1Mixtures_2016_Prior <- sapply(Round1Mixtures_2015_EstimatesStats[1:3], function(Mix) {
  Prior.GCL(groupvec = groupvec10, groupweights = Mix[, "mean"], minval = 0.01)}, simplify = FALSE)
names(Round1Mixtures_2016_Prior) <- gsub(pattern = "C15", replacement = "C16", x = names(Round1Mixtures_2016_Prior))  # This changes the names
str(Round1Mixtures_2016_Prior)
str(KMA211Pops10FlatPrior)

# Add flat prior for Mainland Early as this is the first time we have sampled that spatio-temporal strata
Round1Mixtures_2016_Prior <- c(Round1Mixtures_2016_Prior, list("KMAINC16_1_Early" = KMA211Pops10FlatPrior))
str(Round1Mixtures_2016_Prior)


dput(x = Round1Mixtures_2016_Prior, file = "Objects/Round1Mixtures_2016_Prior.txt")
str(Round1Mixtures_2016_Prior)

# Verify
sapply(Round1Mixtures_2016, function(geomix) {plot(as.vector(Round1Mixtures_2016_Prior[[geomix]]), type = "h", main = geomix)})

## Dumping Mixture files
sapply(Round1Mixtures_2016, function(Mix) {CreateMixture.GCL(sillys = Mix, loci = loci42, IDs = NULL, mixname = Mix, dir = "BAYES/Mixture", type = "BAYES", PT = FALSE)} )

## Dumping Control files
sapply(Round1Mixtures_2016, function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA211Pops, loci = loci42, mixname = Mix, basename = "KMA211Pops42Loci", suffix = "", nreps = 40000, nchains = 5,
                        groupvec = groupvec10, priorvec = Round1Mixtures_2016_Prior[[Mix]], initmat = KMA211PopsInits, dir = "BAYES/Control",
                        seeds = KMA211PopsChinookSeeds, thin = c(1, 1, 100), mixfortran = KMA21142MixtureFormat, basefortran = KMA211Pops42Loci.baseline, switches = "F T F T T T F")
})

## Create output directory
sapply(Round1Mixtures_2016, function(Mix) {dir.create(paste("BAYES/Output/", Mix, sep = ""))})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Go run BAYES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Round 1 2016 Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round1Mixtures_2016_Estimates <- CustomCombineBAYESOutput.GCL(
  groupvec = seq(groups10), groupnames = groups10 ,
  maindir = "BAYES/Output", 
  mixvec = Round1Mixtures_2016, prior = "",  
  ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)

# Dput 1) estimates stats + posterior output & 2) estimates stats
dput(Round1Mixtures_2016_Estimates, file = "Estimates objects/Round1Mixtures_2016_Estimates.txt")
dput(Round1Mixtures_2016_Estimates$Stats, file = "Estimates objects/Round1Mixtures_2016_EstimatesStats.txt")

Round1Mixtures_2016_Estimates <- dget(file = "Estimates objects/Round1Mixtures_2016_Estimates.txt")
Round1Mixtures_2016_EstimatesStats <- dget(file = "Estimates objects/Round1Mixtures_2016_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2
sapply(Round1Mixtures_2016_Estimates$Stats, function(Mix) {Mix[, "GR"]})
sapply(Round1Mixtures_2016_Estimates$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})
require(gplots)
par(mfrow = c(1, 1), mar = c(3.1, 4.6, 3.1, 1.1), oma = rep(0, 4))
sapply(Round1Mixtures_2016, function(Mix) {
  BarPlot <- barplot2(Round1Mixtures_2016_EstimatesStats[[Mix]][, "GR"], col = "blue", 
                      ylim = c(1, pmax(1.5, max(Round1Mixtures_2016_EstimatesStats[[Mix]][, "GR"]))), ylab = "Gelman-Rubin", 
                      xpd = FALSE, main = Mix, names.arg = '', cex.lab = 1.5, cex.main = 2)
  abline(h = 1.2, lwd = 3, xpd = FALSE)
  text(x = BarPlot, y = 1, labels = groups10tworows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
})

# Quick look at raw posterior output
str(Round1Mixtures_2016_Estimates$Output)
Round1Mixtures_2016_Header <- setNames(object = c("Southwest Kodiak / Alitak Early June 1-July 5, 2016",
                                                  "Eastside Kodiak Early June 1-July 5, 2016",
                                                  "Westwide Kodiak Early June 1-July 5, 2016",
                                                  "Mainland Kodiak Early June 1-July 5, 2016"), 
                                       nm = Round1Mixtures_2016)
dput(x = Round1Mixtures_2016_Header, file = "Objects/Round1Mixtures_2016_Header.txt")

PlotPosterior(mixvec = Round1Mixtures_2016, output = Round1Mixtures_2016_Estimates$Output, 
              groups = groups10tworows, colors = colors10, 
              header = Round1Mixtures_2016_Header, set.mfrow = c(5, 2), thin = 10)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Round 1 2016 Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Quick barplot
QuickBarplot(mixvec = Round1Mixtures_2016, estimatesstats = Round1Mixtures_2016_EstimatesStats, groups = groups10, groups2rows = groups10tworows, header = Round1Mixtures_2016_Header)

## Make violin plots of posteriors with RGs sorted
ViolinPlot(mixvec = Round1Mixtures_2016, estimates = Round1Mixtures_2016_Estimates, groups = groups10tworows, colors = colors10, header = Round1Mixtures_2016_Header)
rm(Round1Mixtures_2016_Estimates)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 2 MSA files for BAYES 2016 Late Strata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Indetify Strata to Run
KMA2016Strata_2_Late <- grep(pattern = "2_Late", x = KMA2016Strata, value = TRUE)
Round2Mixtures_2016 <- KMA2016Strata_2_Late
dput(x = Round2Mixtures_2016, file = "Objects/Round2Mixtures_2016.txt")

# Create rolling prior based on Round 1 estimates
Round2Mixtures_2016_Prior <- sapply(Round1Mixtures_2016_EstimatesStats[1:4], function(Mix) {Prior.GCL(groupvec = groupvec10, groupweights = Mix[, "mean"], minval = 0.01)}, simplify = FALSE)
names(Round2Mixtures_2016_Prior) <- gsub(pattern = "1_Early", replacement = "2_Late", x = names(Round2Mixtures_2016_Prior))  # This changes the names
dput(x = Round2Mixtures_2016_Prior, file = "Objects/Round2Mixtures_2016_Prior.txt")
str(Round2Mixtures_2016_Prior)

## Dumping Mixture files
sapply(Round2Mixtures_2016, function(Mix) {CreateMixture.GCL(sillys = Mix, loci = loci42, IDs = NULL, mixname = Mix, dir = "BAYES/Mixture", type = "BAYES", PT = FALSE)} )

## Dumping Control files
sapply(Round2Mixtures_2016, function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA211Pops, loci = loci42, mixname = Mix, basename = "KMA211Pops42Loci", suffix = "", nreps = 40000, nchains = 5,
                        groupvec = groupvec10, priorvec = Round2Mixtures_2016_Prior[[Mix]], initmat = KMA211PopsInits, dir = "BAYES/Control",
                        seeds = KMA211PopsChinookSeeds, thin = c(1, 1, 100), mixfortran = KMA21142MixtureFormat, basefortran = KMA211Pops42Loci.baseline, switches = "F T F T T T F")
})

## Create output directory
sapply(Round2Mixtures_2016, function(Mix) {dir.create(paste("BAYES/Output/", Mix, sep = ""))})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Go run BAYES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Round 2 Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round2Mixtures_2016_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = seq(groups10), groupnames = groups10, 
                                                              maindir = "BAYES/Output", mixvec = Round2Mixtures_2016, prior = "",  
                                                              ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)

# Dput 1) estimates stats + posterior output & 2) estimates stats
dput(Round2Mixtures_2016_Estimates, file = "Estimates objects/Round2Mixtures_2016_Estimates.txt")
dput(Round2Mixtures_2016_Estimates$Stats, file = "Estimates objects/Round2Mixtures_2016_EstimatesStats.txt")

Round2Mixtures_2016_Estimates <- dget(file = "Estimates objects/Round2Mixtures_2016_Estimates.txt")
Round2Mixtures_2016_EstimatesStats <- dget(file = "Estimates objects/Round2Mixtures_2016_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2
sapply(Round2Mixtures_2016_Estimates$Stats, function(Mix) {Mix[, "GR"]})
sapply(Round2Mixtures_2016_Estimates$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})
require(gplots)
sapply(Round2Mixtures_2016, function(Mix) {
  BarPlot <- barplot2(Round2Mixtures_2016_EstimatesStats[[Mix]][, "GR"], col = "blue", ylim = c(1, pmax(1.5, max(Round2Mixtures_2016_EstimatesStats[[Mix]][, "GR"]))), ylab = "Gelman-Rubin", type = "h", xpd = FALSE, main = Mix, names.arg = '')
  abline(h = 1.2, lwd = 3, xpd = FALSE)
  text(x = BarPlot, y = 1, labels = groups10tworows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
})

# Quick look at raw posterior output
str(Round2Mixtures_2016_Estimates$Output)
Round2Mixtures_2016_Header <- setNames(object = c("Southwest Kodiak / Alitak Late July 6-August 5, 2016",
                                                  "Eastside Kodiak Late July 6-August 5, 2016",
                                                  "Westwide Kodiak Late July 6-August 5, 2016",
                                                  "Mainland Late July 6-August 5, 2016"), 
                                       nm = Round2Mixtures_2016)
dput(x = Round2Mixtures_2016_Header, file = "Objects/Round2Mixtures_2016_Header.txt")


PlotPosterior(mixvec = Round2Mixtures_2016, output = Round2Mixtures_2016_Estimates$Output, 
              groups = groups10, colors = colors10, 
              header = Round2Mixtures_2016_Header, set.mfrow = c(5, 2), thin = 10)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Round 2 Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Barplots
QuickBarplot(mixvec = Round2Mixtures_2016, estimatesstats = Round2Mixtures_2016_EstimatesStats, groups = groups10, groups2rows = groups10tworows, header = Round2Mixtures_2016_Header)

## Make violin plots of posteriors with RGs sorted
ViolinPlot(estimates = Round2Mixtures_2016_Estimates, groups = groups10tworows, colors = colors10, header = Round2Mixtures_2016_Header)
rm(Round2Mixtures_2016_Estimates)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Create Final Estimates Objects for Each Temporal Strata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Round1Mixtures_2016_EstimatesStats <- dget("Estimates objects/Round1Mixtures_2016_EstimatesStats.txt")
Round2Mixtures_2016_EstimatesStats <- dget("Estimates objects/Round2Mixtures_2016_EstimatesStats.txt")

# View as tables by year
require(reshape)
samp.df.2016 <- data.frame(t(sapply(KMA2016Strata, function(strata) {
  location <- unlist(strsplit(x = strata, split = "_"))[1]
  temporal.strata <- as.numeric(unlist(strsplit(x = strata, split = "_"))[2])
  n <- get(paste(strata, ".gcl", sep = ""))$n
  c(location = location, temporal.strata = temporal.strata, n = n)
})))
cast(data = samp.df.2016, location ~ temporal.strata)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1_Early
str(Round1Mixtures_2016_EstimatesStats)
KMA2016Strata_1_Early_EstimatesStats <- Round1Mixtures_2016_EstimatesStats

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2_Late
str(Round2Mixtures_2016_EstimatesStats)
KMA2016Strata_2_Late_EstimatesStats <- Round2Mixtures_2016_EstimatesStats

str(KMA2016Strata_1_Early_EstimatesStats)
str(KMA2016Strata_2_Late_EstimatesStats)


dput(x = KMA2016Strata_1_Early_EstimatesStats, file = "Estimates objects/Final/KMA2016Strata_1_Early_EstimatesStats.txt")
dput(x = KMA2016Strata_2_Late_EstimatesStats, file = "Estimates objects/Final/KMA2016Strata_2_Late_EstimatesStats.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot stock composition results 2016 KMA Mixtures ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
KMA2016Strata_1_Early_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2016Strata_1_Early_EstimatesStats.txt")
KMA2016Strata_2_Late_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2016Strata_2_Late_EstimatesStats.txt")

KMA2016Strata_EstimatesStats <- c(KMA2016Strata_1_Early_EstimatesStats, 
                                  KMA2016Strata_2_Late_EstimatesStats)
dput(x = KMA2016Strata_EstimatesStats, file = "Estimates objects/Final/KMA2016Strata_EstimatesStats.txt")
KMA2016Strata_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2016Strata_EstimatesStats.txt")

str(KMA2016Strata_EstimatesStats)


KMA2016 <- c("KALITC16", "KEASTC16", "KWESTC16", "KMAINC16")
dput(KMA2016, file = "Objects/KMA2016.txt")
TempMix16 <- sapply(KMA2016, function(geo) {grep(pattern = geo, x = names(KMA2016Strata_EstimatesStats), value = TRUE)} , simplify = FALSE)

Legend16 <- setNames(object = c("June 1-July 5", "July 6-August 5"), 
                     nm = c("1_Early", "2_Late"))
TempLegend16 <- sapply(KMA2016, function(geo) {
  Legend16[sapply(TempMix16[[geo]], function(strata) {unlist(strsplit(x = strata, split = paste(geo, "_", sep = '')))[2]} )]
} , simplify = FALSE)

GeoHeader <- setNames(object = c(paste0("SW Kodiak/Alitak 255, 256, 257"),
                                 paste0("Eastside Kodiak/Afognak 252, 258, 259"),
                                 paste0("Westside Kodiak/Afognak 251, 253, 254"),
                                 paste0("Mainland 262")),
                      nm = unlist(strsplit(x = KMA2016, split = "16")))

ProportionColors <- colorpanel(n = 2, low = "blue", high = "white")
TempProportionColors16 <- sapply(KMA2016, function(geo) {
  ProportionColors[sapply(TempMix16[[geo]], function(strata) {as.numeric(unlist(strsplit(x = strata, split = "_"))[2])} )]
}, simplify = FALSE)

Estimates <- KMA2016Strata_EstimatesStats
Groups <- groups10
Groups2Rows <- groups10tworows
cex.lab <- 1.5
cex.yaxis <- 1.3
cex.xaxis <- 0.6
cex.main <- 1.7
cex.leg <- 1.3
ci.lwd <- 2.5

dir.create("Figures/2016")
require(devEMF)

sapply(names(TempMix16), function(geomix) {
  emf(file = paste("Figures/2016/", geomix, ".emf", sep = ''), width = 8.5, height = 6.5, family = "serif", bg = "white")
  par(mar = c(2.1, 4.1, 2.6, 0.6))
  
  Barplot1 <- barplot2(height = t(sapply(TempMix16[[geomix]], function(tempmix) {Estimates[[tempmix]][, "median"]})) * 100, 
                       beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                       ci.l = t(sapply(TempMix16[[geomix]], function(tempmix) {Estimates[[tempmix]][, "5%"]})) * 100, 
                       ci.u = t(sapply(TempMix16[[geomix]], function(tempmix) {Estimates[[tempmix]][, "95%"]})) * 100, 
                       ylim = c(0, 100), col = TempProportionColors16[[geomix]], cex.axis = cex.yaxis, yaxt = "n", xaxt = 'n')
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
  legend(legend = TempLegend16[[geomix]], x = "topleft", fill = TempProportionColors16[[geomix]], border = "black", bty = "n", cex = cex.leg, title="2016")
  abline(h = 0, xpd = FALSE)
  
  mtext(text = "Percentage of Catch", side = 2, cex = cex.yaxis, line = 3)
  mtext(text = Groups2Rows, side = 1, line = 1, at = apply(Barplot1, 2, mean), adj = 0.5, cex = cex.xaxis)
  mtext(text = GeoHeader[unlist(strsplit(geomix, split = "16"))], side = 3, cex = cex.main, line = 1)
  dev.off()
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

Stat.Area.Eastside <- unique(harvest16$Stat.Area)[c(
  which(unique(harvest16$Stat.Area) >= 25800 & unique(harvest16$Stat.Area) <=25999),
  which(unique(harvest16$Stat.Area) >= 25200 & unique(harvest16$Stat.Area) <=25299))]

Stat.Area.Westside <- unique(harvest16$Stat.Area)[c(
  which(unique(harvest16$Stat.Area) >= 25100 & unique(harvest16$Stat.Area) <=25199),
  which(unique(harvest16$Stat.Area) >= 25300 & unique(harvest16$Stat.Area) <=25499))]

Stat.Area.SWAlitak <- unique(harvest16$Stat.Area)[c(
  which(unique(harvest16$Stat.Area) >= 25500 & unique(harvest16$Stat.Area) <=25799))]

Stat.Area.Mainland <- unique(harvest16$Stat.Area)[c(
  which(unique(harvest16$Stat.Area) >= 26200))]


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
cast(aggregate(Number ~ Geo + Strata, data = harvest16, sum), Geo ~ Strata, value = "Number")



HarvestByStrata2016 <- data.matrix(cast(aggregate(Number ~ Geo + Strata, data = harvest16, sum), Geo ~ Strata, value = "Number")[c(3, 1, 4, 2), 2:3])
dimnames(HarvestByStrata2016) <- list(KMA2016, c("1_Early", "2_Late"))
HarvestByStrata2016

HarvestByStrata2016["KALITC16", ] <- c(347, 427)
HarvestByStrata2016["KWESTC16", ] <- c(1004, 1094)
HarvestByStrata2016


dput(x = HarvestByStrata2016, file = "Objects/HarvestByStrata2016.txt")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot stock specific harvest results 2016 KMA Mixtures ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Multiply by harvest
HarvestByStrata2016


KMA2016Strata_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2016Strata_EstimatesStats.txt")
str(KMA2016Strata_EstimatesStats)
names(KMA2016Strata_EstimatesStats)
dimnames(KMA2016Strata_EstimatesStats[[1]])


KMA2016Strata_HarvestEstimatesStats <- sapply(names(KMA2016Strata_EstimatesStats), function(strata) {
  strata.split <- unlist(strsplit(x = strata, split = "_"))
  strata.split <- c(strata.split[1], paste(c(strata.split[2], strata.split[3]), collapse = "_"))
  
  cbind(KMA2016Strata_EstimatesStats[[strata]][, c("mean", "sd", "median", "5%", "95%")] * HarvestByStrata2016[strata.split[1], strata.split[2]],
        KMA2016Strata_EstimatesStats[[strata]][, c("P=0", "GR")])
}, simplify = FALSE )

dput(x = KMA2016Strata_HarvestEstimatesStats, file = "Estimates objects/Final/KMA2016Strata_HarvestEstimatesStats.txt")
KMA2016Strata_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2016Strata_HarvestEstimatesStats.txt")
str(KMA2016Strata_HarvestEstimatesStats)

# What should ymax be?
max(sapply(KMA2016Strata_HarvestEstimatesStats, function(strata) strata[, "95%"]))


TempMix16 <- sapply(KMA2016, function(geo) {grep(pattern = geo, x = names(KMA2016Strata_EstimatesStats), value = TRUE)} , simplify = FALSE)

Legend16 <- setNames(object = c("June 1-July 5", "July 6-August 5"), 
                     nm = c("1_Early", "2_Late"))
TempLegend16 <- sapply(KMA2016, function(geo) {
  Legend16[sapply(TempMix16[[geo]], function(strata) {unlist(strsplit(x = strata, split = paste(geo, "_", sep = '')))[2]} )]
} , simplify = FALSE)

GeoHeader <- setNames(object = c(paste0("SW Kodiak/Alitak 255, 256, 257"),
                                 paste0("Eastside Kodiak/Afognak 252, 258, 259"),
                                 paste0("Westside Kodiak/Afognak 251, 253, 254"),
                                 paste0("Mainland 262")),
                      nm = unlist(strsplit(x = KMA2016, split = "16")))


HarvestColors <- colorpanel(n = 2, low = "green", high = "white")
TempHarvestColors16 <- sapply(KMA2016, function(geo) {
  HarvestColors[sapply(TempMix16[[geo]], function(strata) {as.numeric(unlist(strsplit(x = strata, split = "_"))[2])} )]
} , simplify = FALSE)

Estimates <- KMA2016Strata_HarvestEstimatesStats
Groups <- groups10
Groups2Rows <- groups10tworows
cex.lab <- 1.5
cex.yaxis <- 1.3
cex.xaxis <- 0.6
cex.main <- 1.7
cex.leg <- 1.3
ci.lwd <- 2.5
ymax <- 1500

sapply(names(TempMix16), function(geomix) {
  emf(file = paste("Figures/2016/", geomix, "Harvest.emf", sep = ''), width = 8.5, height = 6.5, family = "serif", bg = "white")
  par(mar = c(2.1, 4.1, 2.6, 0.6))
  
  Barplot1 <- barplot2(height = t(sapply(TempMix16[[geomix]], function(tempmix) {Estimates[[tempmix]][, "median"]})), 
                       beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                       ci.l = t(sapply(TempMix16[[geomix]], function(tempmix) {Estimates[[tempmix]][, "5%"]})), 
                       ci.u = t(sapply(TempMix16[[geomix]], function(tempmix) {Estimates[[tempmix]][, "95%"]})), 
                       ylim = c(0, ymax), col = TempHarvestColors16[[geomix]], cex.axis = cex.yaxis, yaxt = "n", xaxt = 'n')
  axis(side = 2, at = seq(0, ymax, 500), labels = formatC(x = seq(0, ymax, 500), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
  legend(legend = TempLegend16[[geomix]], x = "topleft", fill = TempHarvestColors16[[geomix]], border = "black", bty = "n", cex = cex.leg, title="2016")
  abline(h = 0, xpd = FALSE)
  
  mtext(text = "Number of Fish Harvested", side = 2, cex = cex.yaxis, line = 3)
  mtext(text = Groups2Rows, side = 1, line = 1, at = apply(Barplot1, 2, mean), adj = 0.5, cex = cex.xaxis)
  mtext(text = GeoHeader[unlist(strsplit(geomix, split = "16"))], side = 3, cex = cex.main, line = 1)
  dev.off()
})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Stratified Regional Roll-Ups 2016 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
HarvestByStrata2016


# For this to work properly, the output mixture directories all need to be in the same place (can be moved back afterwards)
sapply(KMA2016, function(geomix) {
  assign(x = paste(geomix, "_Annual_Stratified", sep = ''), 
         value = StratifiedEstimator.GCL(
           groupvec = seq(groups10), groupnames = groups10, maindir = "BAYES/Output", 
           mixvec = paste(geomix, colnames(HarvestByStrata2016)[!is.na(HarvestByStrata2016[geomix,])], sep = "_"), 
           catchvec = HarvestByStrata2016[geomix, !is.na(HarvestByStrata2016[geomix,])], 
           newname = paste(geomix, "_Annual_Stratified", sep = ''), priorname = '',
           ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1), 
         pos = 1)
  dput(x = get(paste(geomix, "_Annual_Stratified", sep = '')), file = paste("Estimates objects/", geomix, "_Annual_Stratified.txt", sep = ''))
} ); beep(5)

str(KALITC16_Annual_Stratified)


# Create a list object with all Stratified Annual Rollups for the "Estimates objects/Final" folder
KMA2016_Annual_EstimatesStats <- sapply(KMA2016, function(geomix) {
  Stats <- get(paste(geomix, "_Annual_Stratified", sep = ''))$Stats
  Stats
}, simplify = FALSE)
str(KMA2016_Annual_EstimatesStats)
dput(x = KMA2016_Annual_EstimatesStats, file = "Estimates objects/Final/KMA2016_Annual_EstimatesStats.txt")




# Create a list object with all Stratified Annual Harvest
KMA2016_Annual_HarvestEstimatesStats <- sapply(KMA2016, function(strata) {
  cbind(KMA2016_Annual_EstimatesStats[[strata]][, c("mean", "sd", "median", "5%", "95%")] * sum(HarvestByStrata2016[strata, ], na.rm = TRUE),
        KMA2016_Annual_EstimatesStats[[strata]][, c("P=0", "GR")])
}, simplify = FALSE )

dput(x = KMA2016_Annual_HarvestEstimatesStats, file = "Estimates objects/Final/KMA2016_Annual_HarvestEstimatesStats.txt")
str(KMA2016_Annual_HarvestEstimatesStats)




# Create a matrix of annual means
Annual2016_Stratified_Estimates <- sapply(KMA2016, function(geomix) {
  KMA2016_Annual_EstimatesStats[[geomix]][, "mean"]
})

# Create a matrix of early strata means
KMA2016Strata_1_Early_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2016Strata_1_Early_EstimatesStats.txt")
EarlyStrata2016_Estimates <- sapply(KMA2016Strata_1_Early_EstimatesStats, function(geomix) {
  geomix[, "mean"]
})
colnames(EarlyStrata2016_Estimates) <- KMA2016[c(1, 2, 4, 3)]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot Annual vs. Early Means
par(mar = c(4.1, 5.1, 3.1, 1.1))
sapply(KMA2016, function(geomix) {
  Barplot <- barplot2(height = rbind(Annual2016_Stratified_Estimates[, geomix], EarlyStrata2016_Estimates[, geomix]) * 100, 
                      beside = TRUE, col = c("blue", "white"), yaxt = 'n', xaxt = 'n', main = geomix, ylab = "Precent of Mixture",
                      cex.lab = 2, cex.main = 2, ylim = c(0, 100))
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = ",", digits = 0, format = "f"), cex.axis = 1.5)
  abline(h = 0, xpd = FALSE)
  text(x = colMeans(Barplot), y = -1, labels = groups10tworows, srt = 90, adj = 1, xpd = TRUE, cex = 0.6)
  legend("topleft", legend = c("Annual Means", "Early Strata Means"), fill = c("blue", "white"), bty = 'n', cex = 1.5)
})



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Bubble Charts ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("V:/Analysis/4_Westward/Chinook/CSRI Westward Sport Harvest 2014-2016/Mixtures")

KMARS14_EstimatesStats <- dget(file = "Estimates objects/KMARS14_EstimatesStats.txt")
KMARS15_EstimatesStats <- dget(file = "Estimates objects/KMARS15_EstimatesStats.txt")
KMARS16_EstimatesStats <- dget(file = "Estimates objects/KMARS16_80K_EstimatesStats.txt")

setwd("V:/Analysis/4_Westward/Chinook/CSRI Westward Commercial Harvest 2014-2016/Mixtures")

groups10short <- c("Russia", "Eastern Bering", "North Peninsula", "Chignik", "Kodiak", "Cook Inlet", "Copper", "Southeast AK", "British Columbia", "West Coast US")
dput(x = groups10short, file = "Objects/groups10short.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create a matrix of annual means
#~~~~~~~~~~~~~~~~~~
# 2014
KMA2014_Annual_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2014_Annual_EstimatesStats.txt")
Annual2014_Stratified_Estimates <- sapply(KMA2014, function(geomix) {
  KMA2014_Annual_EstimatesStats[[geomix]][, "mean"]
})
rownames(Annual2014_Stratified_Estimates) <- groups10short


KMA2014_Annual_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2014_Annual_HarvestEstimatesStats.txt")
Annual2014_Stratified_HarvestEstimates <- sapply(KMA2014, function(geomix) {
  round(KMA2014_Annual_HarvestEstimatesStats[[geomix]][, "median"])
})
rownames(Annual2014_Stratified_HarvestEstimates) <- groups10short

Annual2014_Stratified_HarvestEstimates <- cbind(Annual2014_Stratified_HarvestEstimates, "KMARS14" = round(KMARS14_EstimatesStats$KMARS14[, "median"] * 8049))


#~~~~~~~~~~~~~~~~~~
# 2015
KMA2015_Annual_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2015_Annual_EstimatesStats.txt")
Annual2015_Stratified_Estimates <- sapply(KMA2015, function(geomix) {
  KMA2015_Annual_EstimatesStats[[geomix]][, "mean"]
})
rownames(Annual2015_Stratified_Estimates) <- groups10short


KMA2015_Annual_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2015_Annual_HarvestEstimatesStats.txt")
Annual2015_Stratified_HarvestEstimates <- sapply(KMA2015, function(geomix) {
  round(KMA2015_Annual_HarvestEstimatesStats[[geomix]][, "median"])
})
rownames(Annual2015_Stratified_HarvestEstimates) <- groups10short

Annual2015_Stratified_HarvestEstimates <- cbind(Annual2015_Stratified_HarvestEstimates, "KMARS15" = round(KMARS15_EstimatesStats$KMARS15[, "median"] * 6709))


#~~~~~~~~~~~~~~~~~~
# 2016
KMA2016_Annual_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2016_Annual_EstimatesStats.txt")
Annual2016_Stratified_Estimates <- sapply(KMA2016, function(geomix) {
  KMA2016_Annual_EstimatesStats[[geomix]][, "mean"]
})
rownames(Annual2016_Stratified_Estimates) <- groups10short


KMA2016_Annual_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2016_Annual_HarvestEstimatesStats.txt")
Annual2016_Stratified_HarvestEstimates <- sapply(KMA2016, function(geomix) {
  round(KMA2016_Annual_HarvestEstimatesStats[[geomix]][, "median"])
})
rownames(Annual2016_Stratified_HarvestEstimates) <- groups10short

Annual2016_Stratified_HarvestEstimates <- cbind(Annual2016_Stratified_HarvestEstimates, "KMARS16" = rep(0, 10))




require(lattice)
new.colors <- colorRampPalette(c("white", "black"))
levelplot(t(Annual2014_Stratified_Estimates[10:1, c(2, 3, 1, 4)]), col.regions = new.colors, xlab = "Reporting Group", 
          ylab = "Fishery", main = "2014 Proportion", at = seq(0, 1, length.out = 100), aspect = "fill", 
          scales = list(x = list(rot = 45)),
          panel = function(...) {
            panel.levelplot(...)
          })

require(lattice)
new.colors <- colorRampPalette(c("white", "black"))
levelplot(t(Annual2015_Stratified_Estimates[10:1, c(2, 3, 1, 4)]), col.regions = new.colors, xlab = "Reporting Group", 
          ylab = "Fishery", main = "2015 Proportion", at = seq(0, 1, length.out = 100), aspect = "fill", 
          scales = list(x = list(rot = 45)),
          panel = function(...) {
            panel.levelplot(...)
          })

require(lattice)
new.colors <- colorRampPalette(c("white", "black"))
levelplot(t(Annual2016_Stratified_Estimates[10:1, c(2, 3, 1, 4)]), col.regions = new.colors, xlab = "Reporting Group", 
          ylab = "Fishery", main = "2016 Proportion", at = seq(0, 1, length.out = 100), aspect = "fill", 
          scales = list(x = list(rot = 45)),
          panel = function(...) {
            panel.levelplot(...)
          })

require(ggplot2)
require(reshape2)

# 2014
Annual2014_Stratified_HarvestEstimates_df <- melt(Annual2014_Stratified_HarvestEstimates)
names(Annual2014_Stratified_HarvestEstimates_df) <- c("RG", "Fishery", "Harvest")
Annual2014_Stratified_HarvestEstimates_df$RG <- factor(Annual2014_Stratified_HarvestEstimates_df$RG, levels = rev(groups10short))
Annual2014_Stratified_HarvestEstimates_df$Fishery <- factor(Annual2014_Stratified_HarvestEstimates_df$Fishery, levels = KMA2014[c(2,3,1,4)])
Annual2014_Stratified_HarvestEstimates_df$Color <- rep(rev(colors10), 4)
str(Annual2014_Stratified_HarvestEstimates_df)

ggplot(data = Annual2014_Stratified_HarvestEstimates_df, aes(x = Fishery, y = RG, size = Harvest, color = RG)) + geom_point() + scale_size_area(max_size = 25) + scale_color_manual(values = rev(rep(colors10, 4)))

ggplot(data = Annual2014_Stratified_HarvestEstimates_df, aes(x = Fishery, y = RG, size = Harvest, color = RG)) + 
  geom_point() + 
  scale_size_continuous(limits = c(0, 1600), breaks = seq(500, 1500, 500), range = c(0, 25)) + 
  scale_color_manual(values = rev(rep(colors10, 4))) +
  ggtitle("2014 Harvest")


# Figure for Report
zmax = 4000; max(Annual2014_Stratified_HarvestEstimates, Annual2015_Stratified_HarvestEstimates, Annual2016_Stratified_HarvestEstimates)

Groups2RowsBubble <- Groups2Rows
Groups2RowsBubble[c(1,4,5,7,8)] <- gsub(pattern = "\n", replacement = "", x = Groups2Rows[c(1,4,5,7,8)])

colnames(Annual2014_Stratified_HarvestEstimates) <- c("SW/Alitak", "Eastside", "Westside", "Mainland", "Sport")

Annual2014_Stratified_HarvestEstimates_df <- melt(Annual2014_Stratified_HarvestEstimates)
names(Annual2014_Stratified_HarvestEstimates_df) <- c("RG", "Fishery", "Harvest")
Annual2014_Stratified_HarvestEstimates_df$RG <- factor(Annual2014_Stratified_HarvestEstimates_df$RG, levels = groups10short)
Annual2014_Stratified_HarvestEstimates_df$Fishery <- factor(Annual2014_Stratified_HarvestEstimates_df$Fishery, levels = rev(c("Westside", "SW/Alitak", "Eastside", "Mainland", "Sport")))
Annual2014_Stratified_HarvestEstimates_df$Color <- rep(colors10, 5)
Annual2014_Stratified_HarvestEstimates_df$Harvest[Annual2014_Stratified_HarvestEstimates_df$Harvest == 0] <- NA
str(Annual2014_Stratified_HarvestEstimates_df)

require(devEMF)
emf(file ="Figures/All Years/2014 Harvest Bubble Plot.emf", width = 9, height = 5.75, family = "serif", bg = "white")

ggplot(data = Annual2014_Stratified_HarvestEstimates_df, aes(x = RG, y = Fishery, size = Harvest, color = RG)) + 
  geom_point() + 
  scale_size_continuous(limits = c(0, zmax), breaks = c(50, 100, 200, 500, seq(1000, 4000, 1000)), range = c(0, 20), labels = c("50", "100", "200", "500", "1,000", "2,000", "3,000", "4,000")) +   scale_color_manual(values = rep(colors10, 5), guide = FALSE) +
  scale_x_discrete(name = "Reporting Group", labels = Groups2RowsBubble) +
  scale_y_discrete(name = "Sampling Area", labels = rev(c("NW Kodiak\nAfognak", "SW Kodiak\nAlitak", "Eastside\nKodiak", "Mainland", "Marine\nSport"))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = rel(1.8), angle = 90, margin = unit(c(0,0.2,0,0), "cm"))) +
  theme(axis.title.x = element_text(size = rel(1.8), angle = 00, margin = unit(c(0,0,0,0), "cm"))) +
  theme(legend.title = element_text(size = rel(1.8), angle = 00)) +
  theme(text = element_text(family = "times"))

dev.off()


# 2015
Annual2015_Stratified_HarvestEstimates_df <- melt(Annual2015_Stratified_HarvestEstimates)
names(Annual2015_Stratified_HarvestEstimates_df) <- c("RG", "Fishery", "Harvest")
Annual2015_Stratified_HarvestEstimates_df$RG <- factor(Annual2015_Stratified_HarvestEstimates_df$RG, levels = rev(groups10short))
Annual2015_Stratified_HarvestEstimates_df$Fishery <- factor(Annual2015_Stratified_HarvestEstimates_df$Fishery, levels = KMA2015[c(2,3,1,4)])
Annual2015_Stratified_HarvestEstimates_df$Color <- rep(rev(colors10), 4)
str(Annual2015_Stratified_HarvestEstimates_df)

ggplot(data = Annual2015_Stratified_HarvestEstimates_df, aes(x = Fishery, y = RG, size = Harvest, color = RG)) + geom_point() + scale_size_area(max_size = 25) + scale_color_manual(values = rev(rep(colors10, 4)))

ggplot(data = Annual2015_Stratified_HarvestEstimates_df, aes(x = Fishery, y = RG, size = Harvest, color = RG)) + 
  geom_point() + 
  scale_size_continuous(limits = c(0, 1600), breaks = seq(500, 1500, 500), range = c(0, 25)) + 
  scale_color_manual(values = rev(rep(colors10, 4))) +
  ggtitle("2015 Harvest")


# Figure for Report
zmax = 4000; max(Annual2015_Stratified_HarvestEstimates, Annual2015_Stratified_HarvestEstimates, Annual2016_Stratified_HarvestEstimates)

colnames(Annual2015_Stratified_HarvestEstimates) <- c("SW/Alitak", "Eastside", "Westside", "Mainland", "Sport")

Annual2015_Stratified_HarvestEstimates_df <- melt(Annual2015_Stratified_HarvestEstimates)
names(Annual2015_Stratified_HarvestEstimates_df) <- c("RG", "Fishery", "Harvest")
Annual2015_Stratified_HarvestEstimates_df$RG <- factor(Annual2015_Stratified_HarvestEstimates_df$RG, levels = groups10short)
Annual2015_Stratified_HarvestEstimates_df$Fishery <- factor(Annual2015_Stratified_HarvestEstimates_df$Fishery, levels = rev(c("Westside", "SW/Alitak", "Eastside", "Mainland", "Sport")))
Annual2015_Stratified_HarvestEstimates_df$Color <- rep(colors10, 5)
Annual2015_Stratified_HarvestEstimates_df$Harvest[Annual2015_Stratified_HarvestEstimates_df$Harvest == 0] <- NA
str(Annual2015_Stratified_HarvestEstimates_df)

require(devEMF)
emf(file ="Figures/All Years/2015 Harvest Bubble Plot.emf", width = 9, height = 5.75, family = "serif", bg = "white")

ggplot(data = Annual2015_Stratified_HarvestEstimates_df, aes(x = RG, y = Fishery, size = Harvest, color = RG)) + 
  geom_point() + 
  scale_size_continuous(limits = c(0, zmax), breaks = c(50, 100, 200, 500, seq(1000, 4000, 1000)), range = c(0, 20), labels = c("50", "100", "200", "500", "1,000", "2,000", "3,000", "4,000")) +   scale_color_manual(values = rep(colors10, 5), guide = FALSE) +
  scale_x_discrete(name = "Reporting Group", labels = Groups2RowsBubble) +
  scale_y_discrete(name = "Sampling Area", labels = rev(c("NW Kodiak\nAfognak", "SW Kodiak\nAlitak", "Eastside\nKodiak", "Mainland", "Marine\nSport"))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = rel(1.8), angle = 90, margin = unit(c(0,0.2,0,0), "cm"))) +
  theme(axis.title.x = element_text(size = rel(1.8), angle = 00, margin = unit(c(0,0,0,0), "cm"))) +
  theme(legend.title = element_text(size = rel(1.8), angle = 00)) +
  theme(text = element_text(family = "times"))

dev.off()




# 2016
Annual2016_Stratified_HarvestEstimates_df <- melt(Annual2016_Stratified_HarvestEstimates)
names(Annual2016_Stratified_HarvestEstimates_df) <- c("RG", "Fishery", "Harvest")
Annual2016_Stratified_HarvestEstimates_df$RG <- factor(Annual2016_Stratified_HarvestEstimates_df$RG, levels = rev(groups10short))
Annual2016_Stratified_HarvestEstimates_df$Fishery <- factor(Annual2016_Stratified_HarvestEstimates_df$Fishery, levels = KMA2016[c(2,3,1,4)])
Annual2016_Stratified_HarvestEstimates_df$Color <- rep(rev(colors10), 4)
str(Annual2016_Stratified_HarvestEstimates_df)

ggplot(data = Annual2016_Stratified_HarvestEstimates_df, aes(x = Fishery, y = RG, size = Harvest, color = RG)) + geom_point() + scale_size_area(max_size = 25) + scale_color_manual(values = rev(rep(colors10, 4)))

ggplot(data = Annual2016_Stratified_HarvestEstimates_df, aes(x = Fishery, y = RG, size = Harvest, color = RG)) + 
  geom_point() + 
  scale_size_continuous(limits = c(0, 1600), breaks = seq(500, 1500, 500), range = c(0, 25)) + 
  scale_color_manual(values = rev(rep(colors10, 4))) +
  ggtitle("2016 Harvest")


# Figure for Report
zmax = 4000; max(Annual2016_Stratified_HarvestEstimates, Annual2016_Stratified_HarvestEstimates, Annual2016_Stratified_HarvestEstimates)

colnames(Annual2016_Stratified_HarvestEstimates) <- c("SW/Alitak", "Eastside", "Westside", "Mainland", "Sport")

Annual2016_Stratified_HarvestEstimates_df <- melt(Annual2016_Stratified_HarvestEstimates)
names(Annual2016_Stratified_HarvestEstimates_df) <- c("RG", "Fishery", "Harvest")
Annual2016_Stratified_HarvestEstimates_df$RG <- factor(Annual2016_Stratified_HarvestEstimates_df$RG, levels = groups10short)
Annual2016_Stratified_HarvestEstimates_df$Fishery <- factor(Annual2016_Stratified_HarvestEstimates_df$Fishery, levels = rev(c("Westside", "SW/Alitak", "Eastside", "Mainland", "Sport")))
Annual2016_Stratified_HarvestEstimates_df$Color <- rep(colors10, 5)
Annual2016_Stratified_HarvestEstimates_df$Harvest[Annual2016_Stratified_HarvestEstimates_df$Harvest == 0] <- NA
str(Annual2016_Stratified_HarvestEstimates_df)

require(devEMF)
emf(file ="Figures/All Years/2016 Harvest Bubble Plot.emf", width = 9, height = 5.75, family = "serif", bg = "white")

ggplot(data = Annual2016_Stratified_HarvestEstimates_df, aes(x = RG, y = Fishery, size = Harvest, color = RG)) + 
  geom_point() + 
  scale_size_continuous(limits = c(0, zmax), breaks = c(50, 100, 200, 500, seq(1000, 4000, 1000)), range = c(0, 20), labels = c("50", "100", "200", "500", "1,000", "2,000", "3,000", "4,000")) +   scale_color_manual(values = rep(colors10, 5), guide = FALSE) +
  scale_x_discrete(name = "Reporting Group", labels = Groups2RowsBubble) +
  scale_y_discrete(name = "Sampling Area", labels = rev(c("NW Kodiak\nAfognak", "SW Kodiak\nAlitak", "Eastside\nKodiak", "Mainland", "Marine\nSport"))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = rel(1.8), angle = 90, margin = unit(c(0,0.2,0,0), "cm"))) +
  theme(axis.title.x = element_text(size = rel(1.8), angle = 00, margin = unit(c(0,0,0,0), "cm"))) +
  theme(legend.title = element_text(size = rel(1.8), angle = 00)) +
  theme(text = element_text(family = "times"))

dev.off()





Annual2016_Stratified_HarvestEstimates_df$Harvest <- NA
Annual2016_Stratified_HarvestEstimates_df$Harvest[1] <- 0

emf(file ="Figures/All Years/Blank Harvest Bubble Plot.emf", width = 9, height = 5.75, family = "serif", bg = "white")

ggplot(data = Annual2016_Stratified_HarvestEstimates_df, aes(x = RG, y = Fishery, size = Harvest, color = RG)) + 
  geom_point() + 
  scale_size_continuous(limits = c(0, zmax), breaks = c(50, 100, 200, 500, seq(1000, 4000, 1000)), range = c(0, 20), labels = c("50", "100", "200", "500", "1,000", "2,000", "3,000", "4,000")) + 
  scale_color_manual(values = rep(colors10, 5), guide = FALSE) +
  scale_x_discrete(name = "Reporting Group", labels = Groups2RowsBubble) +
  scale_y_discrete(name = "Sampling Area", labels = rev(c("NW Kodiak\nAfognak", "SW Kodiak\nAlitak", "Eastside\nKodiak", "Mainland", "Marine\nSport"))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = rel(1.8), angle = 90, margin = unit(c(0,0.2,0,0), "cm"))) +
  theme(axis.title.x = element_text(size = rel(1.8), angle = 00, margin = unit(c(0,0,0,0), "cm"))) +
  theme(legend.title = element_text(size = rel(1.8), angle = 00)) +
  theme(text = element_text(family = "times"))

dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Percentages for KMA Mixtures ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KMA2014Strata_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2014Strata_EstimatesStats.txt")
KMA2015Strata_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2015Strata_EstimatesStats.txt")
KMA2016Strata_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2016Strata_EstimatesStats.txt")

str(KMA2014Strata_EstimatesStats)




# Three barplot layout
layoutmat <- matrix(data=c(  1, 2,
                             1, 3,
                             1, 4,
                             5, 6), nrow = 4, ncol = 2, byrow = TRUE)



GeoMix <- c("KALITC", "KEASTC", "KWESTC", "KMAINC")

filenames <- setNames(object = c("SW Kodiak Alitak Proportions 2014-2016", 
                                 "Eastside Proportions 2014-2016",
                                 "Westside Proportions 2014-2016",
                                 "Mainland Proportions 2014-2016"), nm = GeoMix)

# If showing proportions (percetages) use blue, otherwise green as "low"
ProportionColors <- colorpanel(n = 2, low = "blue", high = "white")


#~~~~~~~~~~~~~~~~~~
# 2014
TempMix14 <- sapply(KMA2014, function(geo) {grep(pattern = geo, x = names(KMA2014Strata_EstimatesStats), value = TRUE)}, simplify = FALSE)

Legend14 <- setNames(object = c("June 1-July 5", "July 6-August 5"), 
                     nm = c("1_Early", "2_Late"))
TempLegend14 <- sapply(KMA2014, function(geo) {
  Legend14[sapply(TempMix14[[geo]], function(strata) {unlist(strsplit(x = strata, split = paste(geo, "_", sep = '')))[2]} )]
}, simplify = FALSE)


TempProportionColors14 <- sapply(KMA2014, function(geo) {
  ProportionColors[sapply(TempMix14[[geo]], function(strata) {as.numeric(unlist(strsplit(x = strata, split = "_"))[2])} )]
}, simplify = FALSE)

Estimates14 <- KMA2014Strata_EstimatesStats

#~~~~~~~~~~~~~~~~~~
# 2015
TempMix15 <- sapply(KMA2015, function(geo) {grep(pattern = geo, x = names(KMA2015Strata_EstimatesStats), value = TRUE)}, simplify = FALSE)

Legend15 <- setNames(object = c("June 1-July 5", "July 6-August 5"), 
                     nm = c("1_Early", "2_Late"))
TempLegend15 <- sapply(KMA2015, function(geo) {
  Legend15[sapply(TempMix15[[geo]], function(strata) {unlist(strsplit(x = strata, split = paste(geo, "_", sep = '')))[2]} )]
}, simplify = FALSE)


TempProportionColors15 <- sapply(KMA2015, function(geo) {
  ProportionColors[sapply(TempMix15[[geo]], function(strata) {as.numeric(unlist(strsplit(x = strata, split = "_"))[2])} )]
}, simplify = FALSE)

Estimates15 <- KMA2015Strata_EstimatesStats

#~~~~~~~~~~~~~~~~~~
# 2016
TempMix16 <- sapply(KMA2016, function(geo) {grep(pattern = geo, x = names(KMA2016Strata_EstimatesStats), value = TRUE)}, simplify = FALSE)

Legend16 <- setNames(object = c("June 1-July 5", "July 6-August 5"), 
                     nm = c("1_Early", "2_Late"))
TempLegend16 <- sapply(KMA2016, function(geo) {
  Legend16[sapply(TempMix16[[geo]], function(strata) {unlist(strsplit(x = strata, split = paste(geo, "_", sep = '')))[2]} )]
}, simplify = FALSE)


TempProportionColors16 <- sapply(KMA2016, function(geo) {
  ProportionColors[sapply(TempMix16[[geo]], function(strata) {as.numeric(unlist(strsplit(x = strata, split = "_"))[2])} )]
}, simplify = FALSE)

Estimates16 <- KMA2016Strata_EstimatesStats


#~~~~~~~~~~~~~~~~~~
# Size Parameters
Groups <- groups10
Groups2Rows <- groups10tworows
cex.lab <- 1.5
cex.xaxis <- 0.5
cex.yaxis <- 1.3
cex.leg <- 1.1
ci.lwd <- 2.5


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make figures as .emf files

# dir.create("Figures/All Years")
require(devEMF)
require(gplots)

sapply(GeoMix, function(geomix) {
  
  emf(file = paste("Figures/All Years/", filenames[geomix], ".emf", sep = ''), width = 6, height = 5.75, family = "serif", bg = "white")
  
  
  layout(mat = layoutmat, widths = c(0.075, 1, 1), heights = c(0.9, 0.9, 0.9, 0.15))
  par(mar = rep(0, 4))
  par(family = "times")
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Y-axis label
  plot.new()
  text(x = 0.25, y = 0.5, labels = "Percentage of Catch", srt = 90, cex = cex.lab)
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## 2014 Barplot
  geomix14 <- grep(pattern = geomix, x = names(TempMix14), value = TRUE)
  par(mar = c(1, 1, 1.5, 1))
  Barplot14 <- barplot2(height = t(sapply(TempMix14[[geomix14]], function(tempmix) {Estimates14[[tempmix]][, "median"]})) * 100, 
                        beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                        ci.l = t(sapply(TempMix14[[geomix14]], function(tempmix) {Estimates14[[tempmix]][, "5%"]})) * 100, 
                        ci.u = t(sapply(TempMix14[[geomix14]], function(tempmix) {Estimates14[[tempmix]][, "95%"]})) * 100, 
                        ylim = c(0, 100), col = TempProportionColors14[[geomix14]], yaxt = "n", xaxt = 'n')
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
  legend(legend = TempLegend14[[geomix14]], x = "topleft", fill = TempProportionColors14[[geomix14]], border = "black", bty = "n", cex = cex.leg, title="2014")
  abline(h = 0, xpd = FALSE)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## 2015 Barplot
  geomix15 <- grep(pattern = geomix, x = names(TempMix15), value = TRUE)
  par(mar = c(1, 1, 1.5, 1))
  Barplot15 <- barplot2(height = t(sapply(TempMix15[[geomix15]], function(tempmix) {Estimates15[[tempmix]][, "median"]})) * 100, 
                        beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                        ci.l = t(sapply(TempMix15[[geomix15]], function(tempmix) {Estimates15[[tempmix]][, "5%"]})) * 100, 
                        ci.u = t(sapply(TempMix15[[geomix15]], function(tempmix) {Estimates15[[tempmix]][, "95%"]})) * 100, 
                        ylim = c(0, 100), col = TempProportionColors15[[geomix15]], yaxt = "n", xaxt = 'n')
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
  legend(legend = TempLegend15[[geomix15]], x = "topleft", fill = TempProportionColors15[[geomix15]], border = "black", bty = "n", cex = cex.leg, title="2015")
  abline(h = 0, xpd = FALSE)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## 2016 Barplot
  geomix16 <- grep(pattern = geomix, x = names(TempMix16), value = TRUE)
  par(mar = c(1, 1, 1.5, 1))
  Barplot16 <- barplot2(height = t(sapply(TempMix16[[geomix16]], function(tempmix) {Estimates16[[tempmix]][, "median"]})) * 100,
                        beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                        ci.l = t(sapply(TempMix16[[geomix16]], function(tempmix) {Estimates16[[tempmix]][, "5%"]})) * 100,
                        ci.u = t(sapply(TempMix16[[geomix16]], function(tempmix) {Estimates16[[tempmix]][, "95%"]})) * 100,
                        ylim = c(0, 100), col = TempProportionColors16[[geomix16]], cex.axis = cex.yaxis, yaxt = "n", xaxt = 'n')
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
  legend(legend = TempLegend16[[geomix16]], x = "topleft", fill = TempProportionColors16[[geomix16]], border = "black", bty = "n", cex = cex.leg, title="2016")
  abline(h = 0, xpd = FALSE)

  mtext(text = Groups2Rows, side = 1, line = 1, at = apply(Barplot16, 2, mean), adj = 0.5, cex = cex.xaxis)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Blank Corner
  par(mar = rep(0, 4))
  plot.new()
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## x-axis label
  par(mar = rep(0, 4))
  plot.new()
  text(x = 0.5, y = 0.25, labels = "Reporting Group", cex = cex.lab)
  
  
  dev.off()
})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

emf(file = paste("Figures/All Years/Blank.emf", sep = ''), width = 6, height = 5.75, family = "serif", bg = "white")


layout(mat = layoutmat, widths = c(0.075, 1, 1), heights = c(0.9, 0.9, 0.9, 0.15))
par(mar = rep(0, 4))
par(family = "times")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Y-axis label
plot.new()
text(x = 0.25, y = 0.5, labels = "Percentage of Catch", srt = 90, cex = cex.lab)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2014 Barplot
par(mar = c(1, 1, 1.5, 1))
Barplot14 <- barplot2(height = matrix(data = 0, nrow = 2, ncol = 10), 
                      beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                      ci.l = matrix(data = 0, nrow = 2, ncol = 10), 
                      ci.u = matrix(data = 0, nrow = 2, ncol = 10), 
                      ylim = c(0, 100), col = ProportionColors, yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
legend(legend = c("Early", "Late"), x = "topleft", fill = ProportionColors, border = "black", bty = "n", cex = cex.leg, title="2014")
abline(h = 0, xpd = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2015 Barplot
par(mar = c(1, 1, 1.5, 1))
Barplot15 <- barplot2(height = matrix(data = 0, nrow = 2, ncol = 10), 
                      beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                      ci.l = matrix(data = 0, nrow = 2, ncol = 10), 
                      ci.u = matrix(data = 0, nrow = 2, ncol = 10), 
                      ylim = c(0, 100), col = ProportionColors, yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
legend(legend = c("Early", "Late"), x = "topleft", fill = ProportionColors, border = "black", bty = "n", cex = cex.leg, title="2015")
abline(h = 0, xpd = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2016 Barplot
par(mar = c(1, 1, 1.5, 1))
Barplot16 <- barplot2(height = matrix(data = 0, nrow = 2, ncol = 10), 
                      beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                      ci.l = matrix(data = 0, nrow = 2, ncol = 10), 
                      ci.u = matrix(data = 0, nrow = 2, ncol = 10), 
                      ylim = c(0, 100), col = ProportionColors, yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
legend(legend = c("Early", "Late"), x = "topleft", fill = ProportionColors, border = "black", bty = "n", cex = cex.leg, title="2016")
abline(h = 0, xpd = FALSE)

mtext(text = Groups2Rows, side = 1, line = 1, at = apply(Barplot16, 2, mean), adj = 0.5, cex = cex.xaxis)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Blank Corner
par(mar = rep(0, 4))
plot.new()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## x-axis label
par(mar = rep(0, 4))
plot.new()
text(x = 0.5, y = 0.25, labels = "Reporting Group", cex = cex.lab)


dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Harvests for KMA Mixtures ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KMA2014Strata_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2014Strata_HarvestEstimatesStats.txt")
KMA2015Strata_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2015Strata_HarvestEstimatesStats.txt")
KMA2016Strata_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2016Strata_HarvestEstimatesStats.txt")

str(KMA2014Strata_EstimatesStats)




# Three barplot layout
layoutmat <- matrix(data=c(  1, 2,
                             1, 3,
                             1, 4,
                             5, 6), nrow = 4, ncol = 2, byrow = TRUE)



GeoMix <- c("KALITC", "KEASTC", "KWESTC", "KMAINC")

filenames <- setNames(object = c("SW Kodiak Alitak Harvest 2014-2016", 
                                 "Eastside Harvest 2014-2016",
                                 "Westside Harvest 2014-2016",
                                 "Mainland Harvest 2014-2016"), nm = GeoMix)

# If showing proportions (percetages) use blue, otherwise green as "low"
HarvestColors <- colorpanel(n = 2, low = "green", high = "white")


#~~~~~~~~~~~~~~~~~~
# 2014
TempMix14 <- sapply(KMA2014, function(geo) {grep(pattern = geo, x = names(KMA2014Strata_EstimatesStats), value = TRUE)}, simplify = FALSE)

Legend14 <- setNames(object = c("June 1-July 5", "July 6-August 5"), 
                     nm = c("1_Early", "2_Late"))
TempLegend14 <- sapply(KMA2014, function(geo) {
  Legend14[sapply(TempMix14[[geo]], function(strata) {unlist(strsplit(x = strata, split = paste(geo, "_", sep = '')))[2]} )]
}, simplify = FALSE)


TempHarvestColors14 <- sapply(KMA2014, function(geo) {
  HarvestColors[sapply(TempMix14[[geo]], function(strata) {as.numeric(unlist(strsplit(x = strata, split = "_"))[2])} )]
}, simplify = FALSE)

HarvestEstimates14 <- KMA2014Strata_HarvestEstimatesStats

#~~~~~~~~~~~~~~~~~~
# 2015
TempMix15 <- sapply(KMA2015, function(geo) {grep(pattern = geo, x = names(KMA2015Strata_EstimatesStats), value = TRUE)}, simplify = FALSE)

Legend15 <- setNames(object = c("June 1-July 5", "July 6-August 5"), 
                     nm = c("1_Early", "2_Late"))
TempLegend15 <- sapply(KMA2015, function(geo) {
  Legend15[sapply(TempMix15[[geo]], function(strata) {unlist(strsplit(x = strata, split = paste(geo, "_", sep = '')))[2]} )]
}, simplify = FALSE)


TempHarvestColors15 <- sapply(KMA2015, function(geo) {
  HarvestColors[sapply(TempMix15[[geo]], function(strata) {as.numeric(unlist(strsplit(x = strata, split = "_"))[2])} )]
}, simplify = FALSE)

HarvestEstimates15 <- KMA2015Strata_HarvestEstimatesStats

#~~~~~~~~~~~~~~~~~~
# 2016
TempMix16 <- sapply(KMA2016, function(geo) {grep(pattern = geo, x = names(KMA2016Strata_EstimatesStats), value = TRUE)}, simplify = FALSE)

Legend16 <- setNames(object = c("June 1-July 5", "July 6-August 5"), 
                     nm = c("1_Early", "2_Late"))
TempLegend16 <- sapply(KMA2016, function(geo) {
  Legend16[sapply(TempMix16[[geo]], function(strata) {unlist(strsplit(x = strata, split = paste(geo, "_", sep = '')))[2]} )]
}, simplify = FALSE)


TempHarvestColors16 <- sapply(KMA2016, function(geo) {
  HarvestColors[sapply(TempMix16[[geo]], function(strata) {as.numeric(unlist(strsplit(x = strata, split = "_"))[2])} )]
}, simplify = FALSE)

HarvestEstimates16 <- KMA2016Strata_HarvestEstimatesStats

#~~~~~~~~~~~~~~~~~~
# Size Parameters
Groups <- groups10
Groups2Rows <- groups10tworows
cex.lab <- 1.5
cex.xaxis <- 0.5
cex.yaxis <- 1.3
cex.leg <- 1.1
ci.lwd <- 2.5
ymax <- 1500  # max(sapply(c(HarvestEstimates14, HarvestEstimates15, HarvestEstimates15), function(strata) {strata[, "95%"]}))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make figures as .emf files

# dir.create("Figures/All Years")
require(devEMF)
require(gplots)

sapply(GeoMix, function(geomix) {
  
  emf(file = paste("Figures/All Years/", filenames[geomix], ".emf", sep = ''), width = 6, height = 5.75, family = "serif", bg = "white")
  
  
  layout(mat = layoutmat, widths = c(0.075, 1, 1), heights = c(0.9, 0.9, 0.9, 0.15))
  par(mar = rep(0, 4))
  par(family = "times")
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Y-axis label
  plot.new()
  text(x = 0.25, y = 0.5, labels = "Number of Fish Harvested", srt = 90, cex = cex.lab)
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## 2014 Barplot
  geomix14 <- grep(pattern = geomix, x = names(TempMix14), value = TRUE)
  par(mar = c(1, 1, 1.5, 1))
  Barplot14 <- barplot2(height = t(sapply(TempMix14[[geomix14]], function(tempmix) {HarvestEstimates14[[tempmix]][, "median"]})), 
                        beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                        ci.l = t(sapply(TempMix14[[geomix14]], function(tempmix) {HarvestEstimates14[[tempmix]][, "5%"]})), 
                        ci.u = t(sapply(TempMix14[[geomix14]], function(tempmix) {HarvestEstimates14[[tempmix]][, "95%"]})), 
                        ylim = c(0, ymax), col = TempHarvestColors14[[geomix14]], yaxt = "n", xaxt = 'n')
  axis(side = 2, at = seq(0, ymax, 500), labels = formatC(x = seq(0, ymax, 500), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
  legend(legend = TempLegend14[[geomix14]], x = "topleft", fill = TempHarvestColors14[[geomix14]], border = "black", bty = "n", cex = cex.leg, title="2014")
  abline(h = 0, xpd = FALSE)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## 2015 Barplot
  geomix15 <- grep(pattern = geomix, x = names(TempMix15), value = TRUE)
  par(mar = c(1, 1, 1.5, 1))
  Barplot15 <- barplot2(height = t(sapply(TempMix15[[geomix15]], function(tempmix) {HarvestEstimates15[[tempmix]][, "median"]})), 
                        beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                        ci.l = t(sapply(TempMix15[[geomix15]], function(tempmix) {HarvestEstimates15[[tempmix]][, "5%"]})), 
                        ci.u = t(sapply(TempMix15[[geomix15]], function(tempmix) {HarvestEstimates15[[tempmix]][, "95%"]})), 
                        ylim = c(0, ymax), col = TempHarvestColors15[[geomix15]], yaxt = "n", xaxt = 'n')
  axis(side = 2, at = seq(0, ymax, 500), labels = formatC(x = seq(0, ymax, 500), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
  legend(legend = TempLegend15[[geomix15]], x = "topleft", fill = TempHarvestColors15[[geomix15]], border = "black", bty = "n", cex = cex.leg, title="2015")
  abline(h = 0, xpd = FALSE)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## 2016 Barplot
  geomix16 <- grep(pattern = geomix, x = names(TempMix16), value = TRUE)
  par(mar = c(1, 1, 1.5, 1))
  Barplot16 <- barplot2(height = t(sapply(TempMix16[[geomix16]], function(tempmix) {HarvestEstimates16[[tempmix]][, "median"]})), 
                        beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                        ci.l = t(sapply(TempMix16[[geomix16]], function(tempmix) {HarvestEstimates16[[tempmix]][, "5%"]})), 
                        ci.u = t(sapply(TempMix16[[geomix16]], function(tempmix) {HarvestEstimates16[[tempmix]][, "95%"]})), 
                        ylim = c(0, ymax), col = TempHarvestColors16[[geomix16]], yaxt = "n", xaxt = 'n')
  axis(side = 2, at = seq(0, ymax, 500), labels = formatC(x = seq(0, ymax, 500), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
  legend(legend = TempLegend16[[geomix16]], x = "topleft", fill = TempHarvestColors16[[geomix16]], border = "black", bty = "n", cex = cex.leg, title="2016")
  abline(h = 0, xpd = FALSE)
  
  mtext(text = Groups2Rows, side = 1, line = 1, at = apply(Barplot16, 2, mean), adj = 0.5, cex = cex.xaxis)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Blank Corner
  par(mar = rep(0, 4))
  plot.new()
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## x-axis label
  par(mar = rep(0, 4))
  plot.new()
  text(x = 0.5, y = 0.25, labels = "Reporting Group", cex = cex.lab)

  dev.off()
})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot South Pen Chignik ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
KMA2014Strata_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2014Strata_EstimatesStats.txt")

SPENCHIG14_EstimatesStats <- KMA2014Strata_EstimatesStats[["KSPENCHIG14"]]
SPENCHIG14_HarvestEstimatesStats <- cbind(KMA2014Strata_EstimatesStats[["KSPENCHIG14"]][, c("mean", "sd", "median", "5%", "95%")] * 12209,
                                          KMA2014Strata_EstimatesStats[["KSPENCHIG14"]][, c("P=0", "GR")])

TempLegend <- "June 1-August 5"

# Three barplot layout
layoutmat <- matrix(data=c(  1, 2,
                             3, 4,
                             5, 6), nrow = 3, ncol = 2, byrow = TRUE)



GeoMix <- c("KSPENCHIG14")

filenames <- setNames(object = c("South Pen Chignik Proportions Harvest 2014"), nm = GeoMix)

# If showing proportions (percetages) use blue, otherwise green as "low"
ProportionColors <- "blue"
HarvestColors <- "green"

#~~~~~~~~~~~~~~~~~~
# Size Parameters
Groups <- groups10
Groups2Rows <- groups10tworows
cex.lab <- 1.5
cex.xaxis <- 0.5
cex.yaxis <- 1.3
cex.leg <- 1.1
ci.lwd <- 2.5
ymax <- 6000  # max(SPENCHIG14_HarvestEstimatesStats[, "95%"])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make figures as .emf files

# dir.create("Figures/All Years")
require(devEMF)
require(gplots)


emf(file = paste("Figures/Final/", filenames, ".emf", sep = ''), width = 6, height = 5.75, family = "serif", bg = "white")


layout(mat = layoutmat, widths = c(0.075, 1), heights = c(1.35, 1.35, 0.15))
par(mar = rep(0, 4))
par(family = "times")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Y-axis label
plot.new()
par(mar = rep(0, 4))
text(x = 0.25, y = 0.5, labels = "Percentage of Catch", srt = 90, cex = cex.lab)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Proportions Barplot
par(mar = c(1, 1, 1.5, 1))
BarplotProp <- barplot2(height = SPENCHIG14_EstimatesStats[,"median"] * 100, 
                      beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                      ci.l = SPENCHIG14_EstimatesStats[,"5%"] * 100, 
                      ci.u = SPENCHIG14_EstimatesStats[,"95%"] * 100, 
                      ylim = c(0, 100), col = ProportionColors, yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
legend(legend = TempLegend, x = "topleft", fill = ProportionColors, border = "black", bty = "n", cex = cex.leg, title="2014")
abline(h = 0, xpd = FALSE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Y-axis label
par(mar = rep(0, 4))
plot.new()
text(x = 0.25, y = 0.5, labels = "Number of Fish Harvested", srt = 90, cex = cex.lab)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Harvest Barplot
par(mar = c(1, 1, 1.5, 1))
BarplotHarvest <- barplot2(height = SPENCHIG14_HarvestEstimatesStats[,"median"], 
                        beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                        ci.l = SPENCHIG14_HarvestEstimatesStats[,"5%"], 
                        ci.u = SPENCHIG14_HarvestEstimatesStats[,"95%"], 
                        ylim = c(0, ymax), col = HarvestColors, yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, ymax, 2000), labels = formatC(x = seq(0, ymax, 2000), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
legend(legend = TempLegend, x = "topleft", fill = HarvestColors, border = "black", bty = "n", cex = cex.leg, title="2014")
abline(h = 0, xpd = FALSE)


mtext(text = Groups2Rows, side = 1, line = 1, at = BarplotHarvest, adj = 0.5, cex = cex.xaxis)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Blank Corner
par(mar = rep(0, 4))
plot.new()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## x-axis label
par(mar = rep(0, 4))
plot.new()
text(x = 0.5, y = 0.25, labels = "Reporting Group", cex = cex.lab)

dev.off()



























KMA2014Strata_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2014Strata_EstimatesStats.txt")

SPENCHIG14_EstimatesStats <- KMA2014Strata_EstimatesStats[["KSPENCHIG14"]]
SPENCHIG14_HarvestEstimatesStats <- cbind(KMA2014Strata_EstimatesStats[["KSPENCHIG14"]][, c("mean", "sd", "median", "5%", "95%")] * 12209,
                                          KMA2014Strata_EstimatesStats[["KSPENCHIG14"]][, c("P=0", "GR")])

TempLegend <- "June 1-August 5"

# Three barplot layout
layoutmat <- matrix(data=c(  1, 2, 3,
                             4, 5, 6), nrow = 2, ncol = 3, byrow = TRUE)



GeoMix <- c("KSPENCHIG14")

filenames <- setNames(object = c("South Pen Chignik Proportions Harvest 2014 OnePlot"), nm = GeoMix)

# If showing proportions (percetages) use blue, otherwise green as "low"

#~~~~~~~~~~~~~~~~~~
# Size Parameters
Groups <- groups10
Groups2Rows <- groups10tworows
cex.lab <- 1.5
cex.xaxis <- 0.5
cex.yaxis <- 1.3
cex.leg <- 1.1
ci.lwd <- 2.5

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make figures as .emf files

# dir.create("Figures/All Years")
require(devEMF)
require(gplots)


emf(file = paste("Figures/Final/", filenames, ".emf", sep = ''), width = 6, height = 5.75, family = "serif", bg = "white")


layout(mat = layoutmat, widths = c(0.075, 1, 0.075), heights = c(2.7, 0.15))
par(mar = rep(0, 4))
par(family = "times")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Left Y-axis label
par(mar = rep(0, 4))
plot.new()
text(x = 0.25, y = 0.5, labels = "Percentage of Catch", srt = 90, cex = cex.lab)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Proportions Barplot
par(mar = c(1, 1, 1.5, 1))
BarplotProp <- barplot2(height = SPENCHIG14_EstimatesStats[,"median"] * 100, 
                        beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                        ci.l = SPENCHIG14_EstimatesStats[,"5%"] * 100, 
                        ci.u = SPENCHIG14_EstimatesStats[,"95%"] * 100, 
                        ylim = c(0, 100), col = "white", yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
legend(legend = TempLegend, x = "topleft", fill = "white", border = "black", bty = "n", cex = cex.leg, title="2014")
abline(h = 0, xpd = FALSE)
axis(side = 4, at = c(0, 3000, 6000, 9000, 12000) / 12209 * 100, labels = formatC(x = c(0, 3000, 6000, 9000, 12000), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)

mtext(text = Groups2Rows, side = 1, line = 1, at = BarplotHarvest, adj = 0.5, cex = cex.xaxis)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Y-axis label
par(mar = rep(0, 4))
plot.new()
text(x = 0.4, y = 0.5, labels = "Number of Fish Harvested", srt = 90, cex = cex.lab, adj = c(0.5, 1))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Blank Corner
par(mar = rep(0, 4))
plot.new()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## x-axis label
par(mar = rep(0, 4))
plot.new()
text(x = 0.5, y = 0.25, labels = "Reporting Group", cex = cex.lab)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Blank Corner
par(mar = rep(0, 4))
plot.new()

dev.off()


























#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Table results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Dates
DatesStrata2014_Final <- matrix(data = c(rep("June 1-July 5", 3), "", rep("July 6-August 5", 4)), 
                                nrow = 4, ncol = 2, byrow = FALSE, dimnames = list(KMA2014, c("1_Early", "2_Late")))
dput(x = DatesStrata2014_Final, file = "Objects/DatesStrata2014_Final.txt")
DatesStrata2014_Final <- dget(file = "Objects/DatesStrata2014_Final.txt")


DatesStrata2015_Final <- matrix(data = c(rep("June 1-July 5", 3), "", rep("July 6-August 5", 4)), 
                                nrow = 4, ncol = 2, byrow = FALSE, dimnames = list(KMA2015, c("1_Early", "2_Late")))
dput(x = DatesStrata2015_Final, file = "Objects/DatesStrata2015_Final.txt")
DatesStrata2015_Final <- dget(file = "Objects/DatesStrata2015_Final.txt")


DatesStrata2016_Final <- matrix(data = c(rep("June 1-July 5", 4), rep("July 6-August 5", 4)), 
                                nrow = 4, ncol = 2, byrow = FALSE, dimnames = list(KMA2016, c("1_Early", "2_Late")))
dput(x = DatesStrata2016_Final, file = "Objects/DatesStrata2016_Final.txt")
DatesStrata2016_Final <- dget(file = "Objects/DatesStrata2016_Final.txt")


## Sample sizes
# KMA2014_2015Strata_SampleSizes_Final <- KMA2014_2015Strata_SampleSizes[, "Final"]
# dput(x = KMA2014_2015Strata_SampleSizes_Final, file = "Objects/KMA2014_2015Strata_SampleSizes_Final.txt")
# KMA2014_2015Strata_SampleSizes_Final <- dget(file = "Objects/KMA2014_2015Strata_SampleSizes_Final.txt")


KMA2014_2016Strata_SampleSizes_Final <- KMA2014_2016Strata_SampleSizes[, "Final"]
dput(x = KMA2014_2016Strata_SampleSizes_Final, file = "Objects/KMA2014_2016Strata_SampleSizes_Final.txt")
KMA2014_2016Strata_SampleSizes_Final <- dget(file = "Objects/KMA2014_2016Strata_SampleSizes_Final.txt")



## Geographic headers
GeoHeader <- setNames(object = c(paste0("Southwest Kodiak/Alitak (District 255, 256, 257)"),
                                 paste0("Eastside Kodiak/Afognak (District 252, 258, 259)"),
                                 paste0("Northwest Kodiak/Afognak (District 251, 253, 254)"),
                                 paste0("Mainland (District 262)")),
                      nm = unlist(strsplit(x = KMA2015, split = "15")))
dput(x = GeoHeader, file = "Objects/GeoHeader.txt")
GeoHeader <- dget(file = "Objects/GeoHeader.txt")


# Get Final Estimates Objects
KMAfinalestimatesobjects <- list.files(path = "Estimates objects/Final", recursive = FALSE)
invisible(sapply(KMAfinalestimatesobjects, function(objct) {assign(x = unlist(strsplit(x = objct, split = ".txt")), value = dget(file = paste("Estimates objects/Final", objct, sep = "/")), pos = 1) })); beep(2)
KMAfinalestimatesobjects; rm(KMAfinalestimatesobjects)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Defining caption variables

EstimatesStats <- c(KMA2014Strata_EstimatesStats, KMA2014_Annual_EstimatesStats,
                    KMA2015Strata_EstimatesStats, KMA2015_Annual_EstimatesStats,
                    KMA2016Strata_EstimatesStats, KMA2016_Annual_EstimatesStats)

KMA2014Strata_HarvestEstimatesStats <- c(KMA2014Strata_HarvestEstimatesStats,
                                         list("KSPENCHIG14" = 
                                                cbind(KMA2014Strata_EstimatesStats[["KSPENCHIG14"]][, c("mean", "sd", "median", "5%", "95%")] * 12209,
                                                      KMA2014Strata_EstimatesStats[["KSPENCHIG14"]][, c("P=0", "GR")])))


HarvestEstimatesStats <- c(KMA2014Strata_HarvestEstimatesStats, KMA2014_Annual_HarvestEstimatesStats,
                           KMA2015Strata_HarvestEstimatesStats, KMA2015_Annual_HarvestEstimatesStats,
                           KMA2016Strata_HarvestEstimatesStats, KMA2016_Annual_HarvestEstimatesStats)



SheetNames <- sort(names(EstimatesStats))
names(SheetNames) <- SheetNames
mixvec <- SheetNames


# mix <- SheetNames[2]
harvest <- rbind(HarvestByStrata2014, HarvestByStrata2015, HarvestByStrata2016)
dates <- rbind(DatesStrata2014_Final, DatesStrata2015_Final, DatesStrata2016_Final)
sampsize <- KMA2014_2016Strata_SampleSizes_Final

SheetNames <- SheetNames[-which(SheetNames == "KSPENCHIG14")]

SheetNames <- grep(pattern = "16", x = SheetNames, value = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dir.create("Estimates tables")
require(xlsx)

for(mix in SheetNames) {
  
  # String split by "_
  SheetNames.split <- unlist(strsplit(x = mix, split = "_"))
  
  
  # The first element is the geographic area + year
  geomix <- SheetNames.split[1]
  yr <- unlist(strsplit(x = geomix, split = "C"))[2]
  geo <- unlist(strsplit(x = geomix, split = "1"))[1]
  
  
  # If it is not an annual roll-up, then get the strata number + strata name
  if(length(SheetNames.split) > 1) {
    tempmix <- paste(c(SheetNames.split[2], SheetNames.split[3]), collapse = "_")
    
    Caption <- paste("Table X.-Estimates of stock composition (%) and stock-specific harvest for temporal stratum ",
                     SheetNames.split[2], " (", dates[geomix, tempmix],
                     "; Harvest=", formatC(x = harvest[geomix, tempmix], format = "f", digits = 0, big.mark = ","),
                     "; n=", sampsize[mix], ")", " of the ", GeoHeader[geo], ", 20", yr,
                     ". Estimates include median, 90% credibility interval (CI), the probability that the group estimate is equal to zero (P=0), mean, and SD.",
                     sep = '')
  } else {
    Caption <- paste("Table X.-Annual estimates of stock composition (%) and stock-specific harvest for the ", GeoHeader[geo], ", 20", yr,
                     ". Estimates include median, 90% credibility interval (CI), the probability that the group estimate is equal to zero (P=0), mean, and SD.",
                     sep = '')
  }
  
  Disclaimer <- "Note: Stock composition estimates may not sum to 100% and stock-specific harvest estimates may not sum to the total harvest due to rounding error."
  
  
  TableX <- matrix(data = "", nrow = 16, ncol = 13)
  
  TableX[1, 1] <- Caption
  TableX[2, c(2, 9)] <- c("Stock Composition", "Stock-specific Harvest")
  TableX[3, c(3, 10)] <- rep("90% CI", 2)
  TableX[4, c(1, 2:4, 6:7, 9:13, 5)] <- c("Reporting Group", rep(c("Median", "5%", "95%", "Mean", "SD"), 2), "P=0")
  TableX[5:14, 1] <- groups10
  TableX[5:14, c(2:4, 6:7)] <- formatC(x = EstimatesStats[[mix]][, c("median", "5%", "95%", "mean", "sd")] * 100, digits = 1, format = "f")
  TableX[5:14, 5] <- formatC(x = EstimatesStats[[mix]][, "P=0"], digits = 2, format = "f")
  TableX[5:14, 9:13] <- formatC(x = HarvestEstimatesStats[[mix]][, c("median", "5%", "95%", "mean", "sd")], digits = 0, format = "f", big.mark = ",")
  TableX[15, 11:12] <- c("Total", formatC(x = sum(HarvestEstimatesStats[[mix]][, "mean"]), digits = 0, format = "f", big.mark = ","))
  TableX[16, 1] <- Disclaimer
  
  
  write.xlsx(x = as.data.frame(TableX), 
             file = "Estimates tables/KMA Chinook Estimates Tables.xlsx",
             col.names = FALSE, row.names = FALSE, append = TRUE, sheetName = mix)
}; beep(5)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Chignik South Pen on its own, no expansion to harvest
mix <- setNames(object = "KSPENCHIG14", nm = "KSPENCHIG14")

# If it is not an annual roll-up, then get the strata number + strata name
Caption <- paste0("Table X.-Annual estimates of stock composition (%) for the South Peninsula / Chignik Outside (District 282, 283, 284, 285, 272, 273, 275), 2014. Estimates include median, 90% credibility interval (CI), the probability that the group estimate is equal to zero (P=0), mean, and SD.")

Disclaimer <- "Note: Stock composition estimates may not sum to 100% due to rounding error."


TableX <- matrix(data = "", nrow = 16, ncol = 7)

TableX[1, 1] <- Caption
TableX[2, 2] <- "Stock Composition"
TableX[3, 3] <- "90% CI"
TableX[4, c(1, 2:4, 6:7, 5)] <- c("Reporting Group", c("Median", "5%", "95%", "Mean", "SD"), "P=0")
TableX[5:14, 1] <- groups10
TableX[5:14, c(2:4, 6:7)] <- formatC(x = EstimatesStats[[mix]][, c("median", "5%", "95%", "mean", "sd")] * 100, digits = 1, format = "f")
TableX[5:14, 5] <- formatC(x = EstimatesStats[[mix]][, "P=0"], digits = 2, format = "f")
TableX[16, 1] <- Disclaimer


write.xlsx(x = as.data.frame(TableX), 
           file = "Estimates tables/KMA Chinook Estimates Tables.xlsx",
           col.names = FALSE, row.names = FALSE, append = TRUE, sheetName = mix)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### South Pen Chignik 2014 Individual Assignment ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
KSPENCHIG14_IA <- IndividualAssignmentSummary.GCL(GroupNames = groups10short, groupvec = groupvec10, mixnames = "KSPENCHIG14", BAYESoutputDir = "BAYES/Output", nchains = 5, nreps = 40000, burn = 1/2, thin = 100)
KSPENCHIG14_IA <- IndividualAssignmentSummary.GCL(GroupNames = Groups15, groupvec = groupvec15, mixnames = "KSPENCHIG14", BAYESoutputDir = "BAYES/Output", nchains = 5, nreps = 40000, burn = 1/2, thin = 100)
str(KSPENCHIG14_IA)

# Plot Indivdiual Assignment probabilites
require(lattice)
new.colors <- colorRampPalette(c("white", "black"))
levelplot(t(KSPENCHIG14_IA$KSPENCHIG14), 
          col.regions = new.colors, 
          at = seq(from = 0, to = 1, length.out = 100), 
          main = "Individual Assignment", xlab = "Group", ylab = "Ind", 
          scales = list(x = list(rot = 90)), 
          aspect = "fill")  # aspect = "iso" will make squares

apply(KSPENCHIG14_IA$KSPENCHIG14, 2, hist)


# Assign individuals to most likely
KSPENCHIG14.gcl$attributes$maxIA <- groups10short[apply(KSPENCHIG14_IA$KSPENCHIG14, 1, which.max)]
KSPENCHIG14.gcl$attributes$maxIA <- Groups15[apply(KSPENCHIG14_IA$KSPENCHIG14, 1, which.max)]

# Create Temporal Strata
KSPENCHIG14.gcl$attributes$MESH_SIZE_COMMENT <- ifelse(as.Date(KSPENCHIG14.gcl$attributes$CAPTURE_DATE) < as.Date("2014-06-30"), "Early", "Late")

# View IA by Spatio-Temporal Strata
table(KSPENCHIG14.gcl$attributes$maxIA, KSPENCHIG14.gcl$attributes$CAPTURE_LOCATION, KSPENCHIG14.gcl$attributes$MESH_SIZE_COMMENT)

# Write fish ID
writeClipboard(sapply(KSPENCHIG14.gcl$attributes$SillySource, function(x) {unlist(strsplit(x = x, split = "_"))[2]}))

# Pair size data based on fish ID (vial)
KSPENCHIG14.gcl$attributes$Size <- as.numeric(readClipboard())
KSPENCHIG14.gcl$attributes$Size[KSPENCHIG14.gcl$attributes$Size == 0] <- NA


# Pair sex data based on fish ID (vial)
KSPENCHIG14.gcl$attributes$Sex <- as.numeric(readClipboard())
KSPENCHIG14.gcl$attributes$Sex[KSPENCHIG14.gcl$attributes$Sex == 0] <- NA

table(KSPENCHIG14.gcl$attributes$Sex, KSPENCHIG14.gcl$attributes$maxIA)



hist(KSPENCHIG14.gcl$attributes$Size, breaks = seq(300, 1000, 20), col = 8)


hist(KSPENCHIG14.gcl$attributes$Size, breaks = seq(300, 1000, 20), col = 8, main = "Histogram of Size for\nSouthPen Chignik 2014 Mixture", xlab = "Size (mm)")
hist(KSPENCHIG14.gcl$attributes$Size[KSPENCHIG14.gcl$attributes$maxIA == "CWAK / Yukon"], breaks = seq(300, 1000, 20), col = 4, add = TRUE)
legend("topleft", legend = c("All MSA Samples", "CWAK / Yukon"), fill = c(8, 4), bty = "n")

hist(KSPENCHIG14.gcl$attributes$Size, breaks = seq(300, 1000, 20), col = 8, main = "Histogram of Size for\nSouthPen Chignik 2014 Mixture", xlab = "Size (mm)")
hist(KSPENCHIG14.gcl$attributes$Size[KSPENCHIG14.gcl$attributes$maxIA == "CWAK / Yukon"], breaks = seq(300, 1000, 20), col = 2, add = TRUE)
hist(KSPENCHIG14.gcl$attributes$Size[KSPENCHIG14.gcl$attributes$maxIA == "CWAK / Yukon" & KSPENCHIG14.gcl$attributes$MESH_SIZE_COMMENT == "Late"], breaks = seq(300, 1000, 20), col = 4, add = TRUE)
legend("topleft", legend = c("All MSA Samples", "CWAK / Yukon Early", "CWAK / Yukon Late"), fill = c(8, 2, 4), bty = "n")

hist(KSPENCHIG14.gcl$attributes$Size, breaks = seq(300, 1000, 20), col = 8, main = "Histogram of Size for\nSouthPen Chignik 2014 Mixture", xlab = "Size (mm)")
hist(KSPENCHIG14.gcl$attributes$Size[KSPENCHIG14.gcl$attributes$maxIA == "CWAK / Yukon" & KSPENCHIG14.gcl$attributes$MESH_SIZE_COMMENT == "Early"], breaks = seq(300, 1000, 20), col = rgb(1,0,0,0.5), add = TRUE)
hist(KSPENCHIG14.gcl$attributes$Size[KSPENCHIG14.gcl$attributes$maxIA == "CWAK / Yukon" & KSPENCHIG14.gcl$attributes$MESH_SIZE_COMMENT == "Late"], breaks = seq(300, 1000, 20), col = rgb(0,0,1,0.5), add = TRUE)
legend("topleft", legend = c("All MSA Samples", "CWAK / Yukon Early", "CWAK / Yukon Late"), fill = c(8, 2, 4), bty = "n")



hist(KSPENCHIG14.gcl$attributes$Size, breaks = seq(300, 1000, 20), col = 4, main = "Histogram of Size for\nSouthPen Chignik 2014 Mixture", xlab = "Size (mm)")
hist(KSPENCHIG14.gcl$attributes$Size[KSPENCHIG14.gcl$attributes$MESH_SIZE_COMMENT == "Early"], breaks = seq(300, 1000, 20), col = 2, add = TRUE)
legend("topleft", legend = c("Early", "Late"), fill = c(2, 4), bty = "n")






hist(KSPENCHIG14.gcl$attributes$Size[KSPENCHIG14.gcl$attributes$CAPTURE_LOCATION %in% c("SE/South Central", "Unimak/Southwestern")], breaks = seq(300, 1000, 20), col = 8, main = "Histogram of Size for\nonly SouthPen fish in 2014 Mixture", xlab = "Size (mm)")
hist(KSPENCHIG14.gcl$attributes$Size[KSPENCHIG14.gcl$attributes$maxIA == "CWAK / Yukon" & KSPENCHIG14.gcl$attributes$CAPTURE_LOCATION %in% c("SE/South Central", "Unimak/Southwestern")], breaks = seq(300, 1000, 20), col = 2, add = TRUE)
hist(KSPENCHIG14.gcl$attributes$Size[KSPENCHIG14.gcl$attributes$maxIA == "CWAK / Yukon" & KSPENCHIG14.gcl$attributes$CAPTURE_LOCATION %in% c("SE/South Central", "Unimak/Southwestern") & KSPENCHIG14.gcl$attributes$MESH_SIZE_COMMENT == "Late"], breaks = seq(300, 1000, 20), col = 4, add = TRUE)
legend("topleft", legend = c("All MSA Samples", "CWAK / Yukon Early", "CWAK / Yukon Late"), fill = c(8, 2, 4), bty = "n")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Stratified Annual Rollups Across Spatio-Temporal Strata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
HarvestByStrata2014
HarvestByStrata2015
HarvestByStrata2016

KMA2014
KMA2014Strata

#~~~~~~~~~~~~~~~~~~
# For this to work properly, the output mixture directories all need to be in the same place (can be moved back afterwards)
# 2014
KMA2014_Annual_Stratified <- StratifiedEstimator.GCL(
  groupvec = seq(groups10), groupnames = groups10, maindir = "BAYES/Output", 
  mixvec = KMA2014Strata, 
  catchvec = as.vector(na.omit(as.vector(t(HarvestByStrata2014[sort(KMA2014), ])))), 
  newname = "KMA2014_Annual_Stratified", priorname = '',
  ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1
)
str(KMA2014_Annual_Stratified)
dput(x = KMA2014_Annual_Stratified, file = "Estimates objects/KMA2014_Annual_Stratified.txt")
KMA2014_Annual_Stratified_EstimatesStats <- KMA2014_Annual_Stratified$Stats
dput(x = KMA2014_Annual_Stratified_EstimatesStats, file = "Estimates objects/Final/KMA2014_Annual_Stratified_EstimatesStats.txt")

#~~~~~~~~~~~~~~~~~~
#2015
KMA2015_Annual_Stratified <- StratifiedEstimator.GCL(
  groupvec = seq(groups10), groupnames = groups10, maindir = "BAYES/Output", 
  mixvec = KMA2015Strata, 
  catchvec = as.vector(na.omit(as.vector(t(HarvestByStrata2015[sort(KMA2015), ])))), 
  newname = "KMA2015_Annual_Stratified", priorname = '',
  ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1
)
str(KMA2015_Annual_Stratified)
dput(x = KMA2015_Annual_Stratified, file = "Estimates objects/KMA2015_Annual_Stratified.txt")
KMA2015_Annual_Stratified_EstimatesStats <- KMA2015_Annual_Stratified$Stats
dput(x = KMA2015_Annual_Stratified_EstimatesStats, file = "Estimates objects/Final/KMA2015_Annual_Stratified_EstimatesStats.txt")

#~~~~~~~~~~~~~~~~~~
#2016
KMA2016_Annual_Stratified <- StratifiedEstimator.GCL(
  groupvec = seq(groups10), groupnames = groups10, maindir = "BAYES/Output", 
  mixvec = KMA2016Strata, 
  catchvec = as.vector(na.omit(as.vector(t(HarvestByStrata2016[sort(KMA2016), ])))), 
  newname = "KMA2016_Annual_Stratified", priorname = '',
  ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1
)
str(KMA2016_Annual_Stratified)
dput(x = KMA2016_Annual_Stratified, file = "Estimates objects/KMA2016_Annual_Stratified.txt")
KMA2016_Annual_Stratified_EstimatesStats <- KMA2016_Annual_Stratified$Stats
dput(x = KMA2016_Annual_Stratified_EstimatesStats, file = "Estimates objects/Final/KMA2016_Annual_Stratified_EstimatesStats.txt")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Harvest

# Create a list object with all Stratified Annual Harvest
KMA2014_Annual_Stratified_HarvestEstimatesStats <-
  cbind(KMA2014_Annual_Stratified_EstimatesStats[, c("mean", "sd", "median", "5%", "95%")] * sum(HarvestByStrata2014, na.rm = TRUE),
        KMA2014_Annual_Stratified_EstimatesStats[, c("P=0", "GR")])
dput(x = KMA2014_Annual_Stratified_HarvestEstimatesStats, file = "Estimates objects/Final/KMA2014_Annual_Stratified_HarvestEstimatesStats.txt")
str(KMA2014_Annual_Stratified_HarvestEstimatesStats)

KMA2015_Annual_Stratified_HarvestEstimatesStats <-
  cbind(KMA2015_Annual_Stratified_EstimatesStats[, c("mean", "sd", "median", "5%", "95%")] * sum(HarvestByStrata2015, na.rm = TRUE),
        KMA2015_Annual_Stratified_EstimatesStats[, c("P=0", "GR")])
dput(x = KMA2015_Annual_Stratified_HarvestEstimatesStats, file = "Estimates objects/Final/KMA2015_Annual_Stratified_HarvestEstimatesStats.txt")
str(KMA2015_Annual_Stratified_HarvestEstimatesStats)

KMA2016_Annual_Stratified_HarvestEstimatesStats <-
  cbind(KMA2016_Annual_Stratified_EstimatesStats[, c("mean", "sd", "median", "5%", "95%")] * sum(HarvestByStrata2016, na.rm = TRUE),
        KMA2016_Annual_Stratified_EstimatesStats[, c("P=0", "GR")])
dput(x = KMA2016_Annual_Stratified_HarvestEstimatesStats, file = "Estimates objects/Final/KMA2016_Annual_Stratified_HarvestEstimatesStats.txt")
str(KMA2016_Annual_Stratified_HarvestEstimatesStats)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Annual KMA Percentages ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KMA2014_Annual_Stratified_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2014_Annual_Stratified_EstimatesStats.txt")
KMA2015_Annual_Stratified_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2015_Annual_Stratified_EstimatesStats.txt")
KMA2016_Annual_Stratified_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2016_Annual_Stratified_EstimatesStats.txt")



# Three barplot layout
layoutmat <- matrix(data=c(  1, 2,
                             3, 4), nrow = 2, ncol = 2, byrow = TRUE)

ProportionColors <- colorpanel(n = 3, low = "blue", high = "white")

#~~~~~~~~~~~~~~~~~~
# Size Parameters
Groups <- groups10
Groups2Rows <- groups10tworows
cex.lab <- 1.5
cex.xaxis <- 0.5
cex.yaxis <- 1.3
cex.leg <- 1.1
ci.lwd <- 2.5


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make figures as .emf files

# dir.create("Figures/All Years")
require(devEMF)
require(gplots)



emf(file = "Figures/All Years/KMA Proportions 2014-2016.emf", width = 6, height = 5.75, family = "serif", bg = "white")


layout(mat = layoutmat, widths = c(0.075, 1), heights = c(2.7, 0.15))
par(mar = rep(0, 4))
par(family = "times")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Y-axis label
plot.new()
text(x = 0.25, y = 0.5, labels = "Percentage of KMA Harvest", srt = 90, cex = cex.lab)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Barplot
par(mar = c(1, 1, 1.5, 1))
Barplot <- barplot2(height = t(cbind(KMA2014_Annual_Stratified_EstimatesStats[, "median"],
                                   KMA2015_Annual_Stratified_EstimatesStats[, "median"],
                                   KMA2016_Annual_Stratified_EstimatesStats[, "median"])) * 100, 
                    beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                    ci.l = t(cbind(KMA2014_Annual_Stratified_EstimatesStats[, "5%"],
                                 KMA2015_Annual_Stratified_EstimatesStats[, "5%"],
                                 KMA2016_Annual_Stratified_EstimatesStats[, "5%"])) * 100, 
                    ci.u = t(cbind(KMA2014_Annual_Stratified_EstimatesStats[, "95%"],
                                 KMA2015_Annual_Stratified_EstimatesStats[, "95%"],
                                 KMA2016_Annual_Stratified_EstimatesStats[, "95%"])) * 100, 
                    ylim = c(0, 100), col = ProportionColors, yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
legend(legend = 2014:2016, x = "topleft", fill = ProportionColors, border = "black", bty = "n", cex = cex.leg, title="")
abline(h = 0, xpd = FALSE)

mtext(text = Groups2Rows, side = 1, line = 0.66, at = apply(Barplot, 2, mean), adj = 0.5, cex = cex.xaxis)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Blank Corner
par(mar = rep(0, 4))
plot.new()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## x-axis label
par(mar = rep(0, 4))
plot.new()
text(x = 0.5, y = 0.25, labels = "Reporting Group", cex = cex.lab)


dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



emf(file = "Figures/All Years/KMA Proportions 2014-2016 Blank.emf", width = 6, height = 5.75, family = "serif", bg = "white")


layout(mat = layoutmat, widths = c(0.075, 1), heights = c(2.7, 0.15))
par(mar = rep(0, 4))
par(family = "times")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Y-axis label
plot.new()
text(x = 0.25, y = 0.5, labels = "Percentage of KMA Harvest", srt = 90, cex = cex.lab)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Barplot
par(mar = c(1, 1, 1.5, 1))
Barplot <- barplot2(height = matrix(data = 0, nrow = 2, ncol = 10), 
                    beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                    ci.l = matrix(data = 0, nrow = 2, ncol = 10), 
                    ci.u = matrix(data = 0, nrow = 2, ncol = 10), 
                    ylim = c(0, 100), col = ProportionColors, yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
legend(legend = 2014:2016, x = "topleft", fill = ProportionColors, border = "black", bty = "n", cex = cex.leg, title="")
abline(h = 0, xpd = FALSE)

mtext(text = Groups2Rows, side = 1, line = 0.66, at = apply(Barplot, 2, mean), adj = 0.5, cex = cex.xaxis)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Blank Corner
par(mar = rep(0, 4))
plot.new()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## x-axis label
par(mar = rep(0, 4))
plot.new()
text(x = 0.5, y = 0.25, labels = "Reporting Group", cex = cex.lab)


dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Annual KMA Harvest ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KMA2014_Annual_Stratified_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2014_Annual_Stratified_HarvestEstimatesStats.txt")
KMA2015_Annual_Stratified_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2015_Annual_Stratified_HarvestEstimatesStats.txt")
KMA2016_Annual_Stratified_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2016_Annual_Stratified_HarvestEstimatesStats.txt")



# Three barplot layout
layoutmat <- matrix(data=c(  1, 2,
                             3, 4), nrow = 2, ncol = 2, byrow = TRUE)

ProportionColors <- colorpanel(n = 3, low = "green", high = "white")

#~~~~~~~~~~~~~~~~~~
# Size Parameters
Groups <- groups10
Groups2Rows <- groups10tworows
cex.lab <- 1.5
cex.xaxis <- 0.5
cex.yaxis <- 1.3
cex.leg <- 1.1
ci.lwd <- 2.5
ymax <- 4000  # max(sapply(list(KMA2014_Annual_Stratified_HarvestEstimatesStats, KMA2015_Annual_Stratified_HarvestEstimatesStats, KMA2016_Annual_Stratified_HarvestEstimatesStats), function(strata) {strata[, "95%"]}))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make figures as .emf files

# dir.create("Figures/All Years")
require(devEMF)
require(gplots)



emf(file = "Figures/All Years/KMA Harvest 2014-2016.emf", width = 6, height = 5.75, family = "serif", bg = "white")


layout(mat = layoutmat, widths = c(0.075, 1), heights = c(2.7, 0.15))
par(mar = rep(0, 4))
par(family = "times")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Y-axis label
plot.new()
text(x = 0.25, y = 0.5, labels = "Number of Fish Harvested", srt = 90, cex = cex.lab)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Barplot
par(mar = c(1, 1, 1.5, 1))
Barplot <- barplot2(height = t(cbind(KMA2014_Annual_Stratified_HarvestEstimatesStats[, "median"],
                                     KMA2015_Annual_Stratified_HarvestEstimatesStats[, "median"],
                                     KMA2016_Annual_Stratified_HarvestEstimatesStats[, "median"])), 
                    beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                    ci.l = t(cbind(KMA2014_Annual_Stratified_HarvestEstimatesStats[, "5%"],
                                   KMA2015_Annual_Stratified_HarvestEstimatesStats[, "5%"],
                                   KMA2016_Annual_Stratified_HarvestEstimatesStats[, "5%"])), 
                    ci.u = t(cbind(KMA2014_Annual_Stratified_HarvestEstimatesStats[, "95%"],
                                   KMA2015_Annual_Stratified_HarvestEstimatesStats[, "95%"],
                                   KMA2016_Annual_Stratified_HarvestEstimatesStats[, "95%"])), 
                    ylim = c(0, ymax), col = ProportionColors, yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, ymax, 1000), labels = formatC(x = seq(0, ymax, 1000), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
legend(legend = 2014:2016, x = "topleft", fill = ProportionColors, border = "black", bty = "n", cex = cex.leg, title="")
abline(h = 0, xpd = FALSE)

mtext(text = Groups2Rows, side = 1, line = 0.66, at = apply(Barplot, 2, mean), adj = 0.5, cex = cex.xaxis)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Blank Corner
par(mar = rep(0, 4))
plot.new()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## x-axis label
par(mar = rep(0, 4))
plot.new()
text(x = 0.5, y = 0.25, labels = "Reporting Group", cex = cex.lab)


dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Table Annual KMA Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dir.create("Estimates tables")
require(xlsx)

for(yr in 14:16){
  
  EstimatesStats <- dget(file = paste0("Estimates objects/Final/KMA20", yr, "_Annual_Stratified_EstimatesStats.txt"))
  HarvestEstimatesStats <- dget(file = paste0("Estimates objects/Final/KMA20", yr, "_Annual_Stratified_HarvestEstimatesStats.txt"))
  
  Caption <- paste("Table X.-Annual estimates of stock composition (%) and stock-specific harvest for KMA, 20", yr,
                   ". Estimates include median, 90% credibility interval (CI), the probability that the group estimate is equal to zero (P=0), mean, and SD.",
                   sep = '')
  
  
  Disclaimer <- "Note: Stock composition estimates may not sum to 100% and stock-specific harvest estimates may not sum to the total harvest due to rounding error."
  
  
  TableX <- matrix(data = "", nrow = 16, ncol = 13)
  
  TableX[1, 1] <- Caption
  TableX[2, c(2, 9)] <- c("Stock Composition", "Stock-specific Harvest")
  TableX[3, c(3, 10)] <- rep("90% CI", 2)
  TableX[4, c(1, 2:4, 6:7, 9:13, 5)] <- c("Reporting Group", rep(c("Median", "5%", "95%", "Mean", "SD"), 2), "P=0")
  TableX[5:14, 1] <- groups10
  TableX[5:14, c(2:4, 6:7)] <- formatC(x = EstimatesStats[, c("median", "5%", "95%", "mean", "sd")] * 100, digits = 1, format = "f")
  TableX[5:14, 5] <- formatC(x = EstimatesStats[, "P=0"], digits = 2, format = "f")
  TableX[5:14, 9:13] <- formatC(x = HarvestEstimatesStats[, c("median", "5%", "95%", "mean", "sd")], digits = 0, format = "f", big.mark = ",")
  TableX[15, 11:12] <- c("Total", formatC(x = sum(HarvestEstimatesStats[, "mean"]), digits = 0, format = "f", big.mark = ","))
  TableX[16, 1] <- Disclaimer
  
  
  write.xlsx(x = as.data.frame(TableX), 
             file = "Estimates tables/KMA Chinook Estimates Tables.xlsx",
             col.names = FALSE, row.names = FALSE, append = TRUE, sheetName = paste0("KMA20", yr))
}; beep(5)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### GIS Maps ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Which of Andy's 211 do not match up with mine!!!
andy211 <- readClipboard()
which(!sapply(1:211, function(i) {andy211[i] %in% unlist(strsplit(x = KMA211Pops[i], split = "\\."))}))

grep(pattern = andy211[2], x = KMA211Pops)


maporder <- sapply(andy211, function(pop) {
  x <- grep(pattern = pop, x = KMA211Pops)
  ifelse(length(x) == 0, NA, x)
  })

# Confirm that the NAs are still the same ones
maporder[is.na(maporder)] <- which(is.na(maporder))

plot(maporder)

# Confirm that the groupvec doesn't change with map order
table(groups10[groupvec10] == groups10[groupvec10[maporder]])  # it doesn't

# Add colors to GIS map
cbind(groups10, colors10)
t(rbind(groups10, col2rgb(colors10)))
