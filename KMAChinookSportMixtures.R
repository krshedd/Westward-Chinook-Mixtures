#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### KMA Chinook Mixtures 2014-2016 ####
# Kyle Shedd Mon Jul 11 15:39:50 2016
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
date()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Introduction ####
# The goal of this script is to analyze Chinook salmon mixtures from the KMA
# sport harvest from 2014-2016 using a coastwide Chinook baseline
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
setwd("V:/Analysis/4_Westward/Chinook/CSRI Westward Sport Harvest 2014-2016/Mixtures")
source("H:/R Source Scripts/Functions.GCL_KS.R")
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")
username <- "krshedd"
password <- "********"


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Pull genotypes from LOKI 2014/2015 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Get collection SILLYs
dir.create("Objects")

KMASport <- c("KMARS14", "KMARS15")
dput(x = KMASport, file = "Objects/KMASport.txt")

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
LOKI2R.GCL(sillyvec = KMASport, username = username, password = password)

rm(username, password)
objects(pattern = "\\.gcl")

## Save unaltered .gcl's as back-up:
dir.create("Raw genotypes")
dir.create("Raw genotypes/OriginalCollections")
invisible(sapply(KMASport, function(silly) {dput(x = get(paste(silly, ".gcl", sep = '')), file = paste("Raw genotypes/OriginalCollections/" , silly, ".txt", sep = ''))} )); beep(8)

## Original sample sizes by SILLY
collection.size.original <- sapply(KMASport, function(silly) get(paste(silly, ".gcl", sep = ""))$n)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Clean workspace; dget .gcl objects and Locus Control ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls(all = TRUE))
setwd("V:/Analysis/4_Westward/Chinook/CSRI Westward Sport Harvest 2014-2016/Mixtures")
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
invisible(sapply(KMASport, function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/OriginalCollections/", silly, ".txt", sep = "")), pos = 1)} )); beep(2)
objects(pattern = "\\.gcl")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Define strata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# No strata, each silly is its own mixture

## All fish have a capture date?
sapply(KMASport, function(silly) {sum(is.na(get(paste(silly, ".gcl", sep = ''))$attributes$CAPTURE_DATE))} )  # Zeros are good


## Confirm that we have other data for these fish?
vials2014 <- as.numeric(readClipboard())
dput(x = vials2014, file= "Objects/vials2014.txt")

# Vials with data, but not genotyped
vials2014[!vials2014 %in% KMARS14.gcl$attributes$FK_FISH_ID]

# Vials genotyped, but with no data
KMARS14.gcl$attributes$FK_FISH_ID[!KMARS14.gcl$attributes$FK_FISH_ID %in% vials2014]
writeClipboard(paste(KMARS14.gcl$attributes$FK_FISH_ID[!KMARS14.gcl$attributes$FK_FISH_ID %in% vials2014], collapse = ", "))



## Confirm that we have other data for these fish?
vials2015 <- as.numeric(readClipboard())
dput(x = vials2015, file= "Objects/vials2015.txt")

# Vials with data, but not genotyped
vials2015[!vials2015 %in% KMARS15.gcl$attributes$FK_FISH_ID]

# Vials genotyped, but with no data
KMARS15.gcl$attributes$FK_FISH_ID[!KMARS15.gcl$attributes$FK_FISH_ID %in% vials2015]



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Data QC/Massage ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
str(KMARS14.gcl)

# Begin QC
require(xlsx)

KMASport

KMASport_SampleSizes <- matrix(data = NA, nrow = length(KMASport), ncol = 4, 
                                         dimnames = list(KMASport, c("Genotyped", "Missing", "Duplicate", "Final")))

#### Check loci
## Get sample size by locus
Original_KMASport_SampleSizebyLocus <- SampSizeByLocus.GCL(sillyvec = KMASport, loci = loci48)
min(Original_KMASport_SampleSizebyLocus)  ## 265
apply(Original_KMASport_SampleSizebyLocus, 1, min) / apply(Original_KMASport_SampleSizebyLocus, 1, max)  # Good, all > 0.9

Original_KMASport_PercentbyLocus <- apply(Original_KMASport_SampleSizebyLocus, 1, function(row) {row / max(row)} )
which(apply(Original_KMASport_PercentbyLocus, 2, min) < 0.8)  # no re-runs!

require(lattice)
new.colors <- colorRampPalette(c("black", "white"))
levelplot(t(Original_KMASport_PercentbyLocus), 
          col.regions = new.colors, 
          at = seq(from = 0, to = 1, length.out = 100), 
          main = "% Genotyped", xlab = "SILLY", ylab = "Locus", 
          scales = list(x = list(rot = 90)), 
          aspect = "fill")  # aspect = "iso" will make squares


#### Check individuals
## View Histogram of Failure Rate by Strata
invisible(sapply(KMASport, function(mix) {
  my.gcl <- get(paste(mix, ".gcl", sep = ''))
  failure <- apply(my.gcl$scores[, , 1], 1, function(ind) {sum(ind == "0") / length(ind)} )
  hist(x = failure, main = mix, xlab = "Failure Rate", col = 8, xlim = c(0, 1), ylim = c(0, 20), breaks = seq(from = 0, to = 1, by = 0.02))
  abline(v = 0.2, lwd = 3)
}))

## Experimental method to look for a way to remove "missing" loci fish based on a strata specific failure rate (rather than 0.8 globally, as done normally)
# Experimented with loess, smooth.split, polynomial regressions, outlier (Franz), and Cooks Distance methods to model outliers.
# Thought of a clever idea to look at the "shoulder" of the cumulative success rate, but this concept is a bit weird.
# Planning to keep it simple!
invisible(sapply(KMASport, function(mix) {
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
Original_KMASport_ColSize <- sapply(paste(KMASport, ".gcl", sep = ''), function(x) get(x)$n)
KMASport_SampleSizes[, "Genotyped"] <- Original_KMASport_ColSize


### Missing
## Remove individuals with >20% missing data
KMASport_MissLoci <- RemoveIndMissLoci.GCL(sillyvec = KMASport, proportion = 0.8)

## Get number of individuals per silly after removing missing loci individuals
ColSize_KMASport_PostMissLoci <- sapply(paste(KMASport, ".gcl", sep = ''), function(x) get(x)$n)
KMASport_SampleSizes[, "Missing"] <- Original_KMASport_ColSize-ColSize_KMASport_PostMissLoci


### Duplicate
## Check within collections for duplicate individuals.
KMASport_DuplicateCheck95MinProportion <- CheckDupWithinSilly.GCL(sillyvec = KMASport, loci = loci48, quantile = NULL, minproportion = 0.95)
KMASport_DuplicateCheckReportSummary <- sapply(KMASport, function(x) KMASport_DuplicateCheck95MinProportion[[x]]$report)
KMASport_DuplicateCheckReportSummary

## Remove duplicate individuals
KMASport_RemovedDups <- RemoveDups.GCL(KMASport_DuplicateCheck95MinProportion)

## Get number of individuals per silly after removing duplicate individuals
ColSize_KMASport_PostDuplicate <- sapply(paste(KMASport, ".gcl", sep = ''), function(x) get(x)$n)
KMASport_SampleSizes[, "Duplicate"] <- ColSize_KMASport_PostMissLoci-ColSize_KMASport_PostDuplicate


### Final
KMASport_SampleSizes[, "Final"] <- ColSize_KMASport_PostDuplicate
KMASport_SampleSizes

dir.create("Output")
write.xlsx(KMASport_SampleSizes, file = "Output/KMASport_SampleSizes.xlsx")
dput(x = KMASport_SampleSizes, file = "Objects/KMASport_SampleSizes.txt")


# dput postQC mixture sillys
dir.create("Raw genotypes/OriginalCollections_PostQC")
invisible(sapply(KMASport, function(silly) {dput(x = get(paste(silly, ".gcl", sep = '')), file = paste("Raw genotypes/OriginalCollections_PostQC/" , silly, ".txt", sep = ''))} )); beep(8)


## Add dates
writeClipboard(as.character(KMARS14.gcl$attributes$FK_FISH_ID))
KMARS14.gcl$attributes$CAPTURE_DATE <- as.Date(as.numeric(readClipboard()), origin = "1899-12-30")
str(KMARS14.gcl)

writeClipboard(as.character(KMARS15.gcl$attributes$FK_FISH_ID))
KMARS15.gcl$attributes$CAPTURE_DATE <- as.Date(as.numeric(readClipboard()), origin = "1899-12-30")
str(KMARS15.gcl)

# dput postQC mixture sillys
dir.create("Raw genotypes/OriginalCollections_PostQC_Dates")
invisible(sapply(KMASport, function(silly) {dput(x = get(paste(silly, ".gcl", sep = '')), file = paste("Raw genotypes/OriginalCollections_PostQC_Dates/" , silly, ".txt", sep = ''))} )); beep(8)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Clean workspace; dget .gcl objects and Locus Control ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls(all = TRUE))
setwd("V:/Analysis/4_Westward/Chinook/CSRI Westward Sport Harvest 2014-2016/Mixtures")
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
invisible(sapply(KMASport, function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/OriginalCollections_PostQC_Dates/", silly, ".txt", sep = "")), pos = 1)} )); beep(2)
objects(pattern = "\\.gcl")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Get/Create MSA Objects ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("V:/Analysis/4_Westward/Chinook/CSRI Westward Commercial Harvest 2014-2016/Mixtures/Objects")

colors10 <- dget(file = "colors10.txt")
groups10 <- dget(file = "groups10.txt")
groups10short <- dget(file = "groups10short.txt")
groupvec10 <- dget(file = "groupvec10.txt")
KMA211Pops <- dget(file = "KMA211Pops.txt")
KMA211Pops42Loci.baseline <- dget(file = "KMA211Pops42Loci.baseline.txt")
loci42 <- dget(file = "loci42.txt")
popnames211 <- dget(file = "popnames211.txt")
KMA211Pops10FlatPrior <- dget(file = "KMA211Pops10FlatPrior.txt")
KMA211PopsInits <- dget(file = "KMA211PopsInits.txt")
KMA211PopsChinookSeeds <- dget(file = "KMA211PopsChinookSeeds.txt")
groups10tworows <- dget(file = "groups10tworows.txt")

setwd("V:/Analysis/4_Westward/Chinook/CSRI Westward Sport Harvest 2014-2016/Mixtures/Objects")

dput(x = colors10, file = "colors10.txt")
dput(x = groups10, file = "groups10.txt")
dput(x = groups10short, file = "groups10short.txt")
dput(x = groupvec10, file = "groupvec10.txt")
dput(x = KMA211Pops, file = "KMA211Pops.txt")
dput(x = KMA211Pops42Loci.baseline, file = "KMA211Pops42Loci.baseline.txt")
dput(x = loci42, file = "loci42.txt")
dput(x = popnames211, file = "popnames211.txt")
dput(x = KMA211Pops10FlatPrior, file = "KMA211Pops10FlatPrior.txt")
dput(x = KMA211PopsInits, file = "KMA211PopsInits.txt")
dput(x = KMA211PopsChinookSeeds, file = "KMA211PopsChinookSeeds.txt")
dput(x = groups10tworows, file = "groups10tworows.txt")

setwd("V:/Analysis/4_Westward/Chinook/CSRI Westward Sport Harvest 2014-2016/Mixtures")


dir.create("BAYES")
sapply(c("Baseline", "Control", "Mixture", "Output"), function(fldr) {dir.create(paste("BAYES", fldr, sep = "/"))})
file.copy(from = "V:/Analysis/4_Westward/Chinook/CSRI Westward Commercial Harvest 2014-2016/Baseline/BAYES/Baseline/KMA211Pops42Loci.bse", 
          to = "BAYES/Baseline/KMA211Pops42Loci.bse")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 1 MSA files for BAYES 2014 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Dumping Mixture files
KMA21142MixtureFormat <- CreateMixture.GCL(sillys = "KMARS14", loci = loci42, IDs = NULL, mixname = "KMARS14",
                                           dir = "BAYES/Mixture", type = "BAYES", PT = FALSE)
dput(KMA21142MixtureFormat, file = "Objects/KMA21142MixtureFormat.txt")


## Dumping Control files
CreateControlFile.GCL(sillyvec = KMA211Pops, loci = loci42, mixname = "KMARS14", basename = "KMA211Pops42Loci", suffix = "", nreps = 40000, nchains = 5,
                      groupvec = groupvec10, priorvec = KMA211Pops10FlatPrior, initmat = KMA211PopsInits, dir = "BAYES/Control",
                      seeds = KMA211PopsChinookSeeds, thin = c(1, 1, 100), mixfortran = KMA21142MixtureFormat, basefortran = KMA211Pops42Loci.baseline, switches = "F T F T T T F")


## Create output directory
dir.create("BAYES/Output/KMARS14")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Go run BAYES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Round 1 Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
KMARS14_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = seq(groups10), groupnames = groups10, 
                                                  maindir = "BAYES/Output", mixvec = "KMARS14", prior = "",  
                                                  ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)

# Dput 1) estimates stats + posterior output & 2) estimates stats
dir.create("Estimates objects")
dput(KMARS14_Estimates, file = "Estimates objects/KMARS14_Estimates.txt")
dput(KMARS14_Estimates$Stats, file = "Estimates objects/KMARS14_EstimatesStats.txt")

KMARS14_Estimates <- dget(file = "Estimates objects/KMARS14_Estimates.txt")
KMARS14_EstimatesStats <- dget(file = "Estimates objects/KMARS14_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2
KMARS14_Estimates$Stats$KMARS14[, "GR"]
table(KMARS14_Estimates$Stats$KMARS14[, "GR"] > 1.2)
sapply("KMARS14", function(Mix) {
  BarPlot <- barplot2(KMARS14_EstimatesStats[[Mix]][, "GR"], col = "blue", ylim = c(1, pmax(1.5, max(KMARS14_EstimatesStats[[Mix]][, "GR"]))), ylab = "Gelman-Rubin", type = "h", xpd = FALSE, main = Mix, names.arg = '')
  abline(h = 1.2, lwd = 3, xpd = FALSE)
  text(x = BarPlot, y = 1, labels = groups10tworows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
})

# Quick look at raw posterior output
str(KMARS14_Estimates$Output)
KMARS14_Header <- setNames(object = c("Kodiak Sport April 16-August 29, 2014"), 
                                       nm = "KMARS14")
dput(x = KMARS14_Header, file = "Objects/KMARS14_Header.txt")

file.copy(from = "V:/Analysis/4_Westward/Chinook/CSRI Westward Commercial Harvest 2014-2016/Mixtures/Objects/PlotPosterior.txt",
          to = "Objects/PlotPosterior.txt")

PlotPosterior(mixvec = "KMARS14", output = KMARS14_Estimates$Output, 
              groups = groups10, colors = colors10, 
              header = KMARS14_Header, set.mfrow = c(5, 2), thin = 10)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Round 1 Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Barplots
file.copy(from = "V:/Analysis/4_Westward/Chinook/CSRI Westward Commercial Harvest 2014-2016/Mixtures/Objects/QuickBarplot.txt",
          to = "Objects/QuickBarplot.txt")
KMARS14_EstimatesStats <- dget(file = "Estimates objects/KMARS14_EstimatesStats.txt")
QuickBarplot(mixvec = "KMARS14", estimatesstats = KMARS14_EstimatesStats, groups = groups10, groups2rows = groups10tworows, header = KMARS14_Header)


## Make violin plots of posteriors with RGs sorted

file.copy(from = "V:/Analysis/4_Westward/Chinook/CSRI Westward Commercial Harvest 2014-2016/Mixtures/Objects/ViolinPlot.txt",
          to = "Objects/ViolinPlot.txt")

ViolinPlot(estimates = KMARS14_Estimates, groups = groups10tworows, colors = colors10, header = KMARS14_Header)
rm(KMARS14_Estimates)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 2 MSA files for BAYES 2015 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create rolling prior based on 2014 Round 1 estimates
KMARS14_EstimatesStats <- dget(file = "Estimates objects/KMARS14_EstimatesStats.txt")

KMARS15_Prior <- sapply(KMARS14_EstimatesStats, function(Mix) {
  Prior.GCL(groupvec = groupvec10, groupweights = Mix[, "mean"], minval = 0.01)}, simplify = FALSE)
names(KMARS15_Prior) <- gsub(pattern = "S14", replacement = "S15", x = names(KMARS15_Prior))  # This changes the names
dput(x = KMARS15_Prior, file = "Objects/KMARS15_Prior.txt")
str(KMARS15_Prior)

# Verify
plot(as.vector(KMARS15_Prior[["KMARS15"]]), type = "h", main = "KMARS15")

## Dumping Mixture files
CreateMixture.GCL(sillys = "KMARS15", loci = loci42, IDs = NULL, mixname = "KMARS15", dir = "BAYES/Mixture", type = "BAYES", PT = FALSE)

## Dumping Control files
CreateControlFile.GCL(sillyvec = KMA211Pops, loci = loci42, mixname = "KMARS15", basename = "KMA211Pops42Loci", suffix = "", nreps = 40000, nchains = 5,
                      groupvec = groupvec10, priorvec = KMARS15_Prior[["KMARS15"]], initmat = KMA211PopsInits, dir = "BAYES/Control",
                      seeds = KMA211PopsChinookSeeds, thin = c(1, 1, 100), mixfortran = KMA21142MixtureFormat, basefortran = KMA211Pops42Loci.baseline, switches = "F T F T T T F")


## Create output directory
dir.create("BAYES/Output/KMARS15")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Go run BAYES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Round 2 Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
KMARS15_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = seq(groups10), groupnames = groups10, 
                                                  maindir = "BAYES/Output", mixvec = "KMARS15", prior = "",  
                                                  ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)

# Dput 1) estimates stats + posterior output & 2) estimates stats
dput(KMARS15_Estimates, file = "Estimates objects/KMARS15_Estimates.txt")
dput(KMARS15_Estimates$Stats, file = "Estimates objects/KMARS15_EstimatesStats.txt")

KMARS15_Estimates <- dget(file = "Estimates objects/KMARS15_Estimates.txt")
KMARS15_EstimatesStats <- dget(file = "Estimates objects/KMARS15_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2
KMARS15_Estimates$Stats$KMARS15[, "GR"]
table(KMARS15_Estimates$Stats$KMARS15[, "GR"] > 1.2)
sapply("KMARS15", function(Mix) {
  BarPlot <- barplot2(KMARS15_EstimatesStats[[Mix]][, "GR"], col = "blue", ylim = c(1, pmax(1.5, max(KMARS15_EstimatesStats[[Mix]][, "GR"]))), ylab = "Gelman-Rubin", type = "h", xpd = FALSE, main = Mix, names.arg = '')
  abline(h = 1.2, lwd = 3, xpd = FALSE)
  text(x = BarPlot, y = 1, labels = groups10tworows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
})

# Quick look at raw posterior output
str(KMARS15_Estimates$Output)
KMARS15_Header <- setNames(object = c("Kodiak Sport May 17-August 14, 2015"), 
                           nm = "KMARS15")
dput(x = KMARS15_Header, file = "Objects/KMARS15_Header.txt")

PlotPosterior(mixvec = "KMARS15", output = KMARS15_Estimates$Output, 
              groups = groups10, colors = colors10, 
              header = KMARS15_Header, set.mfrow = c(5, 2), thin = 10)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Round 2 Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Barplots
KMARS15_EstimatesStats <- dget(file = "Estimates objects/KMARS15_EstimatesStats.txt")
QuickBarplot(mixvec = "KMARS15", estimatesstats = KMARS15_EstimatesStats, groups = groups10, groups2rows = groups10tworows, header = KMARS15_Header)

## Make violin plots of posteriors with RGs sorted
ViolinPlot(estimates = KMARS15_Estimates, groups = groups10tworows, colors = colors10, header = KMARS15_Header)
rm(KMARS15_Estimates)

# Are 2014 and 2015 different?
KMARS14vs15 <- compare_comps_between.GCL(mixnames = c("KMARS14", "KMARS15"), groupnames = groups10, mixdir = "BAYES/Output")
str(KMARS14vs15)







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Clean workspace; dget .gcl objects and Locus Control ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls(all = TRUE))
setwd("V:/Analysis/4_Westward/Chinook/CSRI Westward Sport Harvest 2014-2016/Mixtures")
# This sources all of the new GCL functions to this workspace
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")
source("H:/R Source Scripts/Functions.GCL_KS.R")

## Get objects
LocusControl <- dget(file = "Objects/OriginalLocusControl48.txt")

KMAobjects <- list.files(path = "Objects", recursive = FALSE)
KMAobjects <- KMAobjects[!KMAobjects %in% c("OriginalLocusControl48.txt")]
KMAobjects

invisible(sapply(KMAobjects, function(objct) {assign(x = unlist(strsplit(x = objct, split = ".txt")), value = dget(file = paste(getwd(), "Objects", objct, sep = "/")), pos = 1) })); beep(2)


## Get collection SILLYs
KMASport <- c("KMARS14", "KMARS15", "KMARS16")
dput(x = KMASport, file = "Objects/KMASport.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Pull genotypes from LOKI 2016 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Pull all data for each silly code and create .gcl objects for each
LOKI2R.GCL(sillyvec = "KMARS16", username = username, password = password)

rm(username, password)
objects(pattern = "\\.gcl")

## Save unaltered .gcl's as back-up:
invisible(sapply("KMARS16", function(silly) {dput(x = get(paste(silly, ".gcl", sep = '')), file = paste("Raw genotypes/OriginalCollections/" , silly, ".txt", sep = ''))} )); beep(8)

## Original sample sizes by SILLY
KMARS16.gcl$n

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Define strata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# No strata, each silly is its own mixture

## All fish have a capture date?
sum(is.na(KMARS16.gcl$attributes$CAPTURE_DATE))  # Zeros are good


## Confirm that we have other data for these fish?
vials2016 <- as.numeric(readClipboard())
dput(x = vials2016, file= "Objects/vials2016.txt")

# Vials with data, but not genotyped
vials2016[!vials2016 %in% KMARS16.gcl$attributes$FK_FISH_ID]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Data QC/Massage ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
str(KMARS16.gcl)

# Begin QC
require(xlsx)

KMASport16 <- "KMARS16"

KMASport16_SampleSizes <- matrix(data = NA, nrow = length(KMASport16), ncol = 4, 
                               dimnames = list(KMASport16, c("Genotyped", "Missing", "Duplicate", "Final")))

#### Check loci
## Get sample size by locus
Original_KMASport16_SampleSizebyLocus <- SampSizeByLocus.GCL(sillyvec = KMASport16, loci = loci48)
min(Original_KMASport16_SampleSizebyLocus)  ## 427
apply(Original_KMASport16_SampleSizebyLocus, 1, min) / apply(Original_KMASport16_SampleSizebyLocus, 1, max)  # Good, all > 0.9

Original_KMASport16_PercentbyLocus <- apply(Original_KMASport16_SampleSizebyLocus, 1, function(row) {row / max(row)} )
which(apply(Original_KMASport16_PercentbyLocus, 2, min) < 0.8)  # no re-runs!

require(lattice)
new.colors <- colorRampPalette(c("black", "white"))
levelplot(t(Original_KMASport16_PercentbyLocus), 
          col.regions = new.colors, 
          at = seq(from = 0, to = 1, length.out = 100), 
          main = "% Genotyped", xlab = "SILLY", ylab = "Locus", 
          scales = list(x = list(rot = 90)), 
          aspect = "fill")  # aspect = "iso" will make squares


#### Check individuals
## View Histogram of Failure Rate by Strata
invisible(sapply(KMASport16, function(mix) {
  my.gcl <- get(paste(mix, ".gcl", sep = ''))
  failure <- apply(my.gcl$scores[, , 1], 1, function(ind) {sum(ind == "0") / length(ind)} )
  hist(x = failure, main = mix, xlab = "Failure Rate", col = 8, xlim = c(0, 1), ylim = c(0, 20), breaks = seq(from = 0, to = 1, by = 0.02))
  abline(v = 0.2, lwd = 3)
}))


### Initial
## Get number of individuals per silly before removing missing loci individuals
Original_KMASport16_ColSize <- sapply(paste(KMASport16, ".gcl", sep = ''), function(x) get(x)$n)
KMASport16_SampleSizes[, "Genotyped"] <- Original_KMASport16_ColSize


### Missing
## Remove individuals with >20% missing data
KMASport16_MissLoci <- RemoveIndMissLoci.GCL(sillyvec = KMASport16, proportion = 0.8)

## Get number of individuals per silly after removing missing loci individuals
ColSize_KMASport16_PostMissLoci <- sapply(paste(KMASport16, ".gcl", sep = ''), function(x) get(x)$n)
KMASport16_SampleSizes[, "Missing"] <- Original_KMASport16_ColSize-ColSize_KMASport16_PostMissLoci


### Duplicate
## Check within collections for duplicate individuals.
KMASport16_DuplicateCheck95MinProportion <- CheckDupWithinSilly.GCL(sillyvec = KMASport16, loci = loci48, quantile = NULL, minproportion = 0.95)
KMASport16_DuplicateCheckReportSummary <- sapply(KMASport16, function(x) KMASport16_DuplicateCheck95MinProportion[[x]]$report)
KMASport16_DuplicateCheckReportSummary

## Remove duplicate individuals
KMASport16_RemovedDups <- RemoveDups.GCL(KMASport16_DuplicateCheck95MinProportion)

## Get number of individuals per silly after removing duplicate individuals
ColSize_KMASport16_PostDuplicate <- sapply(paste(KMASport16, ".gcl", sep = ''), function(x) get(x)$n)
KMASport16_SampleSizes[, "Duplicate"] <- ColSize_KMASport16_PostMissLoci-ColSize_KMASport16_PostDuplicate


### Final
KMASport16_SampleSizes[, "Final"] <- ColSize_KMASport16_PostDuplicate
KMASport16_SampleSizes


write.xlsx(KMASport16_SampleSizes, file = "Output/KMASport16_SampleSizes.xlsx")
dput(x = KMASport16_SampleSizes, file = "Objects/KMASport16_SampleSizes.txt")


# Create final QC matrix with all three years
KMASport_SampleSizes <- rbind(KMASport_SampleSizes, KMASport16_SampleSizes)
write.xlsx(KMASport_SampleSizes, file = "Output/KMASport_SampleSizes.xlsx")
dput(x = KMASport_SampleSizes, file = "Objects/KMASport_SampleSizes.txt")


# dput postQC mixture sillys
invisible(sapply(KMASport16, function(silly) {dput(x = get(paste(silly, ".gcl", sep = '')), file = paste("Raw genotypes/OriginalCollections_PostQC/" , silly, ".txt", sep = ''))} )); beep(8)


## Add dates
writeClipboard(as.character(KMARS16.gcl$attributes$FK_FISH_ID))
KMARS16.gcl$attributes$CAPTURE_DATE <- as.Date(as.numeric(readClipboard()), origin = "1899-12-30")
str(KMARS16.gcl)

# dput postQC mixture sillys
invisible(sapply(KMASport16, function(silly) {dput(x = get(paste(silly, ".gcl", sep = '')), file = paste("Raw genotypes/OriginalCollections_PostQC_Dates/" , silly, ".txt", sep = ''))} )); beep(8)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Clean workspace; dget .gcl objects and Locus Control ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls(all = TRUE))
setwd("V:/Analysis/4_Westward/Chinook/CSRI Westward Sport Harvest 2014-2016/Mixtures")
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
invisible(sapply(KMASport, function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/OriginalCollections_PostQC_Dates/", silly, ".txt", sep = "")), pos = 1)} )); beep(2)
objects(pattern = "\\.gcl")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 3 MSA files for BAYES 2016 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create rolling prior based on 2015 Round 1 estimates
KMARS15_EstimatesStats <- dget(file = "Estimates objects/KMARS15_EstimatesStats.txt")

KMARS16_Prior <- sapply(KMARS15_EstimatesStats, function(Mix) {
  Prior.GCL(groupvec = groupvec10, groupweights = Mix[, "mean"], minval = 0.01)}, simplify = FALSE)
names(KMARS16_Prior) <- gsub(pattern = "S15", replacement = "S16", x = names(KMARS16_Prior))  # This changes the names
dput(x = KMARS16_Prior, file = "Objects/KMARS16_Prior.txt")
str(KMARS16_Prior)

# Verify
plot(as.vector(KMARS16_Prior[["KMARS16"]]), type = "h", main = "KMARS16")

## Dumping Mixture files
CreateMixture.GCL(sillys = "KMARS16", loci = loci42, IDs = NULL, mixname = "KMARS16", dir = "BAYES/Mixture", type = "BAYES", PT = FALSE)

## Dumping Control files
CreateControlFile.GCL(sillyvec = KMA211Pops, loci = loci42, mixname = "KMARS16", basename = "KMA211Pops42Loci", suffix = "", nreps = 40000, nchains = 5,
                      groupvec = groupvec10, priorvec = KMARS16_Prior[["KMARS16"]], initmat = KMA211PopsInits, dir = "BAYES/Control",
                      seeds = KMA211PopsChinookSeeds, thin = c(1, 1, 100), mixfortran = KMA21142MixtureFormat, basefortran = KMA211Pops42Loci.baseline, switches = "F T F T T T F")


## Create output directory
dir.create("BAYES/Output/KMARS16")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Go run BAYES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Round 3 Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
KMARS16_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = seq(groups10), groupnames = groups10, 
                                                  maindir = "BAYES/Output", mixvec = "KMARS16", prior = "",  
                                                  ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)

# Dput 1) estimates stats + posterior output & 2) estimates stats
dput(KMARS16_Estimates, file = "Estimates objects/KMARS16_Estimates.txt")
dput(KMARS16_Estimates$Stats, file = "Estimates objects/KMARS16_EstimatesStats.txt")

KMARS16_Estimates <- dget(file = "Estimates objects/KMARS16_Estimates.txt")
KMARS16_EstimatesStats <- dget(file = "Estimates objects/KMARS16_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2
KMARS16_Estimates$Stats$KMARS16[, "GR"]
table(KMARS16_Estimates$Stats$KMARS16[, "GR"] > 1.2)
sapply("KMARS16", function(Mix) {
  BarPlot <- barplot2(KMARS16_EstimatesStats[[Mix]][, "GR"], col = "blue", ylim = c(1, pmax(1.5, max(KMARS16_EstimatesStats[[Mix]][, "GR"]))), ylab = "Gelman-Rubin", type = "h", xpd = FALSE, main = Mix, names.arg = '')
  abline(h = 1.2, lwd = 3, xpd = FALSE)
  text(x = BarPlot, y = 1, labels = groups10tworows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
})

# Quick look at raw posterior output
str(KMARS16_Estimates$Output)
KMARS16_Header <- setNames(object = c("Kodiak Sport May 22-August 13, 2016"), 
                           nm = "KMARS16")
dput(x = KMARS16_Header, file = "Objects/KMARS16_Header.txt")

PlotPosterior(mixvec = "KMARS16", output = KMARS16_Estimates$Output, 
              groups = groups10, colors = colors10, 
              header = KMARS16_Header, set.mfrow = c(5, 2), thin = 10)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Round 3 Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Barplots
KMARS16_EstimatesStats <- dget(file = "Estimates objects/KMARS16_EstimatesStats.txt")
QuickBarplot(mixvec = "KMARS16", estimatesstats = KMARS16_EstimatesStats, groups = groups10, groups2rows = groups10tworows, header = KMARS16_Header)

## Make violin plots of posteriors with RGs sorted
ViolinPlot(estimates = KMARS16_Estimates, groups = groups10tworows, colors = colors10, header = KMARS16_Header)
rm(KMARS16_Estimates)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Re-do Round 3 MSA files for BAYES 2016 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Cook Inlet Gelman-Rubin is over 1.2
table(KMARS16_Estimates$Stats$KMARS16[, "GR"] > 1.2)


## Dumping Control files
CreateControlFile.GCL(sillyvec = KMA211Pops, loci = loci42, mixname = "KMARS16", basename = "KMA211Pops42Loci", suffix = "", nreps = 80000, nchains = 5,
                      groupvec = groupvec10, priorvec = KMARS16_Prior[["KMARS16"]], initmat = KMA211PopsInits, dir = "BAYES/Control",
                      seeds = KMA211PopsChinookSeeds, thin = c(1, 1, 100), mixfortran = KMA21142MixtureFormat, basefortran = KMA211Pops42Loci.baseline, switches = "F T F T T T F")



KMARS16_80K_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = seq(groups10), groupnames = groups10, 
                                                  maindir = "BAYES/Output/80K/", mixvec = "KMARS16", prior = "",  
                                                  ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)

# Dput 1) estimates stats + posterior output & 2) estimates stats
dput(KMARS16_80K_Estimates, file = "Estimates objects/KMARS16_80K_Estimates.txt")
dput(KMARS16_80K_Estimates$Stats, file = "Estimates objects/KMARS16_80K_EstimatesStats.txt")

KMARS16_80K_Estimates <- dget(file = "Estimates objects/KMARS16_80K_Estimates.txt")
KMARS16_80K_EstimatesStats <- dget(file = "Estimates objects/KMARS16_80K_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2
KMARS16_80K_Estimates$Stats$KMARS16[, "GR"]
table(KMARS16_80K_Estimates$Stats$KMARS16[, "GR"] > 1.2)
sapply("KMARS16", function(Mix) {
  BarPlot <- barplot2(KMARS16_80K_EstimatesStats[[Mix]][, "GR"], col = "blue", ylim = c(1, pmax(1.5, max(KMARS16_80K_EstimatesStats[[Mix]][, "GR"]))), ylab = "Gelman-Rubin", type = "h", xpd = FALSE, main = Mix, names.arg = '')
  abline(h = 1.2, lwd = 3, xpd = FALSE)
  text(x = BarPlot, y = 1, labels = groups10tworows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
})

# Quick look at raw posterior output
str(KMARS16_80K_Estimates$Output)
KMARS16_Header <- setNames(object = c("Kodiak Sport May 22-August 13, 2016"), 
                           nm = "KMARS16")
dput(x = KMARS16_Header, file = "Objects/KMARS16_Header.txt")

PlotPosterior(mixvec = "KMARS16", output = KMARS16_80K_Estimates$Output, 
              groups = groups10, colors = colors10, 
              header = KMARS16_Header, set.mfrow = c(5, 2), thin = 10)




# Plot Annual vs. Early Means
KMARS16_EstimatesStats <- dget(file = "Estimates objects/KMARS16_EstimatesStats.txt")

par(mar = c(4.1, 5.1, 3.1, 1.1))
Barplot <- barplot2(height = rbind(KMARS16_EstimatesStats[[1]][, "median"], KMARS16_80K_EstimatesStats[[1]][, "median"]) * 100, plot.ci = TRUE,
                    ci.l = rbind(KMARS16_EstimatesStats[[1]][, "5%"], KMARS16_80K_EstimatesStats[[1]][, "5%"]) * 100,
                    ci.u = rbind(KMARS16_EstimatesStats[[1]][, "95%"], KMARS16_80K_EstimatesStats[[1]][, "95%"]) * 100,
                    beside = TRUE, col = c("blue", "white"), yaxt = 'n', xaxt = 'n', main = "KMARS16", ylab = "Precent of Mixture",
                    cex.lab = 2, cex.main = 2, ylim = c(0, 100))
axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = ",", digits = 0, format = "f"), cex.axis = 1.5)
abline(h = 0, xpd = FALSE)
text(x = colMeans(Barplot), y = -1, labels = groups10tworows, srt = 90, adj = 1, xpd = TRUE, cex = 0.6)
legend("topleft", legend = c("40K", "80K"), fill = c("blue", "white"), bty = 'n', cex = 1.5)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Percentages for KMA Mixtures ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KMARS14_EstimatesStats <- dget(file = "Estimates objects/KMARS14_EstimatesStats.txt")
KMARS15_EstimatesStats <- dget(file = "Estimates objects/KMARS15_EstimatesStats.txt")
KMARS16_EstimatesStats <- dget(file = "Estimates objects/KMARS16_80K_EstimatesStats.txt")

str(KMARS14_EstimatesStats)

# Three barplot layout
layoutmat <- matrix(data=c(  1, 2,
                             1, 3,
                             1, 4,
                             5, 6), nrow = 4, ncol = 2, byrow = TRUE)

#~~~~~~~~~~~~~~~~~~
# Size Parameters
Groups <- groups10
Groups2Rows <- groups10tworows
cex.lab <- 1.5
cex.xaxis <- 0.7
cex.yaxis <- 1.3
cex.leg <- 1.1
ci.lwd <- 2.5


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make figures as .emf files

# dir.create("Figures")
require(devEMF)
require(gplots)


emf(file = paste("Figures/KMA Sport Proportions 2014-2016.emf", sep = ''), width = 6, height = 5.75, family = "serif", bg = "white")


layout(mat = layoutmat, widths = c(0.075, 1, 1), heights = c(0.9, 0.9, 0.9, 0.15))
par(mar = rep(0, 4))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Y-axis label
plot.new()
par(family = "times")
text(x = 0.25, y = 0.5, labels = "Percentage of Catch", srt = 90, cex = cex.lab)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2014 Barplot
tempmix <- "KMARS14"
par(mar = c(1, 1, 1, 1))
Barplot14 <- barplot2(height = KMARS14_EstimatesStats[[tempmix]][, "median"] * 100, 
                      beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                      ci.l = KMARS14_EstimatesStats[[tempmix]][, "5%"] * 100, 
                      ci.u = KMARS14_EstimatesStats[[tempmix]][, "95%"] * 100, 
                      ylim = c(0, 100), col = "blue", yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
legend(legend = "April 16-August 29", x = "topleft", fill = "blue", border = "black", bty = "n", cex = cex.leg, title="2014")
abline(h = 0, xpd = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2015 Barplot
tempmix <- "KMARS15"
par(mar = c(1, 1, 1, 1))
Barplot15 <- barplot2(height = KMARS15_EstimatesStats[[tempmix]][, "median"] * 100, 
                      beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                      ci.l = KMARS15_EstimatesStats[[tempmix]][, "5%"] * 100, 
                      ci.u = KMARS15_EstimatesStats[[tempmix]][, "95%"] * 100, 
                      ylim = c(0, 100), col = "blue", yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
legend(legend = "May 17-August 14", x = "topleft", fill = "blue", border = "black", bty = "n", cex = cex.leg, title="2015")
abline(h = 0, xpd = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2016 Barplot
tempmix <- "KMARS16"
par(mar = c(1, 1, 1, 1))
Barplot16 <- barplot2(height = KMARS16_EstimatesStats[[tempmix]][, "median"] * 100, 
                      beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                      ci.l = KMARS16_EstimatesStats[[tempmix]][, "5%"] * 100, 
                      ci.u = KMARS16_EstimatesStats[[tempmix]][, "95%"] * 100, 
                      ylim = c(0, 100), col = "blue", yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
legend(legend = "May 22-August 13", x = "topleft", fill = "blue", border = "black", bty = "n", cex = cex.leg, title="2016")
abline(h = 0, xpd = FALSE)

mtext(text = Groups2Rows, side = 1, line = 1.1, at = Barplot16, adj = 0.5, cex = cex.xaxis)

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
#### Table results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require(xlsx)
# dir.create("Estimates tables")


# 2014
Caption <- paste0("Table X.-Annual estimates of stock composition (%) for the Kodiak Area Sport Fishery, April 16-August 29, 2014. Estimates include median, 90% credibility interval (CI), the probability that the group estimate is equal to zero (P=0), mean, and SD.")

Disclaimer <- "Note: Stock composition estimates may not sum to 100% due to rounding error."


TableX <- matrix(data = "", nrow = 16, ncol = 7)

TableX[1, 1] <- Caption
TableX[2, 2] <- "Stock Composition"
TableX[3, 3] <- "90% CI"
TableX[4, c(1, 2:4, 6:7, 5)] <- c("Reporting Group", c("Median", "5%", "95%", "Mean", "SD"), "P=0")
TableX[5:14, 1] <- groups10
TableX[5:14, c(2:4, 6:7)] <- formatC(x = KMARS14_EstimatesStats[["KMARS14"]][, c("median", "5%", "95%", "mean", "sd")] * 100, digits = 1, format = "f")
TableX[5:14, 5] <- formatC(x = KMARS14_EstimatesStats[["KMARS14"]][, "P=0"], digits = 2, format = "f")
TableX[16, 1] <- Disclaimer


write.xlsx(x = as.data.frame(TableX), 
           file = "Estimates tables/KMA Chinook Estimates Tables.xlsx",
           col.names = FALSE, row.names = FALSE, append = TRUE, sheetName = "KMARS14")



# 2015
Caption <- paste0("Table X.-Annual estimates of stock composition (%) for the Kodiak Area Sport Fishery, May 17-August 14, 2015. Estimates include median, 90% credibility interval (CI), the probability that the group estimate is equal to zero (P=0), mean, and SD.")

Disclaimer <- "Note: Stock composition estimates may not sum to 100% due to rounding error."


TableX <- matrix(data = "", nrow = 16, ncol = 7)

TableX[1, 1] <- Caption
TableX[2, 2] <- "Stock Composition"
TableX[3, 3] <- "90% CI"
TableX[4, c(1, 2:4, 6:7, 5)] <- c("Reporting Group", c("Median", "5%", "95%", "Mean", "SD"), "P=0")
TableX[5:14, 1] <- groups10
TableX[5:14, c(2:4, 6:7)] <- formatC(x = KMARS15_EstimatesStats[["KMARS15"]][, c("median", "5%", "95%", "mean", "sd")] * 100, digits = 1, format = "f")
TableX[5:14, 5] <- formatC(x = KMARS15_EstimatesStats[["KMARS15"]][, "P=0"], digits = 2, format = "f")
TableX[16, 1] <- Disclaimer


write.xlsx(x = as.data.frame(TableX), 
           file = "Estimates tables/KMA Chinook Estimates Tables.xlsx",
           col.names = FALSE, row.names = FALSE, append = TRUE, sheetName = "KMARS15")



# 2016
Caption <- paste0("Table X.-Annual estimates of stock composition (%) for the Kodiak Area Sport Fishery, May 22-August 13, 2016. Estimates include median, 90% credibility interval (CI), the probability that the group estimate is equal to zero (P=0), mean, and SD.")

Disclaimer <- "Note: Stock composition estimates may not sum to 100% due to rounding error."


TableX <- matrix(data = "", nrow = 16, ncol = 7)

TableX[1, 1] <- Caption
TableX[2, 2] <- "Stock Composition"
TableX[3, 3] <- "90% CI"
TableX[4, c(1, 2:4, 6:7, 5)] <- c("Reporting Group", c("Median", "5%", "95%", "Mean", "SD"), "P=0")
TableX[5:14, 1] <- groups10
TableX[5:14, c(2:4, 6:7)] <- formatC(x = KMARS16_EstimatesStats[["KMARS16"]][, c("median", "5%", "95%", "mean", "sd")] * 100, digits = 1, format = "f")
TableX[5:14, 5] <- formatC(x = KMARS16_EstimatesStats[["KMARS16"]][, "P=0"], digits = 2, format = "f")
TableX[16, 1] <- Disclaimer


write.xlsx(x = as.data.frame(TableX), 
           file = "Estimates tables/KMA Chinook Estimates Tables.xlsx",
           col.names = FALSE, row.names = FALSE, append = TRUE, sheetName = "KMARS16")
