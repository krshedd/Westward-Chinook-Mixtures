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
QuickBarplot(mixvec = "KMARS15", estimatesstats = KMARS15_EstimatesStats, groups = groups10, groups2rows = groups10tworows, header = KMARS15_Header)

## Make violin plots of posteriors with RGs sorted
ViolinPlot(estimates = KMARS15_Estimates, groups = groups10tworows, colors = colors10, header = KMARS15_Header)
rm(KMARS15_Estimates)

# Are 2014 and 2015 different?
KMARS14vs15 <- compare_comps_between.GCL(mixnames = c("KMARS14", "KMARS15"), groupnames = groups10, mixdir = "BAYES/Output")
str(KMARS14vs15)
