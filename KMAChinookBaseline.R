#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### KMA Chinook Baseline 2014-2016 ####
# Kyle Shedd Mon Jul 11 15:39:50 2016
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
date()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Introduction ####
# The goal of this script is copy Andy's Lower Cook Inlet baseline objects and
# create a new groupvec to specify the reporting groups from the op plan (see
# below). After loading Andy's workspace, I will save the relevant objects,
# create the appropriate groupvec, and quickly look at the Likelihood Profile
# to verify that this baseline will work for us.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Specific Objectives ####
# This script will:
# 1) Import Andy's Lower Cook Inlet baseline data
# 2) Save relevant objects
# 3) Find a common markerset
# 4) Create new groupvec
# 5) Add new collections to existing populations (sample size)
# 6) Likelihood profile + MDS to view structure
# 7) Tree for report
# 7) Dump .bse


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

# Load Andy's workspace
load("V:/Analysis/2_Central/Chinook/Lower Cook Inlet/2015/Baseline/LowerCIChinook2015Baseline.RData")


# rm(list = ls(all = TRUE))
setwd("V:/Analysis/4_Westward/Chinook/CSRI Westward Commercial Harvest 2014-2016/Baseline")
source("H:/R Source Scripts/Functions.GCL_KS.R")
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")


## Create Directories
# dirs <- c("Objects", "Likelihood Profiles", "Raw gentoypes", "BAYES")
# sapply(dirs, dir.create)
# dir.create(path = "BAYES/Baseline")

## Move BAYES baseline .bse file
# file.copy(from = "V:/Analysis/2_Central/Chinook/Lower Cook Inlet/2015/Baseline/BAYES/Baseline/LCI211pops43loci.bse", 
#           to = "BAYES/Baseline/LCI211pops43loci.bse")

## Save objects
dput(x = LocusControl, file = "Objects/OriginalLocusControl_loci52.txt")
dput(x = loci52, file = "Objects/loci52.txt")
dput(x = loci43, file = "Objects/loci43.txt")
dput(x = PooledNames211, file = "Objects/PooledNames211.txt")
# dput(x = PooledNames211ordered, file = "Objects/PooledNames211ordered.txt")  # Don't know why this is here, as Andy did NOT use it
dput(x = popnames211, file = "Objects/popnames211.txt")
dput(x = groups14, file = "Objects/groups14.txt")
dput(x = groupvec14, file = "Objects/groupvec14.txt")
dput(x = colors14, file = "Objects/colors14.txt")

## Confirm how baseline was made
# basefortran=CreateBaseline.GCL(sillyvec=PooledNames211,loci=loci43,dir="BAYES/Baseline",basename="LCI211pops43loci",type="BAYES",groupvec=NULL)
# CreateLocusControl.GCL(markersuite = "Chinook2011_52SNPs", username = username, password = password)
# loci43 <- loci52[-c(7,11,13,21,26,29,33,49,52)]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Locus Control for Mixtures ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Read in data
username <- "krshedd"
password <- "********"

rm(LocusControl)
CreateLocusControl.GCL(markersuite = "Chinook_Kodiak_2016_48SNPs", username = username, password = password)
loci48 <- LocusControl$locusnames

dput(LocusControl, file = "Objects/OriginalLocusControl_loci48.txt")
dput(loci48, file = "Objects/loci48.txt")

## Save unaltered .gcl's as back-up:
dir.create("Raw genotypes/PostQCPooledPops")
invisible(sapply(PooledNames211, function(silly) {dput(x = get(paste(silly, ".gcl", sep = '')), file = paste("Raw genotypes/PostQCPooledPops/" , silly, ".txt", sep = ''))} )); beep(8)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Clean workspace; dget .gcl objects and Locus Control ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls(all = TRUE))
setwd("V:/Analysis/4_Westward/Chinook/CSRI Westward Commercial Harvest 2014-2016/Baseline")
source("H:/R Source Scripts/Functions.GCL_KS.R")
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Get objects
baselineobjects <- sapply(list.files(path = "Objects"), function(file) {unlist(strsplit(x = file, split = ".txt"))} , USE.NAMES = FALSE)
invisible(sapply(baselineobjects, function(file) {assign(x = file, value = dget(file = paste("Objects/", file, ".txt", sep = "")), pos = 1)} ))
LocusControl <- OriginalLocusControl_loci48; rm(OriginalLocusControl_loci48, OriginalLocusControl_loci52)

## Get pops
require(beepr)
invisible(sapply(PooledNames211, function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/PostQCPooledPops/", silly, ".txt", sep = "")), pos = 1)} )); beep(2)
objects(pattern = "\\.gcl")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Final Markerset ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### There are 4 markersets
# loci48 (what mixtures were genotyped for)
loci48

# loci43 (inheretted from Andy's baseline)
loci43

# loci42 (what is in common between loci48 and loci43)
loci42 <- loci43[loci43 %in% loci48]
dput(loci42, file = "Objects/loci42.txt")

# 45 loci templin (what Templin et al. 2011 used)
x <- readClipboard()  # from Table 1 Templin et al. 2011 Genetic differentiation of Alaska Chinook salmon: the missing link for migratory studies
x <- strsplit(x = x, split = " ")
templin.loci <- sapply(x, function(row) {row[1]})
templin.loci <- templin.loci[-length(templin.loci)]  # remove "Overall"
templin.loci <- gsub(pattern = "S71", replacement = "S7-1", x = templin.loci)  # rename Ots_S71 to Ots_S7-1
dput(x = templin.loci, file = "Objects/templin.loci.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Should make sure to pool, etc. with loci48
# Should do any sort of testing with the final markerset
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## Extra markers that mixture samples were genotyped for not in baseline
loci48[!loci48 %in% loci43]

## Marker in baseline, not in mixture samples
loci43[!loci43 %in% loci48]  # Checked Andy's allele frequency plots, and this marker appears fixed in almost all baseline pops

## Investigate Andy's FreqFisPlots
## Ots_IL-1RA and Ots_SERPC1-209 both have rampant HWE issues (- and + Fis respectively)
## Templin et al. 2011 kept both markers in (Bonferroni corrected had minor issues with IL-1RA)
## Wed Aug 10, 2016, Andy keeps both, Chris and Tyler decided to keep both, so that Andy and I have the "same baseline"

loci42
LocusControl$Publishedlocusnames[loci42]

# Fill in published names with locus names where absent
loci42.published <- apply(cbind(LocusControl$Publishedlocusnames[loci42], loci42), 1, function(row) {
  ifelse(is.na(row[1]), row[2], row[1])
})

loci48.published <- apply(cbind(LocusControl$Publishedlocusnames[loci48], loci48), 1, function(row) {
  ifelse(is.na(row[1]), row[2], row[1])
})

## Which markers were in Templin et al. 2011, but not in loci42?
templin.loci[!templin.loci %in% loci42.published]  # Two with crosses were dropped for LD, S71 is the same as S-71

templin.loci[!templin.loci %in% loci48.published]  # We ran Ots_ZNF330-181 (aka Ots_NRP), but not using in baseline?

# WTF happened to Ots_ZNF330-181 (aka Ots_NRP)???
FreqPop.GCL(sillyvec = PooledNames211, loci = "Ots_NRP")

# Does Ots_NRP appear in scores for all of Andy's baseline pops?
Ots_NRP.PooledNames211 <- sapply(PooledNames211, function(silly) {
  my.gcl <- get(paste(silly, ".gcl", sep = ''))
  "Ots_NRP" %in% dimnames(my.gcl$scores)[[2]]
})
table(Ots_NRP.PooledNames211)  # Missing in 1 population!!!

PooledNames211[!Ots_NRP.PooledNames211]  # Missing from the Kenai

KKENA04.KKENAI03.KKENM06.gcl$n

# This marker does exist for these collections in OceanAK, dunno why Andy dropped it?

# Check to make sure that it was in fact genotyped for all of these collections...

hist(FreqPop.GCL(sillyvec = PooledNames211[Ots_NRP.PooledNames211] , loci = "Ots_NRP"), col = 8, breaks = seq(0, 800, 20), xlab = "Allele Count", main = "Count of Ots_NRP per population")

# Thu Aug 11 10:50:08 2016, confirmed with Andy that Ots_NRP was NOT run for his mixtures
# This is because Ots_NRP was not run for lots of the CI baseline (mostly fixed)
# Given that there are lots of holes for Ots_NRP in the CI baseline, Andy dropped Ots_NRP from his baseline
# This represents a departure from the markers used in Templin et al. 2011



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Create New Groupvec ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Make new groupvec10
pops211groups10 <- readClipboard()  # LowerCI BAseline Summary_KS.xlsx, "Pooling" tab, column Q "Groups10"
groups10 <- unique(pops211groups10)
groupvec10 <- as.numeric(factor(x = pops211groups10, levels = groups10))
writeClipboard(as.character(groupvec10))
colors10 <- c("purple", "green", "orange", "red", "pink", "brown", "blue", "cyan", "violet", "goldenrod")

dput(x = groups10, file = "Objects/groups10.txt")
dput(x = groupvec10, file = "Objects/groupvec10.txt")
dput(x = colors10, file = "Objects/colors10.txt")


## Pop sample sizes
pop.n <- sapply(PooledNames211, function(silly) {get(paste(silly, ".gcl", sep = ''))$n})
sum(pop.n)  # 28,243 fish in baseline
summary(pop.n)  # smallest pop is 42 fish, largest is 396
hist(pop.n, col = 8, breaks = seq(0, 400, 10), main = "Histogram of Populations Size", xlab = "Population Size")


## Group sample sizes
group.n <- setNames(object = sapply(seq(groups10), function(group) {sum(sapply(PooledNames211[groupvec10 == group], function(silly) {get(paste(silly, ".gcl", sep = ''))$n})) } ), nm = groups10)
sapply(PooledNames211[groupvec10 == which(groups10 == "Chignik")], function(silly) {get(paste(silly, ".gcl", sep = ''))$n})  # Only 75 fish in Chignik RG, if we add KCHIG12, we get 141 total for this pop, presuming it pools

## Remove individuals not genotyped for at least 1 locus
original.n <- sapply(PooledNames211, function(silly) {get(paste(silly, ".gcl", sep = ''))$n})
RemoveIndMissLoci.GCL(sillyvec = PooledNames211, proportion = 1/42); beep(2)
postmiss.n <- sapply(PooledNames211, function(silly) {get(paste(silly, ".gcl", sep = ''))$n})

PooledNames211[original.n - postmiss.n > 0]
str(KANDR02.KANDR03.gcl)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Likelihood Profiles 42 Loci PooledNames211 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Allele Frequencies
KMA211Pops_42loci_AlleleCounts <- FreqPop.GCL(sillyvec = PooledNames211, loci = loci42)
str(KMA211Pops_42loci_AlleleCounts)

dir.create("Likelihood Profiles")

KMA211Pops_42loci_Likelihood_Profile <- LeaveOneOutDist.GCL(sillyvec = PooledNames211, loci = loci42, groupvec = groupvec10)
dput(x = KMA211Pops_42loci_Likelihood_Profile, file = "Likelihood Profiles/KMA211Pops_42loci_Likelihood_Profile.txt")

KMA211Pops_42loci_Likelihood_Profile <- dget(file = "Likelihood Profiles/KMA211Pops_42loci_Likelihood_Profile.txt")
str(KMA211Pops_42loci_Likelihood_Profile[[1]], max.level = 1)

# Individual baseline genetic likelihood for RGs assigned (i.e. what is the genotype likelihood of all baseline individuals to RG X)
invisible(lapply(KMA211Pops_42loci_Likelihood_Profile[[1]], function(lst) {boxplot(lst, pch=16, cex = 0.5, notch=TRUE, col=colors10[sort(groupvec10)], ylab="Probability", xlim = c(0, length(PooledNames211) * 1.17), bty = "n", axes = FALSE, xlab = "Population", cex.lab = 1.5)
  axis(side = 2, cex.axis = 1.5)
  axis(side = 1, cex.axis = 1.5, at = c(seq(from = 0, to = length(PooledNames211), by = 50), length(PooledNames211)))
  
  # segments(x0 = c(0, cumsum(table(groupvec10))[-max(groupvec10)]) + 0.5, y0 = 0, x1 = c(0, cumsum(table(groupvec10))[-max(groupvec10)]) + 0.5, y1 = 1, col = colors10, lwd = 4)
  
  legend(x = length(PooledNames211), y = 1, legend = groups10, fill = colors10, bty = "n", cex = 1.3)} ))  # Regional Flat Prior
# Confusion Matrices
KMA211Pops_42loci_Confusion <- ConfusionMatrices.GCL(LeaveOneOutDist = KMA211Pops_42loci_Likelihood_Profile, groupnames = groups10, groupvec = groupvec10, sillyvec = PooledNames211)  # Regional Flat Prior
dput(x = KMA211Pops_42loci_Confusion, file = "Objects/KMA211Pops_42loci_Confusion.txt")
KMA211Pops_42loci_Confusion <- dget(file = "Objects/KMA211Pops_42loci_Confusion.txt")

diag(KMA211Pops_42loci_Confusion[[1]])
# dimnames(KMA211Pops_42loci_Confusion[[1]]) <- list(groups, PCGroups15)

require(lattice)
require(devEMF)

new.colors <- colorRampPalette(c("white", "black"))

# emf(file = "Likelihood Profiles/KMA211Pops_42loci_Confusion.emf", width = 6.5, height = 6.5, family = "Times")
levelplot(KMA211Pops_42loci_Confusion[[1]], col.regions = new.colors, xlab = "Known Origin", ylab = "Mean Genotype Likelihood", 
          at = seq(0, 1, length.out = 100), scales = list(x = list(rot = 90)))
# dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# New LeaveOneOutLikeProfile.GCL.R
KMA211Pops_42loci_Likelihood_Profile_NEW <- LeaveOneOutLikeProfile.GCL(popvec = PooledNames211, loci = loci42, groupvec = groupvec10, groupnames = groups10, groupcomps = NULL, ncores = 6)
# dput(x = KMA211Pops_42loci_Likelihood_Profile_NEW, file = "Likelihood Profiles/KMA211Pops_42loci_Likelihood_Profile_NEW.txt")
KMA211Pops_42loci_Likelihood_Profile_NEW <- dget(file = "Likelihood Profiles/KMA211Pops_42loci_Likelihood_Profile_NEW.txt")
str(KMA211Pops_42loci_Likelihood_Profile_NEW)

PlotLikeProfile.GCL(likeprof = KMA211Pops_42loci_Likelihood_Profile_NEW, popvec = PooledNames211, loci = loci42, groupvec = groupvec10, groupnames = groups10, dir = "Likelihood Profiles", filename = "KMA211")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Read in Additional SILLY's to test for Pooling ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Compare SILLYs between Andy's Lower CI (Templin 2011 + CI) and the 319
andy.sillys <- readClipboard()
andy.sillys <- unlist(lapply(andy.sillys, function(silly) {strsplit(x = silly, split = "\\.")}))

nick.sillys <- readClipboard()
nick.sillys <- unlist(lapply(nick.sillys, function(silly) {strsplit(x = silly, split = "\\.")}))

# Which sillys are in the 319, but not in the Andy's Lower CI baseline
nick.sillys[!nick.sillys %in% andy.sillys]


new.sillys.to.pool <- c("KCHIG12", "KAYAK07", "KKARL07", "KKARL12", "KMONA09", "KPILL13", "KBIGCK08", "KPLEN14", "KLAND12", "KSAPSUK12", "KSAPSUK13", "KBLACH07")
LOKI2R.GCL(sillyvec = new.sillys.to.pool, username = "krshedd", password = password); rm (password)
length(new.sillys.to.pool)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Data QC/Massage ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

require(xlsx); require(beepr)

new.sillys.to.pool

new.sillys.to.pool_SampleSizes <- matrix(data = NA, nrow = length(new.sillys.to.pool), ncol = 4, 
                                        dimnames = list(new.sillys.to.pool, c("Genotyped", "Missing", "Duplicate", "Final")))

#### Check loci
## Get sample size by locus
Original_new.sillys.to.pool_SampleSizebyLocus <- SampSizeByLocus.GCL(new.sillys.to.pool, loci42)
min(Original_new.sillys.to.pool_SampleSizebyLocus)  # 30
sort(apply(Original_new.sillys.to.pool_SampleSizebyLocus,1,min)/apply(Original_new.sillys.to.pool_SampleSizebyLocus,1,max))  # Two under 0.8
table(apply(Original_new.sillys.to.pool_SampleSizebyLocus,1,min)/apply(Original_new.sillys.to.pool_SampleSizebyLocus,1,max) < 0.8)  # 2 SILLY's with at least one locus fail
str(Original_new.sillys.to.pool_SampleSizebyLocus)

## Percent by locus
Original_new.sillys.to.pool_PercentbyLocus <- apply(Original_new.sillys.to.pool_SampleSizebyLocus, 1, function(row) {row / max(row)})
which(apply(Original_new.sillys.to.pool_PercentbyLocus, 2, min) < 0.8)
# writeClipboard(as.character(apply(Original_new.sillys.to.pool_PercentbyLocus, 2, min) < 0.8))

require(lattice)
new.colors <- colorRampPalette(c("black", "white"))
levelplot(t(Original_new.sillys.to.pool_PercentbyLocus), col.regions = new.colors, xlab = "SILLY", ylab = "Locus", at = seq(0, 1, length.out = 100), scales = list(x = list(rot = 90)), aspect = "fill")  # aspect = "iso" will make squares
## This looks pretty good, no holes are evident after tweaking ReadLOKIv2.GCL to ReadLOKIv3.GCL with "onlymarkersuitefish" switch turned to TRUE


#### Check individuals
### Initial
## Get number of individuals per silly before removing missing loci individuals
Original_new.sillys.to.pool_ColSize <- sapply(paste(new.sillys.to.pool, ".gcl", sep = ''), function(x) get(x)$n)
new.sillys.to.pool_SampleSizes[, "Genotyped"] <- Original_new.sillys.to.pool_ColSize


### Missing
## Remove individuals with >20% missing data
new.sillys.to.pool_MissLoci <- RemoveIndMissLoci.GCL(sillyvec = new.sillys.to.pool, proportion = 0.8); beep(8)

## Get number of individuals per silly after removing missing loci individuals
ColSize_new.sillys.to.pool_PostMissLoci <- sapply(paste(new.sillys.to.pool, ".gcl", sep = ''), function(x) get(x)$n)
new.sillys.to.pool_SampleSizes[, "Missing"] <- Original_new.sillys.to.pool_ColSize-ColSize_new.sillys.to.pool_PostMissLoci


### Duplicate
## Check within collections for duplicate individuals.
new.sillys.to.pool_DuplicateCheck95MinProportion <- CheckDupWithinSilly.GCL(sillyvec = new.sillys.to.pool, loci = loci42, quantile = NULL, minproportion = 0.95); beep(8)
new.sillys.to.pool_DuplicateCheckReportSummary <- sapply(new.sillys.to.pool, function(x) new.sillys.to.pool_DuplicateCheck95MinProportion[[x]]$report)

## Remove duplicate individuals
# Don't include "SKANA07"  "SKANA10"  "SKANAL13" (inbred/bottlenecked)
new.sillys.to.pool_RemovedDups <- RemoveDups.GCL(new.sillys.to.pool_DuplicateCheck95MinProportion)

## Get number of individuals per silly after removing duplicate individuals
ColSize_new.sillys.to.pool_PostDuplicate <- sapply(paste(new.sillys.to.pool, ".gcl", sep = ''), function(x) get(x)$n)
new.sillys.to.pool_SampleSizes[, "Duplicate"] <- ColSize_new.sillys.to.pool_PostMissLoci-ColSize_new.sillys.to.pool_PostDuplicate


### Final
new.sillys.to.pool_SampleSizes[, "Final"] <- ColSize_new.sillys.to.pool_PostDuplicate
new.sillys.to.pool_SampleSizes

dir.create("Output")
write.xlsx(new.sillys.to.pool_SampleSizes, file = "Output/new.sillys.to.pool_SampleSizes.xlsx")

## Save post-QC .gcl's as back-up:
dir.create(path = "Raw genotypes/PostQCCollections")
invisible(sapply(new.sillys.to.pool, function(silly) {dput(x = get(paste(silly, ".gcl", sep = '')), file = paste("Raw genotypes/PostQCCollections/" , silly, ".txt", sep = ''))} )); beep(8)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### New Pooling ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

invisible(sapply(new.sillys.to.pool, function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/PostQCCollections/", silly, ".txt", sep = "")), pos = 1)} )); beep(2)


length(c(new.sillys.to.pool, PooledNames211))
KMA223Collections_42loci_AlleleCounts <- FreqPop.GCL(sillyvec = c(new.sillys.to.pool, PooledNames211), loci = loci42)

dir.create("Pooling Fishers Test Results")

# Chignik
grep(pattern = "KCHIG", x = objects(pattern = "\\.gcl"), value = TRUE)
ChignikPoolTest <- list(c("KCHIG95.KCHIG06", "KCHIG12"))
ChignikPoolTestResult <- FishersTest.GCL(freq = KMA223Collections_42loci_AlleleCounts, loci = loci42, tests = ChignikPoolTest)
ChignikPoolTestResult$OverallResults  # 0.9976
dput(x = ChignikPoolTestResult, file = "Pooling Fishers Test Results/ChignikPoolTestResult.txt")
hist(sort(ChignikPoolTestResult$ResultsByLocus$KCHIG95.KCHIG06.KCHIG12$pval), breaks = seq(0, 1, 0.05), col = 8, main = "Histogram of p-values by locus", xlab = "P-values by locus")

PoolCollections.GCL(collections = unlist(ChignikPoolTest), loci = loci42, IDs = NULL, newname = paste(unlist(ChignikPoolTest), collapse = "."))
str(KCHIG95.KCHIG06.KCHIG12.gcl)

# Ayakulik
grep(pattern = "KAYA", x = objects(pattern = "\\.gcl"), value = TRUE)
AyakulikPoolTest <- list(c("KIAYA93.KAYAK06", "KAYAK07"))
AyakulikPoolTestResult <- FishersTest.GCL(freq = KMA223Collections_42loci_AlleleCounts, loci = loci42, tests = AyakulikPoolTest)
AyakulikPoolTestResult$OverallResults  # 0.861
dput(x = AyakulikPoolTestResult, file = "Pooling Fishers Test Results/AyakulikPoolTestResult.txt")
hist(sort(AyakulikPoolTestResult$ResultsByLocus$KIAYA93.KAYAK06.KAYAK07$pval), breaks = seq(0, 1, 0.05), col = 8, main = "Histogram of p-values by locus", xlab = "P-values by locus")

PoolCollections.GCL(collections = unlist(AyakulikPoolTest), loci = loci42, IDs = NULL, newname = paste(unlist(AyakulikPoolTest), collapse = "."))
str(KIAYA93.KAYAK06.KAYAK07.gcl)

KIAYA93.KAYAK06.gcl$n; KAYAK07.gcl$n
str(KAYAK07.gcl)
KAYAK07.gcl$counts[, 1, 1]

# Karluk
grep(pattern = "KKARL", x = objects(pattern = "\\.gcl"), value = TRUE)
KarlukPoolTest <- list(c("KIKAR93.KKARL06", "KKARL07", "KKARL12"))
KarlukPoolTestResult <- FishersTest.GCL(freq = KMA223Collections_42loci_AlleleCounts, loci = loci42, tests = KarlukPoolTest)
KarlukPoolTestResult$OverallResults  # 0.6341
dput(x = KarlukPoolTestResult, file = "Pooling Fishers Test Results/KarlukPoolTestResult.txt")
hist(sort(KarlukPoolTestResult$ResultsByLocus$KIKAR93.KKARL06.KKARL07.KKARL12$pval), breaks = seq(0, 1, 0.05), col = 8, main = "Histogram of p-values by locus", xlab = "P-values by locus")

PoolCollections.GCL(collections = unlist(KarlukPoolTest), loci = loci42, IDs = NULL, newname = paste(unlist(KarlukPoolTest), collapse = "."))
str(KIKAR93.KKARL06.KKARL07.KKARL12.gcl)
KIKAR93.KKARL06.gcl$n; KKARL07.gcl$n; KKARL12.gcl$n


# Monashka + Pillar Creek
KMONA09.gcl$n; KPILL13.gcl$n

PillarPoolTest <- list(c("KMONA09", "KPILL13"))
PillarPoolTestResult <- FishersTest.GCL(freq = KMA223Collections_42loci_AlleleCounts, loci = loci42, tests = PillarPoolTest)
PillarPoolTestResult$OverallResults  # 2e-04
dput(x = PillarPoolTestResult, file = "Pooling Fishers Test Results/PillarPoolTestResult.txt")
hist(sort(PillarPoolTestResult$ResultsByLocus$KMONA09.KPILL13$pval), breaks = seq(0, 1, 0.05), col = 8, main = "Histogram of p-values by locus", xlab = "P-values by locus")
sapply(unlist(PillarPoolTest), function(silly) {get(paste(silly, ".gcl", sep = ''))$n})

PoolCollections.GCL(collections = unlist(PillarPoolTest), loci = loci42, IDs = NULL, newname = paste(unlist(PillarPoolTest), collapse = "."))
str(KMONA09.KPILL13.gcl)

PillarPoolTest2 <- list(c("KPILL13", "KIKAR93.KKARL06"))
PillarPoolTest2Result <- FishersTest.GCL(freq = KMA223Collections_42loci_AlleleCounts, loci = loci42, tests = PillarPoolTest2)
PillarPoolTest2Result$OverallResults  # 0
dput(x = PillarPoolTest2Result, file = "Pooling Fishers Test Results/PillarPoolTest2Result.txt")

PillarPoolTest3 <- list(c("KMONA09", "KIKAR93.KKARL06"))
PillarPoolTest3Result <- FishersTest.GCL(freq = KMA223Collections_42loci_AlleleCounts, loci = loci42, tests = PillarPoolTest3)
PillarPoolTest3Result$OverallResults  # 0.0136
dput(x = PillarPoolTest3Result, file = "Pooling Fishers Test Results/PillarPoolTest3Result.txt")

# Big Creek
grep(pattern = "KBIGCK", x = objects(pattern = "\\.gcl"), value = TRUE)
BigCreekPoolTest <- list(c("KBIGCK04", "KBIGCK08"))
BigCreekPoolTestResult <- FishersTest.GCL(freq = KMA223Collections_42loci_AlleleCounts, loci = loci42, tests = BigCreekPoolTest)
BigCreekPoolTestResult$OverallResults  # 0.9611
dput(x = BigCreekPoolTestResult, file = "Pooling Fishers Test Results/BigCreekPoolTestResult.txt")
hist(sort(BigCreekPoolTestResult$ResultsByLocus$KBIGCK04.KBIGCK08$pval), breaks = seq(0, 1, 0.05), col = 8, main = "Histogram of p-values by locus", xlab = "P-values by locus")

PoolCollections.GCL(collections = unlist(BigCreekPoolTest), loci = loci42, IDs = NULL, newname = paste(unlist(BigCreekPoolTest), collapse = "."))
str(KBIGCK04.KBIGCK08.gcl)
KBIGCK04.gcl$n; KBIGCK08.gcl$n

# Meshik
MeshikPoolTest <- list(c("KMESH06", "KPLEN14", "KLAND12"))
MeshikPoolTestResult <- FishersTest.GCL(freq = KMA223Collections_42loci_AlleleCounts, loci = loci42, tests = MeshikPoolTest)
MeshikPoolTestResult$OverallResults  # 0.516
dput(x = MeshikPoolTestResult, file = "Pooling Fishers Test Results/MeshikPoolTestResult.txt")
hist(sort(MeshikPoolTestResult$ResultsByLocus$KMESH06.KPLEN14.KLAND12$pval), breaks = seq(0, 1, 0.05), col = 8, main = "Histogram of p-values by locus", xlab = "P-values by locus")

PoolCollections.GCL(collections = unlist(MeshikPoolTest), loci = loci42, IDs = NULL, newname = paste(unlist(MeshikPoolTest), collapse = "."))
str(KMESH06.KPLEN14.KLAND12.gcl)
KMESH06.gcl$n; KPLEN14.gcl$n; KLAND12.gcl$n

# Black Hills
BlackHillsPoolTest <- list(c("KBLACH06", "KBLACH07"))
BlackHillsPoolTestResult <- FishersTest.GCL(freq = KMA223Collections_42loci_AlleleCounts, loci = loci42, tests = BlackHillsPoolTest)
BlackHillsPoolTestResult$OverallResults  # 0.9895
dput(x = BlackHillsPoolTestResult, file = "Pooling Fishers Test Results/BlackHillsPoolTestResult.txt")
hist(sort(BlackHillsPoolTestResult$ResultsByLocus$KBLACH06.KBLACH07$pval), breaks = seq(0, 1, 0.05), col = 8, main = "Histogram of p-values by locus", xlab = "P-values by locus")

PoolCollections.GCL(collections = unlist(BlackHillsPoolTest), loci = loci42, IDs = NULL, newname = paste(unlist(BlackHillsPoolTest), collapse = "."))
str(KBLACH06.KBLACH07.gcl)
KBLACH06.gcl$n; KBLACH07.gcl$n

# Nelson/Sapsuk
NelsonPoolTest <- list(c("KNELS06", "KSAPSUK12", "KSAPSUK13"))
NelsonPoolTestResult <- FishersTest.GCL(freq = KMA223Collections_42loci_AlleleCounts, loci = loci42, tests = NelsonPoolTest)
NelsonPoolTestResult$OverallResults  # 0.0804
dput(x = NelsonPoolTestResult, file = "Pooling Fishers Test Results/NelsonPoolTestResult.txt")
hist(sort(NelsonPoolTestResult$ResultsByLocus$KNELS06.KSAPSUK12.KSAPSUK13$pval), breaks = seq(0, 1, 0.05), col = 8, main = "Histogram of p-values by locus", xlab = "P-values by locus")
sapply(unlist(NelsonPoolTest), function(silly) {get(paste(silly, ".gcl", sep = ''))$n})
PoolCollections.GCL(collections = unlist(NelsonPoolTest), loci = loci42, IDs = NULL, newname = paste(unlist(NelsonPoolTest), collapse = "."))
str(KNELS06.KSAPSUK12.KSAPSUK13.gcl)


SapsukPoolTest <- list(c("KSAPSUK12", "KSAPSUK13"))
SapsukPoolTestResult <- FishersTest.GCL(freq = KMA223Collections_42loci_AlleleCounts, loci = loci42, tests = SapsukPoolTest)
SapsukPoolTestResult$OverallResults  # 0.9958
dput(x = SapsukPoolTestResult, file = "Pooling Fishers Test Results/SapsukPoolTestResult.txt")
hist(sort(SapsukPoolTestResult$ResultsByLocus$KSAPSUK12.KSAPSUK13$pval), breaks = seq(0, 1, 0.05), col = 8, main = "Histogram of p-values by locus", xlab = "P-values by locus")
sapply(unlist(SapsukPoolTest), function(silly) {get(paste(silly, ".gcl", sep = ''))$n})


PoolCollections.GCL(collections = unlist(SapsukPoolTest), loci = loci42, IDs = NULL, newname = paste(unlist(SapsukPoolTest), collapse = "."))
str(KSAPSUK12.KSAPSUK13.gcl)

NelsonPoolTest2 <- list(c("KNELS06", "KSAPSUK12.KSAPSUK13"))
NelsonPoolTest2_42loci_AlleleCounts <- FreqPop.GCL(sillyvec = unlist(NelsonPoolTest2), loci = loci42)
NelsonPoolTest2Result <- FishersTest.GCL(freq = NelsonPoolTest2_42loci_AlleleCounts, loci = loci42, tests = NelsonPoolTest2)
NelsonPoolTest2Result$OverallResults  # 0.0902
dput(x = NelsonPoolTest2Result, file = "Pooling Fishers Test Results/NelsonPoolTest2Result.txt")
hist(sort(NelsonPoolTest2Result$ResultsByLocus$KNELS06.KSAPSUK12.KSAPSUK13$pval), breaks = seq(0, 1, 0.05), col = 8, main = "Histogram of p-values by locus", xlab = "P-values by locus")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### HWE for New Pooling ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FreqFis to verify pooling
pooling.collections <- c("KCHIG95.KCHIG06", "KCHIG12", "KIAYA93.KAYAK06", "KAYAK07", "KIKAR93.KKARL06", "KKARL07", "KKARL12", "KMONA09", "KPILL13", "KBIGCK04", "KBIGCK08", "KMESH06", "KPLEN14", "KLAND12", "KBLACH06", "KBLACH07", "KNELS06", "KSAPSUK12", "KSAPSUK13")
length(pooling.collections)

PoolingFreqFis <- FreqFisPlot4SNPs.GCL(sillyvec = pooling.collections, loci = loci42, groupvec = c(1, 1, 2, 2, 3, 3, 3, 4, 4, 5, 5, 6, 6, 6, 7, 7, 8, 8, 8), file = "FreqFisPlots/KMA19PoolingCollections_42loci_FreqFisPlot.pdf")
dput(x = PoolingFreqFis, file = "FreqFisPlots/KMA19PoolingCollections_42loci_FreqFis.txt")

# Histogram of Fis values by silly
sapply(rownames(PoolingFreqFis$Fis), function(silly) {hist(PoolingFreqFis$Fis[silly, ], breaks = seq(-1, 1, 0.05), col = 8, xlab = "Fis", main = silly)})

# Ho Fis Fst Table
PoolingHoFisTable <- HoFisFstTable.GCL(sillyvec = pooling.collections, loci = loci42, dir = "Output")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HWE in pooling collections
dir.create("Genepop")
length(pooling.collections)
mito.loci <- which(LocusControl$ploidy[loci42] == 1)
gcl2Genepop.GCL(sillyvec = pooling.collections, loci = loci42[-mito.loci], path = "Genepop/KMA19PoolingCollections_41nuclearloci.txt", VialNums = TRUE)

HWE.pooling.collections <- ReadGenepopHWE.GCL(file = "Genepop/KMA19PoolingCollections_41nuclearloci.txt.P")
str(HWE.pooling.collections)
HWE.pooling.collections$SummaryPValues["Overall Loci", ]  # KMONA09 and KPILL13 are marginal for HWE, overall Fis is -0.04 & -0.03, respectively (.DIV)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HWE in pooled populations
pooled.pops <- c("KCHIG95.KCHIG06.KCHIG12", "KIAYA93.KAYAK06.KAYAK07", "KIKAR93.KKARL06.KKARL07.KKARL12", "KMONA09.KPILL13", "KBIGCK04.KBIGCK08", "KMESH06.KPLEN14.KLAND12", "KBLACH06.KBLACH07", "KNELS06.KSAPSUK12.KSAPSUK13", "KSAPSUK12.KSAPSUK13")
length(pooled.pops)
gcl2Genepop.GCL(sillyvec = pooled.pops, loci = loci42[-mito.loci], path = "Genepop/KMA9PooledPops_41nuclearloci.txt", VialNums = TRUE)

HWE.pooled.pops <- ReadGenepopHWE.GCL(file = "Genepop/KMA9PooledPops_41nuclearloci.txt.P")
str(HWE.pooled.pops)
HWE.pooled.pops$SummaryPValues["Overall Loci", ]  # KMONA09.KPILL13 is marginal for HWE, overall Fis is -0.03
HWE.pooled.pops$DataByPop
unique(HWE.pooled.pops$DataByPop$Pop)[4]


# Fis for each and HWE P-value
PillarFisHWETab <- cbind(sapply(c("KMONA09", "KPILL13"), function(silly) {PoolingFreqFis$Fis[silly, loci42[-mito.loci]]}),
subset(x = HWE.pooled.pops$DataByPop, subset = Pop == unique(HWE.pooled.pops$DataByPop$Pop)[4], select = "WC Fis"),
subset(x = HWE.pooled.pops$DataByPop, subset = Pop == unique(HWE.pooled.pops$DataByPop$Pop)[4], select = "PValue"))

write.table(x = PillarFisHWETab, file = "Output/PillarFisHWETable.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Final Population List ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Decision is to keep things that pool, but not add any new populations

PooledNames211
pooled.pops

poooling.mgsub <- c(Chignik = ChignikPoolTest,
     Ayakulik = AyakulikPoolTest,
     Karluk = KarlukPoolTest,
     BigCreek = BigCreekPoolTest,
     Meshik = MeshikPoolTest,
     BlackHills = BlackHillsPoolTest,
     Nelson = NelsonPoolTest)

pooling.mgsub.pattern <- sapply(poooling.mgsub, function(pop) {pop[1]})
pooling.mgsub.replacement <- sapply(poooling.mgsub, function(pop) {paste(pop, collapse = ".")})

require(qdap)
KMA211Pops <- mgsub(pattern = pooling.mgsub.pattern, replacement = pooling.mgsub.replacement, text.var = PooledNames211)
dput(x = KMA211Pops, file = "Objects/KMA211Pops.txt")

KMA211Pops[!KMA211Pops %in% PooledNames211]
PooledNames211[!PooledNames211 %in% KMA211Pops]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Data QC ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Check loci
## Get sample size by locus
KMA211Pops_SampleSizebyLocus <- SampSizeByLocus.GCL(KMA211Pops, loci42)
min(KMA211Pops_SampleSizebyLocus)  # 8
sort(apply(KMA211Pops_SampleSizebyLocus,1,min)/apply(KMA211Pops_SampleSizebyLocus,1,max))  # Several under 0.8
table(apply(KMA211Pops_SampleSizebyLocus,1,min)/apply(KMA211Pops_SampleSizebyLocus,1,max) < 0.8)  # 36 SILLY's with at least one locus fail
str(KMA211Pops_SampleSizebyLocus)

## Percent by locus
KMA211Pops_PercentbyLocus <- apply(KMA211Pops_SampleSizebyLocus, 1, function(row) {row / max(row)})
which(apply(KMA211Pops_PercentbyLocus, 2, min) < 0.8)
str(KMA211Pops_PercentbyLocus)
# writeClipboard(as.character(apply(KMA211Pops_PercentbyLocus, 2, min) < 0.8))

require(lattice)
new.colors <- colorRampPalette(c("black", "white"))
levelplot(t(KMA211Pops_PercentbyLocus), col.regions = new.colors, xlab = "SILLY", ylab = "Locus", at = seq(0, 1, length.out = 100), scales = list(x = list(rot = 90)), aspect = "fill")  # aspect = "iso" will make squares
## This looks pretty good, no holes are evident after tweaking ReadLOKIv2.GCL to ReadLOKIv3.GCL with "onlymarkersuitefish" switch turned to TRUE

KMA211Pops_SampleSizebyLocus[8, ]
str(KANDR02.KANDR03.gcl)


## Look for weird scores
table(KANDR02.KANDR03.gcl$scores)
unique(names(table(KANDR02.KANDR03.gcl$scores)))

# What are the alleles present in each silly
KMA211Pops.scores <- sapply(KMA211Pops, function(silly) {
  names(table(get(paste(silly, ".gcl", sep = ''))$scores))
})

# Do any silly's have weird alleles?
table(sapply(KMA211Pops.scores, length))

KMA211Pops.scores[which(sapply(KMA211Pops.scores, length) > 6)]



for(silly in names(which(sapply(KMA211Pops.scores, length) > 6))) {
  my.gcl <- get(paste(silly, ".gcl", sep = ''))
  
  counts=my.gcl$counts
  scores=my.gcl$scores
  n=my.gcl$n
  attributes=my.gcl$attributes
  
  scores <- gsub(pattern = "Unk", replacement = "0", scores)

  assign(paste0(silly,".gcl"),list(counts=counts,scores=scores,n=n,attributes=attributes),pos=1)
}

unique(names(table(KANDR02.KANDR03.gcl$scores)))


## Save final .gcl's as back-up:
dir.create("Raw genotypes/PostQCPooledPopsClean")
invisible(sapply(KMA211Pops, function(silly) {dput(x = get(paste(silly, ".gcl", sep = '')), file = paste("Raw genotypes/PostQCPooledPopsClean/" , silly, ".txt", sep = ''))} )); beep(8)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Clean workspace; dget .gcl objects and Locus Control ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls(all = TRUE))
setwd("V:/Analysis/4_Westward/Chinook/CSRI Westward Commercial Harvest 2014-2016/Baseline")
source("H:/R Source Scripts/Functions.GCL_KS.R")
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Get objects
baselineobjects <- sapply(list.files(path = "Objects"), function(file) {unlist(strsplit(x = file, split = ".txt"))} , USE.NAMES = FALSE)
invisible(sapply(baselineobjects, function(file) {assign(x = file, value = dget(file = paste("Objects/", file, ".txt", sep = "")), pos = 1)} ))
LocusControl <- OriginalLocusControl_loci48; rm(OriginalLocusControl_loci48, OriginalLocusControl_loci52)

## Get pops
require(beepr)
invisible(sapply(KMA211Pops, function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/PostQCPooledPopsClean/", silly, ".txt", sep = "")), pos = 1)} )); beep(2)
objects(pattern = "\\.gcl")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Freq Fis Plots Post Pooling New Collections ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dir.create("FreqFisPlots")
KMA211PopsPostPool_42loci_FreqFis <- FreqFisPlot4SNPs.GCL(sillyvec = KMA211Pops, loci = loci42, groupvec = groupvec10, groupcol = colors10, file = "FreqFisPlots/KMA211PopsPostPool_42loci_FreqFisPlot.pdf", dot.cex = 0.8)
str(KMA211PopsPostPool_42loci_FreqFis)
dput(x = KMA211PopsPostPool_42loci_FreqFis, file = "FreqFisPlots/KMA211PopsPostPool_42loci_FreqFis.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### MDS ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Dump a GENEPOP file for genind
require(adegenet)
gcl2Genepop.GCL(sillyvec = KMA211Pops, loci = loci42, path = "Genepop/KMA211Pops_42loci.gen", VialNums = TRUE)
gcl2Genepop.GCL(sillyvec = KMA211Pops, loci = loci42[-mito.loci], path = "Genepop/KMA211Pops_41nuclearloci.gen", VialNums = TRUE)
genind <- read.genepop(file = "Genepop/KMA211Pops_41nuclearloci.gen")

genpop <- genind2genpop(genind)

AdegenetNei211Pop41nuclearloci <- dist.genpop(genpop, method = 1, diag = TRUE, upper = TRUE)
dir.create("Trees")
dput(x = AdegenetNei211Pop41nuclearloci, file = "Trees/AdegenetNei211Pop41nuclearloci.txt")
str(AdegenetNei211Pop41nuclearloci)

require(ape)
Nei211NJtree <- nj(AdegenetNei211Pop41nuclearloci)
str(Nei211NJtree)

par(mar = c(5, 5, 5, 5), oma = c(3, 0, 0, 0))
plot.phylo(x = Nei211NJtree, cex = 0.5, no.margin = TRUE, type = "p")
axisPhylo(1, las = 1, backward = FALSE)

Nei211NJtree$tip.label <- popnames211
par(mar = c(5, 5, 5, 5), oma = c(3, 0, 0, 0))
plot.phylo(x = Nei211NJtree, cex = 0.3, no.margin = TRUE, type = "p")
axisPhylo(1, las = 1, backward = FALSE)

write.tree(Nei211NJtree, "Trees/NJofNei211Popstree.nex", tree.names = popnames211)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## MDS
library('rgl')

MDS <- cmdscale(as.matrix(AdegenetNei211Pop41nuclearloci), k = 3)  # Did in base R as it kept crashing in RStudio...dunno why.
dput(x = MDS, file = "Objects/MDSAdegenetNei211Pop41nuclearloci.txt")
MDS <- dget(file = "Objects/MDSAdegenetNei211Pop41nuclearloci.txt")

x <- as.vector(MDS[, 1])   
y <- as.vector(MDS[, 2])
z <- as.vector(MDS[, 3])

par(family = "serif")
plot3d(x, y, z + abs(range(z)[1]), xlab = '', ylab = '', zlab = '', aspect = FALSE, col = colors10[groupvec10], size = 0.5, type = 's', axes = TRUE, box = TRUE, top = TRUE, cex = 1)
box3d()
# plot3d(x, y, z + abs(range(z)[1]), aspect = FALSE, col = "black", size = 0.5, type = 'h', box = TRUE, axes = FALSE, top = FALSE, add = TRUE)  # adds pins to spheres
texts3d(x, y, z + abs(range(z)[1]), adj = c(-0.2, 0), text = seq(KMA211Pops), font = 2, cex = 0.8, add = TRUE, top = TRUE, axes = FALSE)  # adds numbers to points(adj moves the numbers around the points)

dir.create("MDS")
rgl.snapshot("MDS/MDSAdegenetNei211Pop41nuclearloci.png", fmt="png", top=TRUE )


# Plot individually
par(mar = c(4, 4, 4, 4))
plot(x, col = colors10[groupvec10], pch = 16)
plot(y, col = colors10[groupvec10], pch = 16)
popnames211[which.min(y)]; popnames211[which.max(y)]
plot(z, col = colors10[groupvec10], pch = 16)
popnames211[which.min(z)]; popnames211[which.max(z)]

cbind(groups10, colors10)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Likelihood Profiles 42 Loci Post Pooling New Collections ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
KMA211PopsPostPool_42loci_Likelihood_Profile <- LeaveOneOutDist.GCL(sillyvec = KMA211Pops, loci = loci42, groupvec = groupvec10)
dput(x = KMA211PopsPostPool_42loci_Likelihood_Profile, file = "Likelihood Profiles/KMA211PopsPostPool_42loci_Likelihood_Profile.txt")

KMA211PopsPostPool_42loci_Likelihood_Profile <- dget(file = "Likelihood Profiles/KMA211PopsPostPool_42loci_Likelihood_Profile.txt")
str(KMA211PopsPostPool_42loci_Likelihood_Profile[[1]], max.level = 1)

# Individual baseline genetic likelihood for RGs assigned (i.e. what is the genotype likelihood of all baseline individuals to RG X)
invisible(lapply(KMA211PopsPostPool_42loci_Likelihood_Profile[[1]], function(lst) {boxplot(lst, pch=16, cex = 0.5, notch=TRUE, col=colors10[sort(groupvec10)], ylab="Probability", xlim = c(0, length(KMA211Pops) * 1.17), bty = "n", axes = FALSE, xlab = "Population", cex.lab = 1.5)
  axis(side = 2, cex.axis = 1.5)
  axis(side = 1, cex.axis = 1.5, at = c(seq(from = 0, to = length(KMA211Pops), by = 50), length(KMA211Pops)))
  
  # segments(x0 = c(0, cumsum(table(groupvec10))[-max(groupvec10)]) + 0.5, y0 = 0, x1 = c(0, cumsum(table(groupvec10))[-max(groupvec10)]) + 0.5, y1 = 1, col = colors10, lwd = 4)
  
  legend(x = length(KMA211Pops), y = 1, legend = groups10, fill = colors10, bty = "n", cex = 1.3)} ))  # Regional Flat Prior
# Confusion Matrices
KMA211PopsPostPool_42loci_Confusion <- ConfusionMatrices.GCL(LeaveOneOutDist = KMA211PopsPostPool_42loci_Likelihood_Profile, groupnames = groups10, groupvec = groupvec10, sillyvec = KMA211Pops)  # Regional Flat Prior
dput(x = KMA211PopsPostPool_42loci_Confusion, file = "Objects/KMA211PopsPostPool_42loci_Confusion.txt")
KMA211PopsPostPool_42loci_Confusion <- dget(file = "Objects/KMA211PopsPostPool_42loci_Confusion.txt")

diag(KMA211PopsPostPool_42loci_Confusion[[1]])
# dimnames(KMA211PopsPostPool_42loci_Confusion[[1]]) <- list(groups, PCGroups15)

require(lattice)
require(devEMF)

new.colors <- colorRampPalette(c("white", "black"))

# emf(file = "Likelihood Profiles/KMA211PopsPostPool_42loci_Confusion.emf", width = 6.5, height = 6.5, family = "Times")
levelplot(KMA211PopsPostPool_42loci_Confusion[[1]], col.regions = new.colors, xlab = "Known Origin", ylab = "Mean Genotype Likelihood", 
          at = seq(0, 1, length.out = 100), scales = list(x = list(rot = 90)))
# dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# New LeaveOneOutLikeProfile.GCL.R
KMA211PopsPostPool_42loci_Likelihood_Profile_NEW <- LeaveOneOutLikeProfile.GCL(popvec = KMA211Pops, loci = loci42, groupvec = groupvec10, groupnames = groups10, groupcomps = NULL, ncores = 6)
# dput(x = KMA211PopsPostPool_42loci_Likelihood_Profile_NEW, file = "Likelihood Profiles/KMA211PopsPostPool_42loci_Likelihood_Profile_NEW.txt")
KMA211PopsPostPool_42loci_Likelihood_Profile_NEW <- dget(file = "Likelihood Profiles/KMA211PopsPostPool_42loci_Likelihood_Profile_NEW.txt")
str(KMA211PopsPostPool_42loci_Likelihood_Profile_NEW)

PlotLikeProfile.GCL(likeprof = KMA211PopsPostPool_42loci_Likelihood_Profile_NEW, popvec = KMA211Pops, loci = loci42, groupvec = groupvec10, groupnames = groups10, dir = "Likelihood Profiles", filename = "KMA211")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Pairwise Fst Tree ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
