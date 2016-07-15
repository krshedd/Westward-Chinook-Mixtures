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
# 3) Create new groupvec
# 4) Likelihood profile


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


## Make new groupvec10
pops211groups10 <- readClipboard()
groups10 <- unique(pops211groups10)
groupvec10 <- as.numeric(factor(x = pops211groups10, levels = groups10))
writeClipboard(as.character(groupvec10))
colors10 <- c("purple", "green", "orange", "red", "pink", "brown", "blue", "cyan", "violet", "goldenrod")

dput(x = groups10, file = "Objects/groups10.txt")
dput(x = groupvec10, file = "Objects/groupvec10.txt")
dput(x = colors10, file = "Objects/colors10.txt")


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




## Extra markers that mixture samples were genotyped for not in baseline
loci48[!loci48 %in% loci43]

## Marker in baseline, not in mixture samples
loci43[!loci43 %in% loci48]  # Checked Andy's allele frequency plots, and this marker appears fixed in almost all baseline pops

loci42 <- loci43[loci43 %in% loci48]
dput(loci42, file = "Objects/loci42.txt")

## Need to create loci 42 and modify .bse file


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

## Locus Control
LocusControl <- dget(file = "Objects/OriginalLocusControl_loci48.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PooledNames211 <- dget(file = "Objects/PooledNames211.txt")
PooledNames211

## Get pops
require(beepr)
invisible(sapply(PooledNames211, function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/PostQCPooledPops/", silly, ".txt", sep = "")), pos = 1)} )); beep(2)
objects(pattern = "\\.gcl")


## Get objects
baselineobjects <- sapply(list.files(path = "Objects"), function(file) {unlist(strsplit(x = file, split = ".txt"))} )
invisible(sapply(baselineobjects, function(file) {assign(x = file, value = dget(file = paste("Objects/", file, ".txt", sep = "")), pos = 1)} ))
LocusControl <- OriginalLocusControl_loci48; rm(OriginalLocusControl_loci48, OriginalLocusControl_loci52)





## Pop sample sizes
summary(sapply(PooledNames211, function(silly) {get(paste(silly, ".gcl", sep = ''))$n}))

original.n <- sapply(PooledNames211, function(silly) {get(paste(silly, ".gcl", sep = ''))$n})
RemoveIndMissLoci.GCL(sillyvec = PooledNames211, proportion = 1/42); beep(8)
postmiss.n <- sapply(PooledNames211, function(silly) {get(paste(silly, ".gcl", sep = ''))$n})

PooledNames211[original.n - postmiss.n > 0]


## Look for NA, etc.
KMA211Pops_42loci_AlleleCounts <- FreqPop.GCL(sillyvec = PooledNames211, loci = loci42)
str(KMA211Pops_42loci_AlleleCounts)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Likelihood Profiles 42 Loci ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

emf(file = "Likelihood Profiles/KMA211Pops_42loci_Confusion.emf", width = 6.5, height = 6.5, family = "Times")
levelplot(KMA211Pops_42loci_Confusion[[1]], col.regions = new.colors, xlab = "Known Origin", ylab = "Mean Genotype Likelihood", 
          at = seq(0, 1, length.out = 100), scales = list(x = list(rot = 90)))
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# New LeaveOneOutLikeProfile.GCL.R
source("V:/DATA/R_GEN/JJs GCL/LeaveOneOutlIkeProfile.GCL.R")
# KMA211Pops_42loci_Likelihood_Profile_NEW <- LeaveOneOutLikeProfile.GCL(popvec = KMA473Pops, loci = loci89, groupvec = KMA473PopsGroupVec15, groupnames = Groups15, groupcomps = NULL, ncores = 6)
# dput(x = KMA211Pops_42loci_Likelihood_Profile_NEW, file = "Likelihood Profiles/KMA211Pops_42loci_Likelihood_Profile_NEW.txt")
KMA211Pops_42loci_Likelihood_Profile_NEW <- dget(file = "Likelihood Profiles/KMA211Pops_42loci_Likelihood_Profile_NEW.txt")
str(KMA211Pops_42loci_Likelihood_Profile_NEW)

source("V:/DATA/R_GEN/JJs GCL/PlotLikeProfile.GCL.R")
PlotLikeProfile.GCL(likeprof = KMA211Pops_42loci_Likelihood_Profile_NEW, popvec = KMA473Pops, loci = loci89, groupvec = KMA473PopsGroupVec15, groupnames = Groups15, dir = "Likelihood Profiles")
likeprof = KMA211Pops_42loci_Likelihood_Profile_NEW; popvec = KMA473Pops; loci = loci89; groupvec = KMA473PopsGroupVec15; groupnames = Groups15; dir = "Likelihood Profiles"
col = Colors15

