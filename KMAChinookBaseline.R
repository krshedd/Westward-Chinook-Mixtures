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


# Create Directories
dirs <- c("Objects", "Likelihood Profiles", "Raw gentoypes", "BAYES")
sapply(dirs, dir.create)
dir.create(path = "BAYES/Baseline")

# Move BAYES baseline .bse file
file.copy(from = "V:/Analysis/2_Central/Chinook/Lower Cook Inlet/2015/Baseline/BAYES/Baseline/LCI211pops43loci.bse", 
          to = "BAYES/Baseline/LCI211pops43loci.bse")

# Save objects
dput(x = loci43, file = "Objects/loci43.txt")
dput(x = PooledNames211, file = "Objects/PooledNames211.txt")
# dput(x = PooledNames211ordered, file = "Objects/PooledNames211ordered.txt")  # Don't know why this is here, as Andy did NOT use it
dput(x = popnames211, file = "Objects/popnames211.txt")
dput(x = groups14, file = "Objects/groups14.txt")
dput(x = groupvec14, file = "Objects/groupvec14.txt")
dput(x = colors14, file = "Objects/colors14.txt")

# Confirm how baseline was made
# basefortran=CreateBaseline.GCL(sillyvec=PooledNames211,loci=loci43,dir="BAYES/Baseline",basename="LCI211pops43loci",type="BAYES",groupvec=NULL)
# CreateLocusControl.GCL(markersuite = "Chinook2011_52SNPs", username = username, password = password)
# loci43 <- loci52[-c(7,11,13,21,26,29,33,49,52)]


# Make new groupvec10
loci43pops211groups10 <- readClipboard()
groups10 <- unique(pops211groups10)
groupvec10 <- as.numeric(factor(x = pops211groups10, levels = groups10))
writeClipboard(as.character(groupvec10))

dput(x = groups10, file = "Objects/groups10.txt")
dput(x = groupvec10, file = "Objects/groupvec10.txt")


# Read in data
username <- "krshedd"
password <- "********"

CreateLocusControl.GCL(markersuite = "Chinook_Kodiak_2016_48SNPs", username = username, password = password)
loci48 <- LocusControl$locusnames

# Extra markers that mixture samples were genotyped for not in baseline
loci48[!loci48 %in% loci43]

# Marker in baseline, not in mixture samples
loci43[!loci43 %in% loci48]  # This marker appears fixed in almost all baseline pops

# Need to create loci 42 and modify .bse file
