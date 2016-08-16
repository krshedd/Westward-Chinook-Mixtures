
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sillyvec <- KMA211Pops
loci <- loci42
dir <- "Trees"
nboots <- 1000
ncores <- 16
returnbootstrapFst <- FALSE

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ptm <- proc.time()

fstatdir=paste(dir,"\\", length(sillyvec), "Pops", length(loci), "Loci_","fstatfile.dat",sep="")  

while(!require("ape")){intall.packages("ape")}

while(!require("hierfstat")){intall.packages("hierfstat")}

nsillys=length(sillyvec)

maxsillychar=nchar(nsillys)+1

nloci=length(loci)

ploidy=LocusControl$ploidy[loci]

nalleles=LocusControl$nalleles[loci]

maxchar=nchar(nalleles)+1
names(maxchar)=loci

alleles=LocusControl$alleles[loci]

my.gcl=lapply(sillyvec,function(silly){get(paste(silly,".gcl",sep=""),pos=1)})
names(my.gcl)=sillyvec

n=sapply(sillyvec,function(silly){my.gcl[[silly]]$n})
names(n)=sillyvec

scores=lapply(sillyvec,function(silly){scores0 <- my.gcl[[silly]]$scores[,loci,] ; scores0[scores0 %in% c("Unk", "XXX", "0")] <- NA ; scores0})
names(scores)=sillyvec

counts=lapply(sillyvec,function(silly){
  sapply(1:n[silly],function(i){
    paste(c(match(silly,sillyvec),sapply(loci,function(locus){
      ifelse(is.na(scores[[silly]][i,locus,1]),paste(rep(0,ploidy[locus]*maxchar[locus]),collapse=""),paste(sapply(1:ploidy[locus],function(allele){
        paste(c(rep(0,maxchar[locus]-nchar(match(scores[[silly]][i,locus,allele],alleles[[locus]]))),match(scores[[silly]][i,locus,allele],alleles[[locus]])),collapse="")
      }),collapse=""))
    })),collapse=" ")
  })
})      
names(counts)=sillyvec

fstat=paste(nsillys,nloci,max(nalleles),max(maxchar),sep=" ")

fstat=rbind(fstat,cbind(loci))

fstat=rbind(fstat,cbind(as.vector(unlist(counts))))

write.table(fstat,fstatdir,row.names=FALSE,col.names=FALSE,quote=FALSE) 

dat=read.fstat.data(fstatdir)

# Fst=array(0,c(nsillys,nsillys),dimnames=list(sillyvec,sillyvec))

pairs=combn(sillyvec,2)

pairnames=apply(pairs,2,function(col){paste(col,collapse=".")})

dimnames(pairs)[[2]]=pairnames

# vc=vector("list",choose(nsillys,2))
# names(vc)=pairnames

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

while(!require(foreach)){install.packages("foreach")}

while(!require(doParallel)){install.packages("doParallel")}

while(!require(parallel)){install.packages("parallel")}

if(ncores > detectCores()){stop("'ncores' is greater than the number of cores available on machine")}

cl=makePSOCKcluster(ncores)

registerDoParallel(cl,cores=ncores)  


ptm <- proc.time()
vc <- foreach(pair=pairnames) %dopar% {
  
  sillys=pairs[,pair]
  sillyIND=sapply(sillys,function(silly){match(silly,sillyvec)})
  mydat=dat[dat[,1]==sillyIND[1] | dat[,1]==sillyIND[2],loci]
  myn=n[sillys]
  levels=rep(1:2,myn)
  
  mito.loci <- which(ploidy[loci] == 1)
  
  # Create matrix of variance components with order loci
  vc.pair <- rbind(
    # Diploid loci
    t(sapply(loci[-mito.loci], function(locus) {
      mydata = data.frame(levels,mydat[,locus])
      varcomp(data.matrix(mydata),diploid=TRUE)$overall
    } )),
    
    # Haploid loci
    t(sapply(loci[mito.loci], function(locus) {
      mydata = data.frame(levels,mydat[,locus])
      c(varcomp(data.matrix(mydata),diploid=FALSE)$overall,0)
    } ))
  )[loci, ]
  
  # Replace NAs with 0s
  vc.pair[is.na(vc.pair)] <- 0
  vc.pair
}

proc.time() - ptm

names(vc) <- pairnames


# Calculate Fst from VC for each pair
Fst.ls <- lapply(vc, function(pair) {
  sum(pair[, 1]) / sum(pair)
})


# Create matrix of pairwise Fst
Fst <- array(0,c(nsillys,nsillys),dimnames=list(sillyvec,sillyvec))

Fst[lower.tri(Fst)] <- unlist(Fst.ls)
Fst[upper.tri(Fst)] <- t(Fst)[upper.tri(Fst)]

# Create tree
tree=nj(Fst)

trees.bootstrapFst <- foreach(i=seq(nboots)) %dopar% {
  
  temploci=sample(loci,nloci,replace=TRUE)
  tempFst.ls <- lapply(vc, function(pair) {
    sum(pair[temploci, 1]) / sum(pair[temploci, 1:3])
  })
  tempFst <- array(0,c(nsillys,nsillys),dimnames=list(sillyvec,sillyvec))
  tempFst[lower.tri(tempFst)] <- unlist(tempFst.ls)
  tempFst[upper.tri(tempFst)] <- t(Fst)[upper.tri(tempFst)]
  
  list(trees = nj(tempFst), bootstrapFst = tempFst)
}

# stopCluster(cl)

bootstrapFst <- lapply(trees.bootstrapFst, function(i) {
  i[["bootstrapFst"]]
})

trees <- lapply(trees.bootstrapFst, function(i) {
  i[["trees"]]
})

bootstrap=prop.clades(tree,trees)


if(returnbootstrapFst) {
  PairwiseFstTree <- list(tree=tree,bootstrap=bootstrap,PairwiseFst=Fst,vc=vc,BootstrapFst=bootstrapFst)
} else {
  PairwiseFstTree <- list(tree=tree,bootstrap=bootstrap,PairwiseFst=Fst,vc=vc)
}

proc.time() - ptm

dput(x = PairwiseFstTree, file = paste(dir,"\\", length(sillyvec), "Pops", length(loci), "Loci_","PairwiseFstTree.txt",sep="")  )
# return(PairwiseFstTree)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

