#install.packages(c("ape","phytools","geiger","fields","animation","phangorn","rgdal","sp"))
source("R/source_functions.R")

#install ImageMagick (http://www.imagemagick.org)
#for macs 10.9 and up, make sure to have X11 installed. you can get it here http://xquartz.macosforge.org/landing/

library(ape)
library(phytools)
library(geiger)
library(fields)
library(animation)
library(phangorn)
library(rgdal)
library(sp)
#library(ppgm)

#Set bounds for analysis
#bounds <- list(sigsq = c(min = 0, max = 1000000), SE = c(0, 0.1), alpha = c(min = 0, max = 150), a = c(min = -1, max = 1), slope = c(min = -100, max = 100), lambda = c(min = 0, max = 1), kappa = c(min = 0, max = 1), delta = c(min = 0, max = 10), drift = c(min = -100, max = 100))

load("data/occurrences.RData")
occu <- occurrences
load("data/beastLeache.RData")
which_run <- sample(1:length(beastLeache), 100)
ex_mytree <- beastLeache[which_run]
#This run is with NO FOSSILS
#Runs may take a while, so do a few trial runs before you commit to many permutations.

#extra code to push back divergence dates to MRCA = 70MYA
for(x in 1:length(ex_mytree)) {
  ex_mytree[[x]]$edge.length <- ex_mytree[[x]]$edge.length * 70/max(ex_mytree[[x]]$edge.length)
}

trialestBM <- ppgm(occurrences = occu, trees = ex_mytree, model = "BM", which.biovars = c(1, 4, 6, 15))
trialestOU <- ppgm(occurrences = occu, trees = ex_mytree, model = "OU", which.biovars = c(1, 4, 6, 15))
trialestEB <- ppgm(occurrences = occu, trees = ex_mytree, model = "EB", which.biovars = c(1, 4, 6, 15))

#Make table for results
for(trees in 1:length(trialestBM$model_min)){
    for(traits in 1:length(trialestBM$model_min[[1]][[1]])){
      temp_min <- cbind(trialestBM$model_min[[trees]][[1]][[traits]]$fitted$aicc, trialestOU$model_min[[trees]][[1]][[traits]]$fitted$aicc, trialestEB$model_min[[trees]][[1]][[traits]]$fitted$aicc)
      temp_max <- cbind(trialestBM$model_max[[trees]][[1]][[traits]]$fitted$aicc, trialestOU$model_max[[trees]][[1]][[traits]]$fitted$aicc, trialestEB$model_max[[trees]][[1]][[traits]]$fitted$aicc)
      rownames(temp_min) <- rownames(temp_max) <- paste("aicc", trees, "_", traits, sep = "")
      if(traits == 1 & trees == 1){
        print_table_min <- temp_min
        print_table_max <- temp_max
      } else {
        print_table_min <- rbind(print_table_min, temp_min)
        print_table_max <- rbind(print_table_max, temp_max)
      }
    }
}
colnames(print_table_min) <- colnames(print_table_max) <- c("BM", "OU", "EB")

print_table_min <- t(print_table_min)
print_table_max <- t(print_table_max)

#combine all phylo and permuts to get average aicc by trait
print(cbind(MAT = apply(print_table_min[, seq(1, dim(print_table_min)[2], 3)], 1, min), TS = apply(print_table_min[, seq(2, dim(print_table_min)[2],3)], 1, min), MTCM = apply(print_table_min, 1, min), PS = apply(print_table_min[, seq(3, dim(print_table_min)[2], 3)], 1, min)), digits = 5)
print(cbind(MAT = apply(print_table_max[, seq(1, dim(print_table_max)[2], 3)], 1, min), TS = apply(print_table_max[, seq(2, dim(print_table_max)[2],3)], 1, min), MTCM = apply(print_table_max, 1, min), PS = apply(print_table_max[, seq(3, dim(print_table_max)[2], 3)], 1, min)), digits = 5)

#NOW INCLUDING FOSSILS
#permuts only matter when including fossils!!!
permuts <- 100
load("data/fossils.RData")
fossils$early_age <- ceiling(fossils$early_age)
fossils$late_age <- ceiling(fossils$late_age)
fossilRanges <- (apply(fossils[, 3:4], 1, function(x) {x[1] - x[2]}) + 1)

#get fossil table in proper format
manipulatedFosssils <- array(NA, dim = c(sum(fossilRanges), 4))
z <- 1
for(i in 1:length(fossils[, 1])){
  for(j in 1:fossilRanges[i]){
    manipulatedFosssils[z, ] <- unlist(c(fossils$late_age[i] + j - 1, fossils[i, c(5, 6, 2)]))
    z <- z + 1
  }
}
#Note: if supplying biovarFossils, supply all 19 or will throw an error
biovarFossils <- getBioclimVars(manipulatedFosssils[, 1:3], which.biovars = 1:19)

#find all the generic scelop, remove specific scelop
fossilsedges <- gsub(1, "", manipulatedFosssils[, 4])
fossilsedges[which(fossilsedges == "")] <- NA
fossilsedges[which(fossilsedges != T)] <- "undulatus"

#cut down fossil sample with subsample
ex_fossils <- biovarFossils[is.na(fossilsedges), ]
trialestWFBM <- ppgm(occurrences = occu, fossils = ex_fossils, trees = ex_mytree, model = "BM", permut = permuts, which.biovars = c(1, 4, 6, 15), bounds = list(alpha = c(0, 1)), control = list(niter = 20))
trialestWFOU <- ppgm(occurrences = occu, fossils = ex_fossils, trees = ex_mytree, model = "OU", permut = permuts, which.biovars = c(1, 4, 6, 15), bounds = list(alpha = c(0, 1)), control = list(niter = 20))
trialestWFEB <- ppgm(occurrences = occu, fossils = ex_fossils, trees = ex_mytree, model = "EB", permut = permuts, which.biovars = c(1, 4, 6, 15), bounds = list(alpha = c(0, 1)), control = list(niter = 20))

#make table for results with fossils
models <- c("BM", "OU", "EB")
for(trees in 1:length(trialestWFBM$model_min)){
  for(permut in 1:length(trialestWFBM$model_min[[trees]])){
    for(traits in 1:length(trialestWFBM$model_min[[trees]][[permut]])){
      temp_min <- cbind(trialestWFBM$model_min[[trees]][[permut]][[traits]]$fitted$aicc, trialestWFOU$model_min[[trees]][[permut]][[traits]]$fitted$aicc, trialestWFEB$model_min[[trees]][[permut]][[traits]]$fitted$aicc)
      temp_max <- cbind(trialestWFBM$model_max[[trees]][[permut]][[traits]]$fitted$aicc, trialestWFOU$model_max[[trees]][[permut]][[traits]]$fitted$aicc, trialestWFEB$model_max[[trees]][[permut]][[traits]]$fitted$aicc)
      rownames(temp_min) <- rownames(temp_max) <- paste("aicc", trees, permut, traits, sep = "")
      if(trees == 1 & traits ==1){
        print_table_min <- temp_min
        print_table_max <- temp_max
      }
      else{
        print_table_min <- rbind(print_table_min, temp_min)
        print_table_max <- rbind(print_table_max, temp_max)
      }
    }
  }
}
colnames(print_table_min) <- colnames(print_table_max) <- models

print_table_min <- t(print_table_min)
print_table_max <- t(print_table_max)

#combine all phylo and permuts to get average aicc by trait
print(cbind(MAT = apply(print_table_min[, seq(1, dim(print_table_min)[2], 3)], 1, min), TS = apply(print_table_min[, seq(2, dim(print_table_min)[2], 3)], 1, min), MTCM = apply(print_table_min, 1, min), PS = apply(print_table_min[, seq(3, dim(print_table_min)[2], 3)], 1, min)), digits = 5)
print(cbind(MAT = apply(print_table_max[, seq(1, dim(print_table_max)[2], 3)], 1, min), TS = apply(print_table_max[, seq(2, dim(print_table_max)[2], 3)], 1, min), MTCM = apply(print_table_max, 1, min), PS = apply(print_table_max[, seq(3, dim(print_table_max)[2], 3)], 1, min)), digits = 5)