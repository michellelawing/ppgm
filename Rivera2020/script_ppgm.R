#install.packages(c("ape","phytools","geiger","fields","animation","phangorn","rgdal","sp"))
#install.packages("https://github.com/michellelawing/ppgm")
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
bounds <- list(sigsq = c(min = 0, max = 1000000), SE = c(0, 0.1), alpha = c(min = 0, max = 150), a = c(min = -1, max = 1), slope = c(min = -100, max = 100), lambda = c(min = 0, max = 1), kappa = c(min = 0, max = 1), delta = c(min = 0, max = 10), drift = c(min = -100, max = 100))

data(occurrences)
data(beastLeache)
which_run <- sample(1:length(beastLeache), 100)
ex_mytree <- beastLeache[which_run]

#This run is with NO FOSSILS
#Runs may take a while, so do a few trial runs before you commit to many permutations.
#trialest <- ppgm(occurrences = occu, fossils = F, trees = ex_mytree, fossils.edges = FALSE, model = "estimate", permut = 1, which.biovars = c(1, 4, 15), path = "scratch/p_", plot.TraitGram = T, plot.AnimatedMaps = F, plot.GeoRates = F)
#trialest1 <- ppgm(occurrences = occu, fossils = F, trees = ex_mytree, fossils.edges = FALSE, model = "estimate", permut = 1, which.biovars = c(6), path = "scratch/p_", plot.TraitGram = T, plot.AnimatedMaps = F, plot.GeoRates = F)
#for one trait and one model, use code below to get the aicc mean for all 100 trees
#mean(unlist(sapply(1:length(trialest1$model_min), function(trees) trialest1$model_min[[trees]][[1]][[1]]$fitted['aicc'])))
#skip run, save time, and load data
data(NoFoss100treesB1B4B15)
data(NoFossilsResults100treesB6)

#Make table for results
models <- c("BM", "OU", "EB", "trend", "lambda", "kappa", "delta", "drift", "white")
for(trees in 1:length(trialest$model_min)){
    for(traits in 1:length(trialest$model_min[[1]][[1]])){
      temp_min <- cbind(unlist(sapply(models, function(z) trialest$model_min[[trees]][[1]][[traits]]$fitted[[z]]['aicc'])))
      temp_max <- cbind(unlist(sapply(models, function(z) trialest$model_max[[trees]][[1]][[traits]]$fitted[[z]]['aicc'])))
      colnames(temp_min) <- colnames(temp_max) <- paste("aicc", trees, 1, traits, sep = "")
      if(trees == 1){
        print_table_min <- temp_min
        print_table_max <- temp_max
      }
      else {
        print_table_min <- cbind(print_table_min, temp_min)
        print_table_max <- cbind(print_table_max, temp_max)
      }
      rownames(print_table_min) <- rownames(print_table_max) <- models
    }
}
print_table_min <- print_table_min[-4,]
print_table_max <- print_table_max[-4,]

for(trees in 1:length(trialest1$model_min)){
      temp_min <- cbind(unlist(sapply(models, function(z) trialest1$model_min[[trees]][[1]][[1]]$fitted[[z]]['aicc'])))
      temp_max <- cbind(unlist(sapply(models, function(z) trialest1$model_max[[trees]][[1]][[1]]$fitted[[z]]['aicc'])))
      colnames(temp_min) <- colnames(temp_max) <- paste("aicc", trees, 1, 1, sep = "")
      if(trees == 1){
        print_table1_min <- temp_min
        print_table1_max <- temp_max
      }
      else {
        print_table1_min <- cbind(print_table1_min, temp_min)
        print_table1_max <- cbind(print_table1_max, temp_max)
      }
      rownames(print_table_min) <- rownames(print_table_max) <- models[-4]
}
#combine all phylo and permuts to get average aicc by trait
print(cbind(MAT = apply(print_table_min[, seq(1, dim(print_table_min)[2], 3)], 1, min), TS = apply(print_table_min[, seq(2, dim(print_table_min)[2],3)], 1, min), MTCM = apply(print_table1_min, 1, min), PS = apply(print_table_min[, seq(3, dim(print_table_min)[2], 3)], 1, min)), digits = 5)
print(cbind(MAT = apply(print_table_max[, seq(1, dim(print_table_max)[2], 3)], 1, min), TS = apply(print_table_max[, seq(2, dim(print_table_max)[2],3)], 1, min), MTCM = apply(print_table1_max, 1, min), PS = apply(print_table_max[, seq(3, dim(print_table_max)[2], 3)], 1, min)), digits = 5)

#NOW INCLUDING FOSSILS
#permuts only matter when including fossils!!!
#permuts <- 100
data(fossils)
ex_fossils <- fossils
fossils$early_age <- ceiling(fossils$early_age)
fossils$late_age <- ceiling(fossils$late_age)
fossilRanges <- (apply(fossils[, 3:4], 1, function(x) {x[1] - x[2]}) + 1)

#get fossil table in proper format
manipulatedFosssils <- array(NA, dim = c(sum(fossilRanges), 4))
z <- 1
for(i in 1:length(ex_fossils[, 1])){
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
#trialestWF <- ppgm(occurrences = occurrences, fossils = ex_fossils,trees = ex_mytree, fossils.edges = F, model = "estimate", permut = permuts, which.biovars = c(1, 4, 15), path = "scratch/q_", plot.TraitGram = TRUE, plot.AnimatedMaps = FALSE,plot.GeoRates = FALSE, bounds = list(alpha = c(0, 1)), control = list(niter = 20))
#trialestWF2 <- ppgm(occurrences = occurrences, fossils = ex_fossils,trees = ex_mytree, fossils.edges = F, model = "estimate", permut = permuts, which.biovars = c(6), path = "scratch/q_", plot.TraitGram = TRUE, plot.AnimatedMaps = FALSE, plot.GeoRates = FALSE, bounds = list(alpha = c(0, 1)), control = list(niter = 20))
#skip run, save time, and load data
data(Fossils100treesB1B4B15)
data(Fossils100treesB6)

#make table for results with fossils
models <- c("BM", "OU", "EB", "trend", "lambda", "kappa", "delta", "drift", "white")
for(trees in 1:length(trialestWF$model_min)){
  for(permut in 1:length(trialestWF$model_min[[trees]])){
    for(traits in 1:length(trialestWF$model_min[[trees]][[permut]])){
      temp_min <- cbind(unlist(sapply(models,function(z) trialestWF$model_min[[trees]][[permut]][[traits]]$fitted[[z]]['aicc'])))
      temp_max <- cbind(unlist(sapply(models,function(z) trialestWF$model_max[[trees]][[permut]][[traits]]$fitted[[z]]['aicc'])))
      colnames(temp_min) <- colnames(temp_max) <- paste("aicc", trees, permut, traits, sep = "")
      if(trees == 1){
        print_table_min <- temp_min
        print_table_max <- temp_max
      }
      else{
        print_table_min <- cbind(print_table_min, temp_min)
        print_table_max <- cbind(print_table_max, temp_max)
      }
      rownames(print_table_min) <- rownames(print_table_max) <- models
    }
  }
}

for(trees in 1:length(trialestWF2$model_min)){
  for(permut in 1:length(trialestWF2$model_min[[trees]])){
  temp_min <- cbind(unlist(sapply(models, function(z) trialestWF2$model_min[[trees]][[permut]][[1]]$fitted[[z]]['aicc'])))
  temp_max <- cbind(unlist(sapply(models, function(z) trialestWF2$model_max[[trees]][[permut]][[1]]$fitted[[z]]['aicc'])))
  colnames(temp_min) <- paste("aicc", trees, permut, 1, sep = "")
  colnames(temp_max) <- paste("aicc", trees, permut, 1, sep = "")
  if(trees == 1){
    print_table1_min <- temp_min
    print_table1_max <- temp_max
  } else {
    print_table1_min <- cbind(print_table1_min, temp_min)
    print_table1_max <- cbind(print_table1_max, temp_max)
  }
  rownames(print_table1_min) <- rownames(print_table1_max) <- models
  }
}
#combine all phylo and permuts to get average aicc by trait
print(cbind(MAT = apply(print_table_min[, seq(1, dim(print_table_min)[2], 3)], 1, min), TS = apply(print_table_min[, seq(2, dim(print_table_min)[2], 3)], 1, min), MTCM = apply(print_table1_min, 1, min), PS = apply(print_table_min[, seq(3, dim(print_table_min)[2], 3)], 1, min)), digits = 5)
print(cbind(MAT = apply(print_table_max[, seq(1, dim(print_table_max)[2], 3)], 1, min), TS = apply(print_table_max[, seq(2, dim(print_table_max)[2], 3)], 1, min), MTCM = apply(print_table1_max, 1, min), PS = apply(print_table_max[, seq(3, dim(print_table_max)[2], 3)], 1, min)), digits = 5)

#below are the average parameter estimates for the "best model" for the runs
parameters <- array(NA, dim = c(4, 4))
colnames(parameters) <- c("B1", "B4", "B6", "B15")
rownames(parameters) <- c("min", "max", "minWF", "maxWF")
parameters[1, c(1, 2, 4)] <- rowMeans(sapply(1:length(ex_mytree), function(trees) sapply(1:3, function(traits) trialest$model_min[[trees]][[1]][[traits]]$fitted$OU$alpha)), na.rm = T)
parameters[1, 3] <- mean(sapply(1:length(ex_mytree), function(trees) trialest1$model_min[[trees]][[1]][[1]]$fitted$OU$alpha), na.rm = T)
parameters[2, c(1, 2, 4)] <- rowMeans(sapply(1:length(ex_mytree), function(trees) sapply(1:3, function(traits) trialest$model_max[[trees]][[1]][[traits]]$fitted$OU$alpha)), na.rm = T)
parameters[2, 3] <- mean(sapply(1:length(ex_mytree), function(trees) trialest1$model_max[[trees]][[1]][[1]]$fitted$OU$alpha), na.rm = T)
parameters[3, c(1, 2, 4)] <- rowMeans(sapply(1:length(ex_mytree), function(trees) sapply(1:3, function(traits) trialestWF$model_min[[trees]][[1]][[traits]]$fitted$OU$alpha)), na.rm = T)
parameters[3, 3] <- mean(sapply(1:length(ex_mytree), function(trees) trialestWF2$model_min[[trees]][[1]][[1]]$fitted$OU$alpha), na.rm = T)
parameters[4, c(1, 2, 4)] <- rowMeans(sapply(1:length(ex_mytree), function(trees) sapply(1:3, function(traits) trialest$model_max[[trees]][[1]][[traits]]$fitted$lambda$lambda)))
parameters[4, 3] <- mean(sapply(1:length(ex_mytree), function(trees) trialestWF2$model_max[[trees]][[1]][[1]]$fitted$delta$delta), na.rm = T)
print(parameters)

#plot ppgmMESS
par(mar = c(2, 1, 1, 1))

#get envelopes 1,4,15 then 6 ... put appropriate order
cem_min <- cbind(trialest$cem[, 1], trialest$cem[, 2], trialest1$cem[, 1],trialest$cem[, 3])
cem_max <- cbind(trialest$cem[, 7], trialest$cem[, 8], trialest1$cem[, 3],trialest$cem[, 9])
rownames(cem_min) <- rownames(cem_max) <- rownames(trialest$cem)

#extract one variable from the results of multivar run
relist1 <- lapply(lapply(1:length(ex_mytree), function(x) array(unlist(trialest$node_est), dim = c(2, 52, 3, 100))[,,1,x]), list)
relist2 <- lapply(lapply(1:length(ex_mytree), function(x) array(unlist(trialest$node_est), dim = c(2, 52, 3, 100))[,,2,x]), list)
relist3 <- lapply(lapply(1:length(ex_mytree), function(x) array(unlist(trialest1$node_est), dim = c(2, 52, 1, 100))[,,1,x]), list)
relist4 <- lapply(lapply(1:length(ex_mytree), function(x) array(unlist(trialest$node_est), dim = c(2, 52, 3, 100))[,,3,x]), list)

#MESS, will produce plots in the working directory
mess <- ppgmMESS(cem_min, cem_max, est = list(relist1, relist2, relist3, relist4), tree = ex_mytree, fossils = ex_fossils, timeslice = c(2, 5, 13, 20), which.biovars = c(1, 4, 6, 15), which.plot = "none")

#PHYSIOLOGICAL model
load("data/paleoclimate.Rdata")
timeslice <- c(2, 5, 13, 20)
Tb <- c(28, 32, 35, 38)
par(mar = c(2, 1, 1, 1))
colorscheme <- c("lightgray", "red", "#990000", "#330000")

for(tb in 1:length(Tb)){
  for(p in 1:length(timeslice)){
    hr <- 6.12 + (0.74 * (paleoclimate[[(timeslice[p])]][, 11] / 10 - Tb[tb]))
    hr[hr < 4] <- 1
    hr[hr >= 4 & hr < 7] <- 2
    hr[hr >= 7 & hr < 10] <- 3
    hr[hr >= 10] <- 4
  
    spdata <- SpatialPoints(paleoclimate[[(timeslice[p] + 1)]][, 2:3])
    proj4string(spdata)  <- CRS("+init=epsg:4326")
    spdata <- spTransform(spdata, CRS("+init=epsg:26978"))
    if(sum(ex_fossils[, 1] == (timeslice[p] + 1)) != 0){
      spfossils <- SpatialPoints(ex_fossils[, 2:3])
      proj4string(spfossils)  <- CRS("+init=epsg:4326")
      spfossils <- spTransform(spfossils, CRS("+init=epsg:26978"))
      spfossils <- spfossils[ex_fossils[, 1] == (timeslice[p] + 1), ]
    }
    #pdf(paste("Physio", Tb[tb], "Time", timeslice[p], ".pdf", sep = ""), width = 80, height = 80, pointsize = 100, useDingbats = F)
    plot(spdata, cex = 1, xlab = "", ylab = "", axes = FALSE, pch = 16, col = colorscheme[hr])
    if(sum(ex_fossils[,1] == (timeslice[p] + 1)) != 0){
      points(spfossils, cex = 2, pch = 16, col = "black") #change cex to 4 to get large points if saving to pdf
    }
    #dev.off()
  }
}