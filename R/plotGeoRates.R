#' @title plotGeoRates
#' @usage plotGeoRates(geo_center, geo_size, time_int, trees, path="")
#' @param geo_center change in geographic center of suitable climate envelope, see
#' @param geo_size change in geographic size of suitable climate envelope
#' @param time_int time intervals to plot
#' @param trees distribution of phylogenies
#' @param path path to the directory where the results to be saved
#' @details Creates plot with gray background of all pairwise comparisons of change in geo center and area through time. Blue points on top show the sequential change in geo center and expansion/contraction for all lineages
#' @return plots of geo rate
#' @seealso \code{getGeoRates}
#' @export
#' @author A. Michelle Lawing, Alexandra F. C. Howard
#' @examples


plotGeoRates <- function(geo_center, geo_size, time_int, trees, path=""){
  jpeg(paste(path,"geo_rates.jpg",sep=""),width=960,height=960,quality=100,pointsize=20)
  par(mfrow=c(2,1),mar=c(4,4,1,1),mgp=c(2.25,1,0))
  if(is.finite(max(geo_center))){ymax=max(geo_center)}
  if(!is.finite(max(geo_center))){ymax=4000}
  plot(1,1,type="l",xlim=c(0,20),ylim=c(0,ymax),xlab="",ylab="Change in Geo Center")
  for(t in 1:length(trees)){
    for(l in 1:length(geo_center[t,,1,1])){
      points(time_int[,],geo_center[t,l,,],col="lightgray")
    }
  }
  for(t in 1:length(trees)){
    for(l in 1:length(geo_center[1,,1,1])){
      for(i in 1:20){
        points(time_int[1,(i+1)],geo_center[t,l,i,(i+1)],col="blue")
      }
    }
  }
  plot(1,1,type="l",xlim=c(0,20),ylim=c(min(geo_size),max(geo_size)),xlab="Time",ylab="Change in Geo Area")
  for(t in 1:length(trees)){
    for(l in 1:length(geo_size[t,,1,1])){
      points(time_int[,],geo_size[t,l,,],col="lightgray")
    }
  }
  for(t in 1:length(trees)){
    for(l in 1:length(geo_size[1,,1,1])){
      for(i in 1:20){
        points(time_int[1,(i+1)],geo_size[t,l,i,(i+1)],col="blue")
      }
    }
  }
  #reset par
  par(mfrow=c(1,1),mar=c(5,4,4,2),mgp=c(3,1,0))
  dev.off()
}

#' @title plotGeoRatesCon
#' @usage plotGeoRatesCon(geo_center, geo_size, time_int, trees, path="")
#' @param geo_center change in geographic center of suitable climate envelope, see
#' @param geo_size change in geographic size of suitable climate envelope
#' @param time_int time intervals to plot
#' @param trees distribution of phylogenies
#' @param path path to the directory where the results to be saved
#' @details Creates plot with gray background of all pairwise comparisons of change in geo center and area through time. Blue points on top show the sequential change in geo center and expansion/contraction for all lineages
#' @return plots of geo rate
#' @seealso \code{getGeoRates}
#' @export
#' @author A. Michelle Lawing, Alexandra F. C. Howard
#' @examples
#' data(beastLeache)
#' data(occurrences)
#' tree <- beastLeache[[25]]
#' sp_data_min<- tapply(occurrences[,4],occurrences$Species,min)
#' sp_data_max<- tapply(occurrences[,4],occurrences$Species,max)
#' treedata_min <- treedata(tree,sp_data_min,sort=TRUE,warnings=F)
#' treedata_max <- treedata(tree,sp_data_max,sort=TRUE,warnings=F)
#' full_est <- nodeEstimateFossils(treedata_min,treedata_max)
#' node_est <- full_est$est
#' example_getEnvelopes <- getEnvelopes(treedata_min, treedata_max, node_est)
#' example_getGeoRate <- getGeoRate(example_getEnvelopes, tree,which.biovars=1)
#' plotGeoRatesCon(example_getGeoRate$geo_center,example_getGeoRate$geo_size,example_getGeoRate$time_int, trees = trees[[1]])

plotGeoRatesCon <- function(geo_center, geo_size, time_int, trees, path=""){
  jpeg(paste(path,"geo_rates.jpg",sep=""),width=960,height=960,quality=100,pointsize=20)
  par(mfrow=c(2,1),mar=c(4,4,1,1),mgp=c(2.25,1,0))
  if(is.finite(max(geo_center))){ymax=max(geo_center)}
  if(!is.finite(max(geo_center))){ymax=4000}
  plot(1,1,type="l",xlim=c(0,20),ylim=c(0,ymax),xlab="",ylab="Change in Geo Center")
    for(l in 1:length(geo_center[,1,1])){
      points(time_int[,],geo_center[l,,],col="lightgray")
    }
    for(l in 1:length(geo_center[,1,1])){
      for(i in 1:20){
        points(time_int[1,(i+1)],geo_center[l,i,(i+1)],col="blue")
      }
    }
  plot(1,1,type="l",xlim=c(0,20),ylim=c(min(geo_size),max(geo_size)),xlab="Time",ylab="Change in Geo Area")
    for(l in 1:length(geo_size[,1,1])){
      points(time_int[,],geo_size[l,,],col="lightgray")
    }
    for(l in 1:length(geo_size[,1,1])){
      for(i in 1:20){
        points(time_int[1,(i+1)],geo_size[l,i,(i+1)],col="blue")
      }
    }
  #reset par
  par(mfrow=c(1,1),mar=c(5,4,4,2),mgp=c(3,1,0))
  dev.off()
}
