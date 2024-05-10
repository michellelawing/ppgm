#' @title plotAnimatedPPGMMultiPhylo
#' @description This function creates an animated gif showing the change in modeled suitable habitat through time in geographic space. It requires ImageMagick or GraphicsMagick to be previously installed in the operating system. This is easy to do if you have macports. Just type sudo port install ImageMagick into terminal.
#' @usage plotAnimatedPPGMMultiPhylo(envelope, tree, filename="ppgm.gif", 
#' which.biovars, path="", use.paleoclimate=TRUE, paleoclimateUser=NULL, layerAge=c(0:20))
#' @param envelope the min and max envelope of each lineage for each time slice
#' @param tree the phylogeny or multiple phylogenies that show the relationship between species
#' @param filename filename of output
#' @param which.biovars A vector of the numbers of the bioclimate variables that should be returned. The bioclimate variables number correspond to the Hijmans table at (https://www.worldclim.org/data/bioclim.html).
#' @param path path to the directory where the results should be saved
#' @param use.paleoclimate if left blank, default North America paleoclimate data is used. If FALSE, user submitted paleoclimate must be provided
#' @param paleoclimateUser list of data frames with paleoclimates, must be dataframes with columns: GlobalID, Longitude, Latitude, bio1, bio2,...,bio19. (see \code{getBioclimvars()}).
#' @param layerAge vector with the ages of the paleoclimate dataframes, if using user submitted paleoclimate data
#' @details Animated gif of species richness through time
#' @author A. Michelle Lawing, Alexandra F. C. Howard
#' @importFrom utils data
#' @importFrom animation saveGIF
#' @importFrom graphics points
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics text
#' @export
#' @examples
#' data(sampletrees)
#' data(occurrences)
#' sampletrees <- sample(sampletrees,5)
#' biooccu <- getBioclimVars(occurrences, which.biovars=1)
#' sp_data_min<- tapply(biooccu[,4],biooccu$Species,min)
#' sp_data_max<- tapply(biooccu[,4],biooccu$Species,max)
#' treedata_min <- treedata_max <- node_est <- envelope <- list()
#' \dontrun{for (tr in 1:length(sampletrees)){
#'   treedata_min[[tr]] <- geiger::treedata(sampletrees[[tr]],sp_data_min,sort=TRUE,warnings=F)
#'   treedata_max[[tr]] <- geiger::treedata(sampletrees[[tr]],sp_data_max,sort=TRUE,warnings=F)
#'   full_est <- nodeEstimateEnvelopes(treedata_min[[tr]],treedata_max[[tr]])
#'   node_est[[tr]] <- full_est$est
#'   envelope[[tr]] <- getEnvelopes(treedata_min[[tr]], treedata_max[[tr]], node_est[[tr]])
#' }
#' animatedplot <- plotAnimatedPPGMMultiPhylo(envelope,sampletrees,which.biovars=1, path=tempdir())}


plotAnimatedPPGMMultiPhylo <- function(envelope, tree, filename="ppgm.gif", which.biovars, path="", use.paleoclimate=TRUE, paleoclimateUser=NULL, layerAge=c(0:20)){
  if(path==""){out=getwd()}
  if(path!=""){out=paste(getwd(),"/",substr(path,1,nchar(path)-1),sep="")}
  #load paleoclimate data
  if(use.paleoclimate) {
    paleoclimate <- paleoclimate #uses paleoclimate data from package
  } else {
    if(is.null(paleoclimateUser)) {
      stop("paleoclimateUser argument must be provided when use.paleoclimate is FALSE.") #uses user inputted paleoclimate
    } else {
      paleoclimate <- paleoclimateUser
    }
  }
  temp_min<-as.list(array(NA,dim=length(tree)))
  temp_max<-as.list(array(NA,dim=length(tree)))
  richnesscount<-as.list(array(NA,dim=length(tree)))
  for(tr in 1:length(tree)){
    temp_min[[tr]]<-lapply(1:length(paleoclimate),function(i){
      temp<-lapply(1:length(which.biovars),function(j){getTimeSlice(layerAge[[i]],tree[[tr]],envelope[[tr]][,2,j])})
      temp<-t(array(unlist(temp),dim=c(length(unlist(temp[[1]]$edge)),2*length(which.biovars))))
      return(temp)})
    temp_max[[tr]]<-lapply(1:length(paleoclimate),function(i){
      temp<-lapply(1:length(which.biovars),function(j){getTimeSlice(layerAge[[i]],tree[[tr]],envelope[[tr]][,5,j])})
      temp<-t(array(unlist(temp),dim=c(length(unlist(temp[[1]]$edge)),2*length(which.biovars))))
      return(temp)})
    for (j in 1:length(paleoclimate)){
      hld<-array(0,dim=length(paleoclimate[[j]][,1]))
      for(i in 1:length(temp_min[[tr]][[j]][1,])){
        matching <- sapply(1:length(which.biovars),function(x){
          paleoclimate[[j]][,which.biovars[x]+3]>temp_min[[tr]][[j]][1:length(which.biovars)*2,i][x] & 
            paleoclimate[[j]][,which.biovars[x]+3]<temp_max[[tr]][[j]][1:length(which.biovars)*2,i][x]
          })
        matching<-which(rowSums(matching)==length(which.biovars),arr.ind=TRUE)
        hld[matching]<-hld[matching]+1
      }
      hld[which(hld==0,arr.ind=TRUE)]=NA
    }
    richnesscount[[tr]] <-hld
  }
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  animation::saveGIF(for(i in 1:length(paleoclimate)){
    richnesscountMean<-rowMeans(array(unlist(lapply(1:length(tree),function(tr) richnesscount[[tr]][[i]])),dim=c(length(richnesscount[[1]][[i]]),length(tree))),na.rm=TRUE)
    lq<-array(unlist(lapply(1:length(tree),function(tr) richnesscount[[tr]][[i]])),dim=c(length(richnesscount[[1]][[i]]),length(tree)))
    richnesscountMIN<-suppressWarnings(unlist(lapply(1:length(lq[,1]),function(r) min(lq[r,],na.rm=TRUE))))
    richnesscountMAX<-suppressWarnings(unlist(lapply(1:length(lq[,1]),function(r) max(lq[r,],na.rm=TRUE))))
    par(mar=c(3,0,3,1))
    plot(paleoclimate[[i]][,2:3],cex=0.5,xlab="",ylab="",axes=FALSE,pch=16,col="lightgray",xlim=c(-200,0),ylim=c(0,90))
    graphics::points(paleoclimate[[i]][,2:3],cex=0.5,pch=16,col=grDevices::colorRampPalette(c("#FFE5CC", "#FF8000", "#990000"))(length(temp_min[[tr]][[1]][1,]))[richnesscountMean],xlim=c(-200,0),ylim=c(0,90))
    center<-colMeans(paleoclimate[[i]][,2:3])
    scale<-(cbind(paleoclimate[[i]][,2]-center[1],paleoclimate[[i]][,3]-center[2]))/3
    translate<-cbind(scale[,1]-160,scale[,2]+25)
    graphics::points(translate,cex=0.5,pch=16,col="lightgray")
    graphics::points(translate,cex=0.5,pch=16,col=colorRampPalette(c("#FFE5CC", "#FF8000", "#990000"))(length(temp_min[[tr]][[1]][1,]))[richnesscountMIN])
    translate<-cbind(scale[,1]-35,scale[,2]+25)
    graphics::points(translate,cex=0.5,pch=16,col="lightgray")
    graphics::points(translate,cex=0.5,pch=16,col=colorRampPalette(c("#FFE5CC", "#FF8000", "#990000"))(length(temp_min[[tr]][[1]][1,]))[richnesscountMAX])
    graphics::text(-160,40,"Min")
    graphics::text(-35,40,"Max")
    graphics::text(-90,85,"Mean")
    graphics::text(-90,90,"Modeled Richness over a distribution of trees")
  },movie.name=filename,outdir=out)
}