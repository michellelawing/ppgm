#' @title plotAnimatedPPGM
#' @description This function creates an animated gif showing the change in modelled suitable habitat through time in geographic space.
#' @usage plotAnimatedPPGM(envelope, tree, filename="ppgm.gif", which.biovars, 
#' path="", use.paleoclimate=TRUE, paleoclimateUser=NULL, layerAge=c(0:20))
#' @param envelope the min and max envelope of each lineage for each time slice
#' @param tree the phylogeny or multiple phylogenies that show the relationship between species
#' @param filename desired filename of output
#' @param which.biovars A vector of the numbers of the bioclimate variables that should be returned. The bioclimate variables number correspond to the Hijmans table at (https://www.worldclim.org/data/bioclim.html).
#' @param path path to the directory where the results should be saved
#' @param use.paleoclimate if left blank, default North America paleoclimate data is used. If FALSE, user submitted paleoclimate must be provided
#' @param paleoclimateUser list of data frames with paleoclimates, must be dataframes with columns: GlobalID, Longitude, Latitude, bio1, bio2,...,bio19. (see \code{getBioclimvars()}).
#' @param layerAge vector with the ages of the paleoclimate dataframes, if using user submitted paleoclimate data
#' @details Requires ImageMagick or GraphicsMagick to be installed on the operating system. This is easy to do if you have macports. Just type sudo port install ImageMagick into terminal.
#' @return An animated gif of species through time
#' @author A. Michelle Lawing, Alexandra F. C. Howard, Maria-Aleja Hurtado-Materon
#' @importFrom utils data
#' @importFrom animation saveGIF
#' @importFrom graphics points
#' @importFrom grDevices colorRampPalette
#' @export
#' @examples
#' data(sampletrees)
#' data(occurrences)
#' tree <- sampletrees[[25]]
#' biooccu <- getBioclimVars(occurrences, which.biovars=1)
#' sp_data_min<- tapply(biooccu[,4],biooccu$Species,min)
#' sp_data_max<- tapply(biooccu[,4],biooccu$Species,max)
#' treedata_min <- geiger::treedata(tree,sp_data_min,sort=TRUE,warnings=F)
#' treedata_max <- geiger::treedata(tree,sp_data_max,sort=TRUE,warnings=F)
#' \dontrun{full_est <- nodeEstimateEnvelopes(treedata_min,treedata_max)
#' node_est <- full_est$est
#' example_getEnvelopes <- getEnvelopes(treedata_min, treedata_max, node_est)
#' animatedplot <- plotAnimatedPPGM(example_getEnvelopes,tree,which.biovars=1,path=tempdir())}


plotAnimatedPPGM<-function(envelope, tree, filename="ppgm.gif", which.biovars, path="", use.paleoclimate=TRUE, paleoclimateUser=NULL, layerAge=c(0:20)){
  #load paleoclimate data: isotopically scaled paleoclimate bioclimate variables for North America
  if(use.paleoclimate) {
    paleoclimate <- paleoclimate #uses paleoclimate data from package
  } else {
    if(is.null(paleoclimateUser)) {
      stop("paleoclimateUser argument must be provided when use.paleoclimate is FALSE.") #uses user inputted paleoclimate
    } else {
      paleoclimate <- paleoclimateUser
    }
  }
  temp_min<-lapply(1:length(paleoclimate),function(i){
    temp<-lapply(1:length(which.biovars),function(j){getTimeSlice(layerAge[[i]],tree,envelope[,2,j])})
    temp<-t(array(unlist(temp),dim=c(length(unlist(temp[[1]]$edge)),2*length(which.biovars))))
    return(temp)})
  temp_max<-lapply(1:length(paleoclimate),function(i){
    temp<-lapply(1:length(which.biovars),function(j){getTimeSlice(layerAge[[i]],tree,envelope[,5,j])})
    temp<-t(array(unlist(temp),dim=c(length(unlist(temp[[1]]$edge)),2*length(which.biovars))))
    return(temp)})
  richnesscount<-lapply(1:length(paleoclimate), function(j){
    hld<-array(0,dim=length(paleoclimate[[j]][,1]))
    for(i in 1:length(temp_min[[j]][1,])){
      matching<-
        sapply(1:length(which.biovars),function(x){paleoclimate[[j]][,which.biovars[x]+3]>temp_min[[j]][1:length(which.biovars)*2,i][x] & paleoclimate[[j]][
          ,which.biovars[x]+3]<temp_max[[j]][1:length(which.biovars)*2,i][x]})
      matching<-which(rowSums(matching)==length(which.biovars),arr.ind=TRUE)
      hld[matching]<-hld[matching]+1
    }
    hld[which(hld==0,arr.ind=TRUE)]=NA
    return(hld)
  })
  animation::saveGIF(for(i in 1:length(paleoclimate)){
    plot(paleoclimate[[i]][,2:3],cex=0.5,pch=16,col="lightgray",xlim=c(-180,0),ylim=c(0,90))
    graphics::points(paleoclimate[[i]][,2:3],cex=0.5,pch=16,col=grDevices::colorRampPalette(c("#FFE5CC", "#FF8000", "#990000"))(length(temp_min[[1]][1,]))[richnesscount[[i]]],xlim=c(-180,0),ylim=c(0,90))
  },movie.name=filename,outdir=getwd())
}
