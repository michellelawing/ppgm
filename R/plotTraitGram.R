#' @title plotTraitGram
#' @description Combine the node estimates based on random or specified fossil placement and plot them on a phylotraitgram in a specified directory.
#' @usage plotTraitGram(treedata_min, treedata_max, node_est, fossils=FALSE, 
#' which.biovars, path="", use.paleoclimate=TRUE, paleoclimateUser=NULL,
#' layerAge=c(0:20))
#' @param treedata_min a tree data object with the min estimate of the climate envelope
#' @param treedata_max a tree data object with the max estimate of the climate envelope
#' @param node_est the estimate of all the nodes, both min and max
#' @param fossils a matrix with four columns of min age, max age, longitude, and latitude, in that order, and rows that are entries for fossil occurrences.
#' @param which.biovars A vector of the numbers of the bioclimate variables that should be returned. The bioclimate variables number correspond to the Hijmans table at (https://www.worldclim.org/data/bioclim.html).
#' @param path path to the directory where the results should be saved
#' @param use.paleoclimate if left blank, default North America paleoclimate data is used. If FALSE, user submitted paleoclimate must be provided
#' @param paleoclimateUser list of data frames with paleoclimates, must be dataframes with columns: GlobalID, Longitude, Latitude, bio1, bio2,...,bio19. (see \code{getBioclimvars()}).
#' @param layerAge vector with the ages of the paleoclimate dataframes, if using user submitted paleoclimate data
#' @return a trait gram for minimum and maximum of biovariables
#' @seealso plotTraitGramMultiPhylo
#' @importFrom utils data
#' @importFrom ape dist.nodes
#' @importFrom grDevices pdf
#' @importFrom graphics lines
#' @importFrom graphics points
#' @importFrom grDevices dev.off
#' @export
#' @author A. Michelle Lawing, Alexandra F. C. Howard
#' @examples
#' data(sampletrees)
#' sampletrees <- sample(sampletrees,5)
#' data(occurrences)
#' occu <- getBioclimVars(occurrences, which.biovars=c(1,2))
#' sp_data_min<-sapply(4:5,function(x) tapply(occu[,x],occu$Species,min)) 
#' sp_data_max<-sapply(4:5,function(x) tapply(occu[,x],occu$Species,max))
#' ex_min <- geiger::treedata(sampletrees[[1]], sp_data_min, sort=TRUE)
#' ex_max <- geiger::treedata(sampletrees[[1]], sp_data_max, sort=TRUE)
#' colnames(ex_min$data)<- colnames(ex_max$data)<-c("bio1","bio2")  #labels biovars
#' \donttest{nodeest<- nodeEstimateEnvelopes(treedata_min=ex_min,treedata_max=ex_max, 
#' model="BM",which.biovars=c(1,2),
#' bounds=list(sigsq = c(min = 0, max = 1000000)))
#' plotTraitGram(ex_min, ex_max, nodeest$est, which.biovars=c(1,2), 
#' path=paste(tempdir(),"/",sep=""))}

plotTraitGram <- function(treedata_min, treedata_max, node_est, fossils=FALSE, which.biovars, path="", use.paleoclimate=TRUE, paleoclimateUser=NULL, layerAge=c(0:20)){
  num_traits <- length(treedata_min$data[1,])
  num_species <- length(treedata_min$data[,1])
  nodes <- (2*num_species)-2
  if(use.paleoclimate) {
    utils::data(paleoclimate) #uses paleoclimate data from package
  } else {
    if(is.null(paleoclimateUser)) {
      stop("paleoclimateUser argument must be provided when use.paleoclimate is FALSE.") #uses user inputted paleoclimate
    } else {
      paleoclimate <- paleoclimateUser
    }
  }
  for(i in 1:num_traits){
    M <- ape::dist.nodes(treedata_min$phy)[1,num_species+1]-ape::dist.nodes(treedata_min$phy)[,num_species+1]
    temp <- array(unlist(node_est),dim=c(2,length(node_est[[1]][1,,1]),num_traits,length(node_est)))
    traitgram_min_min <- cbind(c(treedata_min$data[,i],sapply(1:(num_species-1),function(x) min(temp[1,x,i,]))),M)
    traitgram_min_max <- cbind(c(treedata_min$data[,i],sapply(1:(num_species-1),function(x) max(temp[1,x,i,]))),M)
    traitgram_max_min <- cbind(c(treedata_max$data[,i],sapply(1:(num_species-1),function(x) min(temp[2,x,i,]))),M)
    traitgram_max_max <- cbind(c(treedata_max$data[,i],sapply(1:(num_species-1),function(x) max(temp[2,x,i,]))),M)
    traitgram_min <- cbind(c(treedata_min$data[,i],colMeans(t(temp[1,1:(num_species-1),i,]))),M)
    traitgram_max <- cbind(c(treedata_max$data[,i],colMeans(t(temp[2,1:(num_species-1),i,]))),M)
    #check on the size of graph and size of text
    grDevices::pdf(paste(path,colnames(treedata_min$data)[i],".pdf",sep=""),width=7,height=7,pointsize=10,useDingbats=F)
    plot(traitgram_min,ylim=rev(range(c(traitgram_min_min[,2]))),xlim=range(c(traitgram_min_min-10, traitgram_max_max+10)),type="n",xlab=colnames(treedata_min$data)[i],ylab="Time (mya)")
    for(j in 1:nodes) {lines(c(traitgram_max_min[j,1],traitgram_max_max[j,1]),c(traitgram_max_min[j,2],traitgram_max_min[j,2]),col="skyblue1",lwd=6)}
    for(j in 1:nodes) {lines(c(traitgram_min_min[j,1],traitgram_min_max[j,1]),c(traitgram_min_min[j,2],traitgram_min_min[j,2]),col="grey",lwd=4)}
    for(j in 1:nodes) {lines(traitgram_max[treedata_max$phy$edge[j,],1],traitgram_max[treedata_max$phy$edge[j,],2],col="skyblue2")}
    for(j in 1:nodes) {lines(traitgram_min[treedata_min$phy$edge[j,],1],traitgram_min[treedata_min$phy$edge[j,],2],col="darkgrey")}
    if(length(fossils)!=1){
      minbiofossils <- getBioclimVars(fossils[,c(1,3:4),drop=F], which.biovars = which.biovars, use.paleoclimate = use.paleoclimate, paleoclimateUser = paleoclimateUser, layerAge = layerAge)
      for(j in 1:length(minbiofossils[,1])) {points(minbiofossils[j,i+3],minbiofossils[j,1],col="black",pch=1)}
      maxbiofossils <- getBioclimVars(fossils[,c(2:4),drop=F], which.biovars = which.biovars, use.paleoclimate = use.paleoclimate, paleoclimateUser = paleoclimateUser, layerAge = layerAge)
      for(j in 1:length(maxbiofossils[,1])) {points(maxbiofossils[j,i+3],maxbiofossils[j,1],col="black",pch=16)}
    }
    dev.off()
  }
}
