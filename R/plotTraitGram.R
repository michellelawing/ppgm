#' @title plotTraitGram
#' @description Combine the node estimates based on random or specified fossil placement and plot them on a phylotraitgram in a specified directory.
#' @usage plotTraitGram(treedata_min, treedata_max, node_est, fossils=FALSE, which.biovars, path="")
#' @param treedata_min a tree data object with the min estimate of the climate envelope
#' @param treedata_max a tree data object with the max estimate of the climate envelope
#' @param node_est the estimate of all the nodes, both min and max
#' @param which.biovars A vector of the numbers of the bioclimate variables that should be returned. The bioclimate variables number correspond to the Hijmans table at (https://www.worldclim.org/data/bioclim.html).
#' @param path path to the directory where the results should be saved
#' @param use.paleoclimate if left blank, default North America paleoclimate data is used. If FALSE, user submitted paleoclimate must be provided
#' @param paleoclimateUser list of data frames with paleoclimates, must be dataframes with columns: GlobalID, Longitude, Latitude, bio1, bio2,...,bio19. (see \code{getBioclimvars()}).
#' @details
#' @return
#' @seealso plotTraitGramMultiPhylo
#' @export
#' @author A. Michelle Lawing, Alexandra F. C. Howard
#' @examples


plotTraitGram <- function(treedata_min, treedata_max, node_est, fossils=FALSE, which.biovars, path="", use.paleoclimate=TRUE, paleoclimateUser=NULL){
  require(ape)
  num_traits <- length(treedata_min$data[1,])
  num_species <- length(treedata_min$data[,1])
  if(use.paleoclimate) {
    data(paleoclimate) #uses paleoclimate data from package
  } else {
    if(is.null(paleoclimateUser)) {
      stop("paleoclimateUser argument must be provided when use.paleoclimate is FALSE.") #uses user inputted paleoclimate
    } else {
      paleoclimate <- paleoclimateUser
    }
  }
  for(i in 1:num_traits){
    M <- dist.nodes(treedata_min$phy)[1,num_species+1]-dist.nodes(treedata_min$phy)[,num_species+1]
    temp <- array(unlist(node_est),dim=c(2,length(node_est[[1]][1,,1]),num_traits,length(node_est)))
    traitgram_min_min <- cbind(c(treedata_min$data[,i],sapply(1:(num_species-1),function(x) min(temp[1,x,i,]))),M)
    traitgram_min_max <- cbind(c(treedata_min$data[,i],sapply(1:(num_species-1),function(x) max(temp[1,x,i,]))),M)
    traitgram_max_min <- cbind(c(treedata_max$data[,i],sapply(1:(num_species-1),function(x) min(temp[2,x,i,]))),M)
    traitgram_max_max <- cbind(c(treedata_max$data[,i],sapply(1:(num_species-1),function(x) max(temp[2,x,i,]))),M)
    traitgram_min <- cbind(c(treedata_min$data[,i],colMeans(t(temp[1,1:(num_species-1),i,]))),M)
    traitgram_max <- cbind(c(treedata_max$data[,i],colMeans(t(temp[2,1:(num_species-1),i,]))),M)
    #check on the size of graph and size of text
    pdf(paste(path,"bio",i,".pdf",sep=""),width=7,height=7,pointsize=10,useDingbats=F)
    plot(traitgram_min,ylim=rev(range(c(traitgram_min_min[,2]))),xlim=range(c(paleoclimate[[1]][,which.biovars[i]+3],
                                                                              paleoclimate[[11]][,which.biovars[i]+3], paleoclimate[[16]][,which.biovars[i]+3], traitgram_min_min,traitgram_max_max))
         ,type="n",xlab=colnames(treedata_min$data)[i],ylab="Time (mya)")
    lines(c(min(paleoclimate[[1]][,which.biovars[i]+3],traitgram_min),max(paleoclimate[[1]][,which.biovars[i]+3],traitgram_max)),
          c(0,0),col="antiquewhite",lwd=10)
    lines(c(min(paleoclimate[[6]][,which.biovars[i]+3]),max(paleoclimate[[6]][,which.biovars[i]+3])),c(5,5),col="antiquewhite",lwd=10)
    lines(c(min(paleoclimate[[11]][,which.biovars[i]+3]),max(paleoclimate[[11]][,which.biovars[i]+3])),c(10,10),col="antiquewhite",lwd=10)
    lines(c(min(paleoclimate[[16]][,which.biovars[i]+3]),max(paleoclimate[[16]][,which.biovars[i]+3])),c(15,15),col="antiquewhite",lwd=10)
    for(j in 1:104) {lines(c(traitgram_max_min[j,1],traitgram_max_max[j,1]),c(traitgram_max_min[j,2],traitgram_max_min[j,2]),col="skyblue1",lwd=6)}
    for(j in 1:104) {lines(c(traitgram_min_min[j,1],traitgram_min_max[j,1]),c(traitgram_min_min[j,2],traitgram_min_min[j,2]),col="grey",lwd=4)}
    for(j in 1:104) {lines(traitgram_max[treedata_max$phy$edge[j,],1],traitgram_max[treedata_max$phy$edge[j,],2],col="skyblue2")}
    for(j in 1:104) {lines(traitgram_min[treedata_min$phy$edge[j,],1],traitgram_min[treedata_min$phy$edge[j,],2],col="darkgrey")}
    if(length(fossils)!=1){for(j in 1:length(fossils[,1])) {points(fossils[j,i+3],fossils[j,1],col="black",pch=16)}}
    dev.off()
  }
}