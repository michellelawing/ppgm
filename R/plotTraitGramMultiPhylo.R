#' @title plotTraitGramMultiPhylo
#' @description Combine the node estimates based on random or specified fossil placement and plot them on a phylotrait gram in a specified directory.
#' @usage plotTraitGramMultiPhylo(treedata_min, treedata_max, node_est, 
#' fossils=FALSE, use.paleoclimate=TRUE, paleoclimateUser=NULL, 
#' which.biovars, path="")
#' @param treedata_min tree data object with min estimate of the climate envelope
#' @param treedata_max tree data object with max estimate of the climate envelope
#' @param node_est the estimate of all the nodes, both min and max. Must be in format [[trees]][[permut]][2,species,trait]
#' @param fossils a matrix with four columns of age to the closest million year integer, longitude, and latitude, in that order, and rows that are entries for fossil occurrences.
#' @param use.paleoclimate if left blank, default North America paleoclimate data is used. If FALSE, user submitted paleoclimate must be provided
#' @param paleoclimateUser list of data frames with paleoclimates, must be dataframes with columns: GlobalID, Longitude, Latitude, bio1, bio2,...,bio19. (see \code{getBioclimvars()}).
#' @param which.biovars A vector of the numbers of the bioclimate variables that should be returned. The bioclimate variables number correspond to the Hijmans table at (https://www.worldclim.org/data/bioclim.html).
#' @param path path to the directory where the results should be saved
#' @details plots a traitgram over multiple phylogenetic trees
#' @return a trait gram for minimum and maximum of biovariables over a distribution of phylogenetic trees
#' @importFrom utils data
#' @importFrom ape dist.nodes
#' @importFrom grDevices rgb
#' @seealso plotTraitGram
#' @export
#' @author A. Michelle Lawing, Alexandra F. C. Howard
#' @examples
#' data(sampletrees)
#' data(occurrences)
#' sampletrees <- sample(sampletrees,5)
#' biooccu <- getBioclimVars(occurrences, which.biovars=1)
#' sp_data_min<- tapply(biooccu[,4],biooccu$Species,min)
#' sp_data_max<- tapply(biooccu[,4],biooccu$Species,max)
#' treedata_min <- treedata_max <- node_est <- envelope <- list()
#' \donttest{for (tr in 1:length(sampletrees)){
#'   treedata_min[[tr]] <- geiger::treedata(sampletrees[[tr]],sp_data_min,sort=TRUE,warnings=F)
#'   treedata_max[[tr]] <- geiger::treedata(sampletrees[[tr]],sp_data_max,sort=TRUE,warnings=F)
#'   full_est <- nodeEstimateEnvelopes(treedata_min[[tr]],treedata_max[[tr]])
#'   node_est[[tr]]<-lapply(1, function(p) full_est$est)
#'   envelope[[tr]] <- getEnvelopes(treedata_min[[tr]], treedata_max[[tr]], node_est[[tr]])
#' }
#' plotTraitGramMultiPhylo(treedata_min,treedata_max, node_est, which.biovars=1, 
#' path=paste(tempdir(),"/",sep=""))
#' }


plotTraitGramMultiPhylo <- function(treedata_min, treedata_max, node_est, fossils=FALSE, use.paleoclimate=TRUE, paleoclimateUser=NULL, which.biovars, path=""){
  alpha.trans=as.integer(255/(1+log(length(node_est))))
  num_traits<-length(treedata_min[[1]]$data[1,])
  num_species<-length(treedata_min[[1]]$data[,1])
  if(use.paleoclimate) {
    paleoclimate <- paleoclimate #uses paleoclimate data from package
  } else {
    if(is.null(paleoclimateUser)) {
      stop("paleoclimateUser argument must be provided when use.paleoclimate is FALSE.") #uses user inputted paleoclimate
    } else {
      paleoclimate <- paleoclimateUser
    }
  }
  traitgram_min_min<-as.list(array(NA,dim=length(node_est)))
  traitgram_min_max<-as.list(array(NA,dim=length(node_est)))
  traitgram_max_min<-as.list(array(NA,dim=length(node_est)))
  traitgram_max_max<-as.list(array(NA,dim=length(node_est)))
  traitgram_min<-as.list(array(NA,dim=length(node_est)))
  traitgram_max<-as.list(array(NA,dim=length(node_est)))
  for(i in 1:num_traits){
    for(j in 1:length(node_est)){
      M<-ape::dist.nodes(treedata_min[[j]]$phy)[1,num_species+1]-ape::dist.nodes(treedata_min[[j]]$phy)[,num_species+1]
      temp <- array(unlist(node_est[[j]]),dim=c(2,length(node_est[[j]][[1]][1,,1]),num_traits,length(node_est[[j]])))
      traitgram_min_min[[j]]<-cbind(c(treedata_min[[j]]$data[,i],sapply(1:(num_species-1),function(x) min(temp[1,x,i,]))),M)
      traitgram_min_max[[j]]<-cbind(c(treedata_min[[j]]$data[,i],sapply(1:(num_species-1),function(x) max(temp[1,x,i,]))),M)
      traitgram_max_min[[j]]<-cbind(c(treedata_max[[j]]$data[,i],sapply(1:(num_species-1),function(x) min(temp[2,x,i,]))),M)
      traitgram_max_max[[j]]<-cbind(c(treedata_max[[j]]$data[,i],sapply(1:(num_species-1),function(x) max(temp[2,x,i,]))),M)
      traitgram_min[[j]]<-cbind(c(treedata_min[[j]]$data[,i],colMeans(t(temp[1,1:(num_species-1),i,]))),M)
      traitgram_max[[j]]<-cbind(c(treedata_max[[j]]$data[,i],colMeans(t(temp[2,1:(num_species-1),i,]))),M)
    }
    #check on the size of graph and size of text
    grDevices::pdf(paste(path,colnames(treedata_min[[1]]$data)[i],".pdf",sep=""),width=7,height=7,pointsize=10,useDingbats=F)
    plot(traitgram_min[[1]],ylim=rev(range(c(traitgram_min_min[[1]][,2]))),xlim=range(c(paleoclimate[[1]][,which.biovars[i]+3], paleoclimate[[11]][,which.biovars[i]+3], paleoclimate[[16]][,which.biovars[i]+3], traitgram_min_min,traitgram_max_max)), type="n",xlab=colnames(treedata_min[[1]]$data)[i],ylab="Time (mya)")
    graphics::lines(c(min(paleoclimate[[1]][,which.biovars[i]+3],traitgram_min[[j]]),max(paleoclimate[[1]][,which.biovars[i]+3],traitgram_max[[j]])),
          c(0,0),col="antiquewhite",lwd=10)
    graphics::lines(c(min(paleoclimate[[6]][,which.biovars[i]+3]),max(paleoclimate[[6]][,which.biovars[i]+3])),c(5,5),col="antiquewhite",lwd=10)
    graphics::lines(c(min(paleoclimate[[11]][,which.biovars[i]+3]),max(paleoclimate[[11]][,which.biovars[i]+3])),c(10,10),col="antiquewhite",lwd=10)
    graphics::lines(c(min(paleoclimate[[16]][,which.biovars[i]+3]),max(paleoclimate[[16]][,which.biovars[i]+3])),c(15,15),col="antiquewhite",lwd=10)
    for(j in 1:length(node_est)){
      for(k in 1:104) {lines(c(traitgram_max_min[[j]][k,1],traitgram_max_max[[j]][k,1]),c(traitgram_max_min[[j]][k,2],traitgram_max_min[[j]][k,2]),col=grDevices::rgb(135,206,235,alpha=alpha.trans,maxColorValue=255),lwd=6)}
      for(k in 1:104) {lines(c(traitgram_min_min[[j]][k,1],traitgram_min_max[[j]][k,1]),c(traitgram_min_min[[j]][k,2],traitgram_min_min[[j]][k,2]),col=grDevices::rgb(128,128,128,alpha=alpha.trans,maxColorValue=255),lwd=4)}
      for(k in 1:104) {lines(traitgram_max[[j]][treedata_max[[j]]$phy$edge[k,],1],traitgram_max[[j]][treedata_max[[j]]$phy$edge[k,],2],col=grDevices::rgb(135,206,235,alpha=alpha.trans,maxColorValue=255))}
      for(k in 1:104) {lines(traitgram_min[[j]][treedata_min[[j]]$phy$edge[k,],1],traitgram_min[[j]][treedata_min[[j]]$phy$edge[k,],2],col=grDevices::rgb(128,128,128,alpha=alpha.trans,maxColorValue=255))}
    }
    if(length(fossils)!=1){for(k in 1:length(fossils[,1])) {points(fossils[k,i+3],fossils[k,1],col="black",pch=16)}}
    dev.off()
  }
}
