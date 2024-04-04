#' @title getEnvelopes
#' @description This function gets the bioclimate envelopes of species and nodes.
#' @usage getEnvelopes(treedata_min, treedata_max, node_est)
#' @param treedata_min tree data object with min estimate of the climate envelope for each species
#' @param treedata_max tree data object with max estimate of the climate envelope for each species
#' @param node_est the estimate of all the nodes, both min and max
#' @details Function derives the minimum, and maximum of each climate variable
#' @return An array containing climate envelopes for each node
#' @author A. Michelle Lawing, Alexandra F. C. Howard
#' @seealso \code{ppgmMESS()}, \code{nodeEstimate}
#' @importFrom ape dist.nodes
#' @export
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

getEnvelopes <- function(treedata_min,treedata_max,node_est){
  num_traits <- length(treedata_min$data[1,])
  num_species <- length(treedata_min$data[,1])
  envelope <- array(NA, dim=c(2*num_species-1,5,num_traits))
  for(i in 1:num_traits){
    M <- ape::dist.nodes(treedata_min$phy)[1,num_species+1]-ape::dist.nodes(treedata_min$phy)[,num_species+1]
    temp <- array(unlist(node_est),dim=c(2,num_species-1,num_traits,length(node_est)))
    traitgram_min_min <- cbind(c(treedata_min$data[,i],sapply(1:(num_species-1),function(x) min(temp[1,x,i,]))),M)
    traitgram_min_max <- cbind(c(treedata_min$data[,i],sapply(1:(num_species-1),function(x) max(temp[1,x,i,]))),M)
    traitgram_max_min <- cbind(c(treedata_max$data[,i],sapply(1:(num_species-1),function(x) min(temp[2,x,i,]))),M)
    traitgram_max_max <- cbind(c(treedata_max$data[,i],sapply(1:(num_species-1),function(x) max(temp[2,x,i,]))),M)
    traitgram_min <- cbind(c(treedata_min$data[,i],colMeans(t(temp[1,1:(num_species-1),i,]))),M)
    traitgram_max <- cbind(c(treedata_max$data[,i],colMeans(t(temp[2,1:(num_species-1),i,]))),M)
    envelope[,,i] <- cbind(traitgram_min_min[,c(2,1)],traitgram_min_max[,1],traitgram_max_min[,1],traitgram_max_max[,1])
  }
  return(envelope)
}
