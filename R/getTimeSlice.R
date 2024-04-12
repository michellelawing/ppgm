#' @title getTimeSlice
#' @description This function extracts estimated ancestral reconstructions for continuous characters any time specified along a phylogeny for all lineages present at the specified time.
#' @usage getTimeSlice(timeSlice, tree, trait, model = "BM", plot.est = FALSE)
#' @param timeSlice single numeric or a vector with the time (or times) to extract the estimated ancestor reconstructions.
#' @param tree an object of the class "phylo" that should be dated
#' @param trait a vector of both tip values and node estimates that correspond to tree
#' @param model if model = "estimate", the best fit model of evolution. If the model was specified, then model is the specified model, passes to \code{fitContinuous()}. Model options currently supported are: "BM", "OU", "EB", "lambda", "kappa", "delta"
#' @param plot.est a conditional stating whether or not to plot the results
#' @details The estimated reconstruction relies on an interpolation between node or between tip and node estimates of the trait. This method assumes a constant rate of evolution along the lineage where the interpolation is taking place.
#' @return \code{edge} for each time specified, a vector of edges that are present during that time are returned
#' @return \code{est} for each time specified, a vector of estimates of the ancestral reconstruction along each edge
#' @seealso \code{fitContinuous()}, \code{nodeEstimate()}
#' @author A. Michelle Lawing, Alexandra F. C. Howard
#' @importFrom ape dist.nodes
#' @importFrom geiger treedata
#' @importFrom graphics points
#' @export
#' @examples
#' data(sampletrees)
#' data(occurrences)
#' occurrences <- getBioclimVars(occurrences, which.biovars=1)
#' sp_data_min<- tapply(occurrences[,4],occurrences$Species,min)
#' treedata_min <- geiger::treedata(sampletrees[[1]], sp_data_min)
#' \dontrun{ex_est <- nodeEstimate(treedata_min, 1, model = 'BM') #runs BM model
#' ex_timeSlice <- getTimeSlice(10,treedata_min$phy,c(treedata_min$data[,1],ex_est$est))}


getTimeSlice <- function(timeSlice, tree, trait, model = "BM", plot.est=FALSE){
  M <- dist.nodes(tree)
  treeage <- max(M)/2
  maxedge <- (as.numeric(treeage - M[tree$edge[,1],tree$edge[1,1]]))
  minedge <- (as.numeric(treeage - M[tree$edge[,2],tree$edge[1,1]]))
  edgesample <- list()
  for(i in 1:length(timeSlice)){
    edgesample[[i]] <- which(maxedge>=timeSlice[i] & minedge<=timeSlice[i]+0.001)
  }
  hldtime <- list()
  for(i in 1:length(timeSlice)){
    hld <- array(NA, dim=length(edgesample[[i]]))
    for(j in 1:length(edgesample[[i]])){
      percentchange <- (maxedge[edgesample[[i]][j]]-timeSlice[i])/(maxedge[edgesample[[i]][j]]-minedge[edgesample[[i]][j]])
      hld[j] <- trait[tree$edge[edgesample[[i]][j],1]] - percentchange*(trait[tree$edge[edgesample[[i]][j],1]]-trait[tree$edge[edgesample[[i]][j],2]])
    }
    hldtime[[i]] <- hld
  }
  if(plot.est){
    xx <- nodeEstimate(geiger::treedata(tree, trait[1:(tree$edge[1,1]-1)]),1, model = model, plot.est = plot.est)
    for(i in 1:length(timeSlice)){
      points(array(max(M)/2-timeSlice[i],dim=length(hldtime[[i]])),hldtime[[i]])
    }
  }
  return(list(edge=edgesample,est=hldtime))
}

