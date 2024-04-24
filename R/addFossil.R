#' @title addFossil
#' @description Adds a fossil as a tip to a specified phylogeny given either an age range that the fossil occurs in, a specific edge that the fossil diverged from, or both. If the specific edge placement for the fossil is unknown, then this function randomly places the fossil on any edge that is within the age range.
#' @usage addFossil(tree, mintime = 0, maxtime = NA, name = "fossil", edge = NA)
#' @param tree An object of the class "phylo"
#' @param mintime The minimum age of the fossil. If no minimum time is specified, the default value is 0.
#' @param maxtime The maximum age of the fossil. If no maximum time is specified, the default value is the maximum tree age.
#' @param name The name of the fossil to appear as a tip.label.
#' @param edge The edge on the tree where the fossil presumably diverged. If no edge is specified, then the function randomly selects an edge within the age range of the fossil.
#' @details There are several random components to this function. First, if an edge is not specified to place a fossil, then an edge is randomly selected that is within the age range of the fossil. Second, the exact placement of the node leading to the fossil is randomly selected within the age range specified. Third, the length of the edge leading to the fossil is randomly selected with constraints on the maximum length of the edge, where the maximum length of the edge cannot render the fossil younger than the minimum time of occurrence as specified in the mintime argument.
#' @return An object of the class "phylo".
#' @author A. Michelle Lawing, Alexandra F. C. Howard
#' @importFrom ape dist.nodes
#' @importFrom stats runif
#' @importFrom ape bind.tree
#' @importFrom phytools pbtree
#' @export
#' @examples
#' mytree <- phytools::pbtree(n=20)
#' newtree <- addFossil(mytree, mintime = max(mytree$edge.length)/2, maxtime= max(mytree$edge.length))
#' plot(newtree)

addFossil <- function(tree, mintime = 0, maxtime = NA, name = "fossil", edge = NA) {
  if(is.na(maxtime)){maxtime=max(ape::dist.nodes(tree))/2}
  tree$node.label <- ((length(tree$tip)+1):((length(tree$tip)*2)-1))
  treeage <- max(ape::dist.nodes(tree))/2
  M <- ape::dist.nodes(tree)
  maxedge <- (as.numeric(treeage - M[tree$edge[,1],tree$edge[1,1]]))
  minedge <- (as.numeric(treeage - M[tree$edge[,2],tree$edge[1,1]]))
  if(!is.na(edge)){edgesample <- edge}
  if(is.na(edge)){edgesample <- sample(which(maxedge>mintime & minedge<maxtime),1)}
  dedge <- tree$edge[edgesample,2]
  place <- stats::runif(1,max(c(minedge[edgesample],mintime)),min(c(maxtime,maxedge[edgesample])))
  fossil <- list(edge=matrix(c(2,1),1,2), tip.label=name, edge.length=stats::runif(1,min=0.0000000001,max=(place-max(c(minedge[edgesample],mintime)))), Nnode=1)
  class(fossil) <- "phylo"
  tree <- ape::bind.tree(tree,fossil,where=dedge,position=place-minedge[edgesample])
  tree$node.label <- as.numeric(tree$node.label)+1
  newnode <- which(is.na(tree$node.label))
  tree$node.label[(newnode+1):length(tree$node.label)] <- as.numeric(tree$node.label[(newnode+1):length(tree$node.label)])+1
  tree$node.label[newnode] <- as.numeric(tree$node.label[newnode-1])+1
  return(tree)
}
