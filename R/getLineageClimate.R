#' @title getLineageClimate
#' @description This function calculates the suitable climate for each specific lineage, starting at the tips and going back through time to the root.
#' @usage getLineageClimate(envelope, tree, which.biovars, 
#' use.paleoclimate=TRUE, paleoclimateUser=NULL)
#' @param envelope the min and max climate envelope of each lineage for each time slice, as outputted by \code{getEnvelopes()}
#' @param tree the phylogeny of all species. An object of class phylo
#' @param which.biovars a vector of the numbers of the bioclimate variables to be included. The bioclimate variables number correspond to the table at (https://www.worldclim.org/data/bioclim.html).
#' @param use.paleoclimate if left blank, default North America paleoclimate data is used. If FALSE, user submitted paleoclimate must be provided
#' @param paleoclimateUser list of data frames with paleoclimates, must be dataframes with columns: GlobalID, Longitude, Latitude, bio1, bio2,...,bio19. (see \code{getBioclimvars()}).
#' @details Calculates rate of geographic change of all lineages. Outputs both the geographic center change, and the geographic size change.
#' @return \code{matchedClim} list of occurrences points for each lineage, for each time slice of paleoclimate data
#' @return \code{lineage} list of lineage specific nodes, as output from phangorn::Ancestors
#' @author A. Michelle Lawing, Alexandra F. C. Howard, Maria A. Hurtado-Materon
#' @seealso \code{getEnvelopes()} \code{getGeoRate()}
#' @importFrom phangorn Ancestors
#' @export
#' @examples
#' data(sampletrees)
#' data(occurrences)
#' data(paleoclimate)
#' occu <- getBioclimVars(occurrences, which.biovars=1)
#' tree <- sampletrees[[25]]
#' #species minimum for biovariable 1
#' sp_data_min<- tapply(occu[,4],occu$Species,min)
#' #species maximum for biovariable 1
#' sp_data_max<- tapply(occu[,4],occu$Species,max)
#' #convert to treedata object
#' treedata_min <- geiger::treedata(tree,sp_data_min,sort=TRUE,warnings=F) 
#' treedata_max <- geiger::treedata(tree,sp_data_max,sort=TRUE,warnings=F)
#' #estimate node values using Brownian Motion
#' \dontrun{full_est <- nodeEstimateEnvelopes(treedata_min,treedata_max)
#' node_est <- full_est$est #extract only node estimates
#' #calculate climate envelopes
#' example_getEnvelopes <- getEnvelopes(treedata_min, treedata_max, node_est)
#' #calculate lineage specific climate
#' example_getLinClim <- getLineageClimate(example_getEnvelopes, tree, which.biovars=1)}

getLineageClimate <- function(envelope, tree, which.biovars, use.paleoclimate=TRUE, paleoclimateUser=NULL){
  if(use.paleoclimate) {
    paleoclimate <- paleoclimate #uses paleoclimate data from package
  } else {
    if(is.null(paleoclimateUser)) {
      stop("paleoclimateUser argument must be provided when use.paleoclimate is FALSE.") #uses user inputted paleoclimate
    } else {
      paleoclimate <- paleoclimateUser
    }
  }
  temp_min <- lapply(1:length(paleoclimate),function(i){
    temp <- lapply(1:length(which.biovars),function(j){getTimeSlice(i-1,tree,envelope[,2,j])})
    temp <- t(array(unlist(temp),dim=c(length(unlist(temp[[1]]$edge)),2*length(which.biovars))))
    return(temp)})
  temp_max <- lapply(1:length(paleoclimate),function(i){
    temp <- lapply(1:length(which.biovars),function(j){getTimeSlice(i-1,tree,envelope[,5,j])})
    temp <- t(array(unlist(temp),dim=c(length(unlist(temp[[1]]$edge)),2*length(which.biovars))))
    return(temp)})
  lineage <- lapply(1:length(temp_min[[1]][1,]),function(i) {c(i,phangorn::Ancestors(tree,i))})
  matchedClim <- list()
  #for every lineage, for every time slice, which lineages are present and then which points is lineage present
  #matchedClim: what points is that lineage climate present
  for(i in 1:length(temp_min[[1]][1,])){
    tempmatch <- list()
    for(j in 1:length(paleoclimate)){
      temp<-match(match(lineage[[i]], tree$edge[,2]),temp_min[[j]][1,])
      temp<-temp[!is.na(temp)]
      matching<-
        sapply(1:length(which.biovars),function(x){paleoclimate[[j]][,which.biovars[x]+3,drop=F]>=temp_min[[j]][1:length(which.biovars)*2,temp][x]&paleoclimate[[j]][ ,which.biovars[x]+3,drop=F]<=temp_max[[j]][1:length(which.biovars)*2,temp][x]})
      matching<-which(rowSums(matching)==length(which.biovars),arr.ind=TRUE)
      tempmatch[[j]] <- matching
    }
    matchedClim[[i]] <- tempmatch
  }
  #name outputs for clarity
  names(matchedClim) <- names(lineage) <- tree$tip.label
  return(list(matchedClim=matchedClim, lineage=lineage))
}
