#' @title getGeoRate
#' @description This function calculates the change in suitable habitat through time in geographic space.
#' @usage getGeoRate(envelope, tree, which.biovars, use.paleoclimate=TRUE, 
#' paleoclimateUser=NULL)
#' @param envelope the min and max climate envelope of each lineage for each time slice, as outputted by \code{getEnvelopes()}
#' @param tree the phylogeny of all species. An object of class phylo
#' @param which.biovars a vector of the numbers of the bioclimate variables to be included. The bioclimate variables number correspond to the table at (https://www.worldclim.org/data/bioclim.html).
#' @param use.paleoclimate if left blank, default North America paleoclimate data is used. If FALSE, user submitted paleoclimate must be provided
#' @param paleoclimateUser list of data frames with paleoclimates, must be dataframes with columns: GlobalID, Longitude, Latitude, bio1, bio2,...,bio19. (see \code{getBioclimvars()}).
#' @details Calculates rate of geographic change of all lineages. Outputs both the geographic center change, and the geographic size change.
#' @return \code{geo_center} change in geographic center of suitable climate envelope
#' @return \code{geo_size} change in geographic size of suitable climate envelope
#' @return \code{time_int} time intervals
#' @author A. Michelle Lawing, Alexandra F. C. Howard, Maria A. Hurtado-Materon
#' @importFrom fields rdist.earth
#' @importFrom phangorn Ancestors
#' @importFrom stats dist
#' @export
#' @seealso \code{getEnvelopes()}
#' @examples
#' data(sampletrees)
#' data(occurrences)
#' data(paleoclimate)
#' tree <- sampletrees[[25]]
#' occu <- getBioclimVars(occurrences, which.biovars=1)
#' sp_data_min<- tapply(occu[,4],occu$Species,min)
#' sp_data_max<- tapply(occu[,4],occu$Species,max)
#' treedata_min <- geiger::treedata(tree,sp_data_min,sort=TRUE,warnings=F)
#' treedata_max <- geiger::treedata(tree,sp_data_max,sort=TRUE,warnings=F)
#' \dontrun{full_est <- nodeEstimateEnvelopes(treedata_min,treedata_max)
#' node_est <- full_est$est
#' example_getEnvelopes <- getEnvelopes(treedata_min, treedata_max, node_est)
#' example_getGeoRate <- getGeoRate(example_getEnvelopes, tree,which.biovars=1)}

getGeoRate <- function(envelope, tree, which.biovars, use.paleoclimate=TRUE, paleoclimateUser=NULL){
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
  geo_center <- array(NA, dim=c(length(paleoclimate),length(lineage),2))
  geo_size <- array(NA, dim=c(length(paleoclimate),length(lineage)))
  #for every lineage, for every time slice, which lineages are present and then which points is lineage present
  for(i in 1:length(temp_min[[1]][1,])){
    for(j in 1:length(paleoclimate)){
      temp<-match(match(lineage[[i]], tree$edge[,2]),temp_min[[j]][1,])
      temp<-temp[!is.na(temp)]
      matching<-
        sapply(1:length(which.biovars),function(x){paleoclimate[[j]][,which.biovars[x]+3,drop=F]>=temp_min[[j]][1:length(which.biovars)*2,temp][x]&paleoclimate[[j]][ ,which.biovars[x]+3,drop=F]<=temp_max[[j]][1:length(which.biovars)*2,temp][x]})
      matching<-which(rowSums(matching)==length(which.biovars),arr.ind=TRUE)
      geo_center[j,i,]<-colMeans(paleoclimate[[j]][matching,2:3])
      geo_size[j,i]<-dim(paleoclimate[[j]][matching,2:3])[1]*50*50
    }
  }
  #matching: what points is that lineage present
  #calculate the geo dist from centers and expansions/contractions from areas
  C <- array(NA, dim=c(length(lineage),length(paleoclimate),length(paleoclimate)))
  A <- array(NA, dim=c(length(lineage),length(paleoclimate),length(paleoclimate)))
  for(i in 1:length(lineage)){
    C[i,,] <- fields::rdist.earth(geo_center[,i,],geo_center[,i,],miles=FALSE)
    A[i,,] <- sapply(1:length(paleoclimate),function(x){geo_size[x,i]-geo_size[,i]})
  }
  T <- as.matrix(stats::dist(1:length(paleoclimate),diag=TRUE))
  return(list(geo_center=C, geo_size=A, time_int=T))
}

