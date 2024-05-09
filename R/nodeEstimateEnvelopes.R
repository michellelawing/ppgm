#' @title nodeEstimateEnvelopes
#' @description This function estimates climate envelopes at nodes with the optional placement of fossils on randomly assigned or specified edges on a tree.
#' @usage nodeEstimateEnvelopes(treedata_min, treedata_max, fossils=FALSE, 
#' fossils.edges=FALSE, model="BM", bounds=list(), control=list(), 
#' use.paleoclimate = TRUE, paleoclimateUser = NULL, layerAge = c(0:20), 
#' which.biovars = which.biovars)
#' @param treedata_min tree data object with min estimate of the climate envelope – list where first object is phylogeny, and second object is array of species with climate data variables (species must match)
#' @param treedata_max tree data object with max estimate of the climate envelope
#' @param fossils a matrix with three columns of age, longitude, and latitude, in that order, and rows that are entries for fossil occurrences.
#' @param fossils.edges the edge number that the fossil occurs on
#' @param model the model of evolution to use in the ancestral state reconstruction. Options are "estimate", "BM", "OU", "EB", "lambda", "kappa", "delta".
#' @param bounds bounds used for the model, passes to \code{fitContinuous()}, uses default if none specified.
#' @param control setting used for optimization of the model likelihood. Passes to \code{fitContinuous()}.
#' @param use.paleoclimate if left blank, default North America paleoclimate data is used. If FALSE, user submitted paleoclimate must be provided
#' @param paleoclimateUser list of data frames with paleoclimates, must be dataframes with columns: GlobalID, Longitude, Latitude, bio1, bio2,...,bio19. (see \code{getBioclimvars()}).
#' @param layerAge vector with the ages of the paleoclimate dataframes, if using user submitted paleoclimate data
#' @param which.biovars A vector of the numbers of the bioclimate variables that should be returned. The bioclimate variables number correspond to the Hijmans table at (https://www.worldclim.org/data/bioclim.html).
#' @details function adds fossils to trees according to \code{addFossil()}, then passes to \code{nodeEstimate()}.
#' @return an object of the class "nodeEstimate".
#' @return \code{model}    if model = "estimate", the best fit model of evolution. If the model was specified, then model is the specified model.
#' @return \code{est}      the ancestral node estimates of the continuous character.
#' @return \code{phy}      the phylogeny used for the estimates, which might be transformed depending on the evolutionary model.
#' @return \code{BM}       if model = "BM", returned values from \code{fitContinuous()} where the model is "BM"
#' @return \code{OU}       if model = "OU", returned values from \code{fitContinuous()} where the model is "OU"
#' @return \code{EB}       if model = "EB", returned values from \code{fitContinuous()} where the model is "EB"
#' @return \code{lambda}   if model = "lambda", returned values from \code{fitContinuous()} where the model is "lambda"
#' @return \code{kappa}    if model = "kappa", returned values from \code{fitContinuous()} where the model is "kappa"
#' @return \code{delta}    if model = "delta", returned values from \code{fitContinuous()} where the model is "delta"
#' @return \code{fitted}   if model = "estimate", returned values from the best fit model of evolution.
#' @seealso \code{nodeEstimate}, \code{fitContinuous}
#' @author A. Michelle Lawing, Alexandra F. C. Howard
#' @importFrom ape which.edge
#' @importFrom ape node.depth.edgelength
#' @importFrom geiger treedata
#' @importFrom ape drop.tip
#' @importFrom stringi stri_detect_fixed
#' @export
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
#' bounds=list(sigsq = c(min = 0, max = 1000000)))}


nodeEstimateEnvelopes <- function(treedata_min, treedata_max, fossils=FALSE, fossils.edges=FALSE, model="BM", bounds=list(), control=list(), use.paleoclimate = TRUE, paleoclimateUser = NULL, layerAge = c(0:20), which.biovars = which.biovars){
  num_species<-length(treedata_min$data[,1])
  num_traits<-length(treedata_min$data[1,])
  node<-array(NA,dim=c(2,num_species-1,num_traits))
  min_model<-as.list(array(NA,dim=num_traits))
  max_model<-as.list(array(NA,dim=num_traits))
  c.names<-colnames(treedata_min$data)
  #include a range around the fossil climate that is similar to extant species ranges
  if(length(treedata_max$data[1,])==1){
    addRangeFossils <- as.integer(mean((treedata_max$data - treedata_min$data)/4))
  } else {
    addRangeFossils <- as.integer(colMeans((treedata_max$data - treedata_min$data)/4))
  }
  #randomly assigns fossils and extracts biovariables for time period
  if(!is.logical(fossils)){
    biofossils <- array(NA,dim=c(length(fossils[,1]),3))
    colnames(biofossils) <- c("Age","Long","Lat")
    for(i in 1:length(fossils[,1])){
      if(is.logical(fossils.edges)) {f.edge <- NA}
      if(!is.logical(fossils.edges)) {
        if(is.na(fossils.edges[i])){
          f.edge <- NA
        } else{
          f.edge <- ape::which.edge(treedata_min$phy,fossils.edges[i])
        }
      }
      fossiltree <- addFossil(treedata_min$phy, mintime = fossils[i,1], maxtime = fossils[i,2], name = paste0("fossil",i), edge = f.edge)
      fossilnodeage <- as.integer(max(node.depth.edgelength(fossiltree)) - node.depth.edgelength(fossiltree)[which(stringi::stri_detect_fixed(fossiltree$tip.label,paste0("fossil",i)))])
      biofossils[i,] <- c(fossilnodeage,fossils[i,3:4,drop=F])
      newbiofossils <- getBioclimVars(biofossils[i,,drop=F], which.biovars = which.biovars, use.paleoclimate = use.paleoclimate, paleoclimateUser = paleoclimateUser, layerAge = layerAge)
      extract <-newbiofossils[,-c(1:3),drop=F]
      rownames(extract) <- paste0("fossil",i)
      treedata_min <- geiger::treedata(fossiltree, rbind(treedata_min$data,extract-addRangeFossils),sort=TRUE,warnings=F)
      treedata_max <- geiger::treedata(fossiltree, rbind(treedata_max$data,extract+addRangeFossils),sort=TRUE,warnings=F)
    }
  }
  treedata_min$phy$node.label<-treedata_max$phy$node.label<-((length(treedata_min$phy$tip)+1):((length(treedata_min$phy$tip)*2)-1))
  #Here is the trait loop start
  for(trai in 1:num_traits){
    min_model[[trai]]<-nodeEstimate(treedata_min,trai,model=model,bounds=bounds,control=control,plot.est=FALSE)
    tmin<-min_model[[trai]]$est
    max_model[[trai]]<-nodeEstimate(treedata_max,trai,model=model,bounds=bounds,control=control,plot.est=FALSE)
    tmax<-max_model[[trai]]$est
    if(length(fossils)!=1){
      node[1,,trai]<-tmin[as.numeric(rownames(tmin)) %in% ape::drop.tip(treedata_min$phy,rownames(fossils))$node.label]
      node[2,,trai]<-tmax[as.numeric(rownames(tmax)) %in% ape::drop.tip(treedata_max$phy,rownames(fossils))$node.label]
    }
    if(length(fossils)==1){
      node[1,,trai]<-tmin
      node[2,,trai]<-tmax
    }
  }
  return(list(est=node,min_model=min_model,max_model=max_model))
}

