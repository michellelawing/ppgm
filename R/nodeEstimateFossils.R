#' @title nodeEstimateFossils
#' @description This function estimates nodes with the placement of fossils on randomly assigned or specified edges on a tree.
#' @usage nodeEstimateFossils(treedata_min, treedata_max, fossils=FALSE, fossils.edges=FALSE, model="BM", bounds=list())
#' @param treedata_min tree data object with min estimate of the climate envelope – list where first object is phylogeny, and second object is array of species with climate data variables (species must match)
#' @param treedata_max tree data object with max estimate of the climate envelope
#' @param fossils the estimate of the climate envelope of the fossil occurrences
#' @param fossils.edges the edge number that the fossil occurs on
#' @param model the model of evolution to use in the ancestral state reconstruction. Options are "estimate", "BM", "OU", "EB", "lambda", "kappa", "delta".
#' @param bounds bounds used for the model, passes to \code{fitContinuous()}, uses default if none specified.
#' @param control setting used for optimization of the model likelihood. Passes to \code{fitContinuous()}.
#' @param … arguments to pass to \code{fitContinuous}.
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
#' @export
#' @examples
#' data(beastLeache)
#' data(occurrences)
#' data(new_fossils)
#' biofossils <- getBioclimVars(new_fossils,c(1,4,15))
#' rownames(biofossils)<-paste("fossil",1:length(biofossils[,1]),sep="")
#' sp_data_min<- tapply(occurrences[,4],occurrences$Species,min)
#' sp_data_max<- tapply(occurrences[,4],occurrences$Species,max)
#' ex_min <- treedata(beastLeache[[1]], sp_data_min)
#' ex_max <- treedata(beastLeache[[1]], sp_data_max)
#' example_nodeest<- nodeEstimateFossils(treedata_min=ex_min,treedata_max=ex_max,model="BM", fossils=biofossils, bounds=list(sigsq = c(min = 0, max = 1000000), SE = c(0, 0.1)))


nodeEstimateFossils <- function(treedata_min, treedata_max, fossils=FALSE, fossils.edges=FALSE, model="BM", bounds=list(), control=list(), ...){
  require(geiger)
  require(phangorn)
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
  if(!is.logical(fossils)){
    for(i in 1:length(fossils[,1])){
      if(is.logical(fossils.edges)) {f.edge <- NA}
      if(!is.logical(fossils.edges)) {
        if(is.na(fossils.edges[i])){
          f.edge <- NA
        } else{
          f.edge <- which.edge(treedata_min$phy,fossils.edges[i])
        }
      }
      treedata_min <- treedata(addFossil(treedata_min$phy, mintime = fossils[i,1], maxtime = fossils[i,1] + 1, name = rownames(fossils)[i],edge = f.edge), rbind(treedata_min$data,fossils[i,colnames(treedata_min$data),drop=F] - addRangeFossils),sort=TRUE,warnings=F)
      rownames(treedata_min$data)<-treedata_min$phy$tip.label
      colnames(treedata_min$data)<-c.names
      treedata_max <- treedata(addFossil(treedata_max$phy, mintime = fossils[i,1], maxtime = fossils[i,1] + 1, name = rownames(fossils)[i],edge = f.edge), rbind(treedata_max$data,fossils[i,colnames(treedata_max$data),drop=F] + addRangeFossils),sort=TRUE,warnings=F)
      rownames(treedata_max$data)<-treedata_max$phy$tip.label
      colnames(treedata_max$data)<-c.names
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
      node[1,,trai]<-tmin[as.numeric(rownames(tmin)) %in% drop.tip(treedata_min$phy,rownames(fossils))$node.label]
      node[2,,trai]<-tmax[as.numeric(rownames(tmax)) %in% drop.tip(treedata_max$phy,rownames(fossils))$node.label]
    }
    if(length(fossils)==1){
      node[1,,trai]<-tmin
      node[2,,trai]<-tmax
    }
  }
  return(list(est=node,min_model=min_model,max_model=max_model))
}

