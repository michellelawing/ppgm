#' @title addfossil
#' @description This function adds a fossil as a tip to a specified phylogeny given either an age range that the fossil occurs in, a specific edge that the fossil diverged from, or both. If the specific edge placement for the fossil is unknown, then this function randomly places the fossil on any edge that is within the age range.
#' @usage addfossil(tree, mintime = 0, maxtime = NA, name = "fossil", edge = NA)
#' @param tree An object of the class "phylo".
#' @param mintime The minimum age of the fossil. If no minimum time is specified, the default value is 0.
#' @param maxtime The maximum age of the fossil. If no maximum time is specified, the default value is the maximum tree age.
#' @param name The name of the fossil to appear as a tip.label.
#' @param edge The edge on the tree where the fossil presumably divergered. If no edge is specified, then the funciton randomly selects an edge within the age range of the fossil.
#' @details There are several random components to this function. First, if an edge is not specified to place a fossil, then an edge is randomly selected that is within the age range of the fossil. Second, the exact placement of the node leading to the fossil is randomly selected within the age range specified. Third, the length of the edge leading to the fossil is randomly selected with constraints on the maximum length of the edge, where the maximum length of the edge cannot render the fossil younger than the minimum time of occurrence as specified in the mintime arguement.
#' @return An object of the class "phylo".
#' @author A. Michelle Lawing
#' @seealso bind.tree
# @examples require(ape) #function and example depend on ape
# mytree <- rcoal(3)
# newtree <- addfossil(mytree, mintime = max(mytree$edge.length)/2, maxtime= max(mytree$edge.length))
# plot(newtree)

addfossil<- function(tree,mintime=0,maxtime=NA,name="fossil",edge=NA) {
  #function depends on ape
  require(ape)
  if(is.na(maxtime)){maxtime=max(dist.nodes(tree))/2}
  tree$node.label<-((length(tree$tip)+1):((length(tree$tip)*2)-1))
  treeage<-max(dist.nodes(tree))/2
  M<-dist.nodes(tree)
  maxedge<-(as.numeric(treeage - M[tree$edge[,1],tree$edge[1,1]]))
  minedge<-(as.numeric(treeage - M[tree$edge[,2],tree$edge[1,1]]))
  if(!is.na(edge)){edgesample<-edge}
  if(is.na(edge)){edgesample<-sample(which(maxedge>mintime & minedge<maxtime),1)}
  dedge<-tree$edge[edgesample,2]
  place<-runif(1,max(c(minedge[edgesample],mintime)),min(c(maxtime,maxedge[edgesample])))
  fossil<-list(edge=matrix(c(2,1),1,2), tip.label=name, edge.length=runif(1,min=0.0000000001,max=(place-max(c(minedge[edgesample],mintime)))), Nnode=1)
  class(fossil)<-"phylo"
  tree<-bind.tree(tree,fossil,where=dedge,position=place-minedge[edgesample])
  tree$node.label<-as.numeric(tree$node.label)+1
  newnode=which(is.na(tree$node.label))
  tree$node.label[(newnode+1):length(tree$node.label)]<-as.numeric(tree$node.label[(newnode+1):length(tree$node.label)])+1
  tree$node.label[newnode]<-as.numeric(tree$node.label[newnode-1])+1
  return(tree)
}

#' @title getBioclimVars
#' @description This function retrieves the bioclimate variables described in Nix (1986) for the specified variables and the specified time period. Modern time period uses the Hijmans et al. (2005) high resolution climate interpolations. The time period 10 Ma uses the GCM by XX for the Tortonian and the time period 15 Ma uses the GCM for the Lanhgian by XX. For the one million year intervals outside the modern and past GCMs, the climate was interpolated based on the benthic marine foram stable oxygen isotope ratio curve from XX. The scale of these variables is at a 50 km equidistant point grain size corresponding to Polly XX.
#' @usage getBioclimVars(occurrences, which.biovars=c(1,12))
#' @param occurrences A matrix or data.frame with three columns and rows to represent idividuals. The first column is either the species name or the age of the fossil. 
#' @param which.biovars A vector of the numbers of the bioclimate variables that should be returned. The bioclimate variables number correspond to the Hijmans table at (www.XX).
#' @details The occurrences argument should contain all extant or all fossils.
#' @return Returns a data.frame with the original occurrences input appended with columns of bioclimate variables specified. If fossils are included, the returned bioclimate variables are from the closest 1 Ma interval of isotopically scaled climate.
#' @author A. Michelle Lawing
#' @references Hijmans; Polly; GCMs
# @examples require(fields)

getBioclimVars<-function(occurrences,which.biovars=c(1,12)){
  require (fields)
  data(paleoclimate) #isotopically scaled paleoclimate bioclimate variables for North America
  if(dim(occurrences)[2]<3){
    return ("ERROR: There are less than three columns for the occurrences input. You need to have at least three (Species Name, Longitude, Latitude).")
  }
  if(dim(occurrences)[2]>3){
    occurrences<-occurrences[,c(1,2,3,which.biovars+3)]
    return (occurrences)
  }
  if(dim(occurrences)[2]==3){
    if(!is.numeric(occurrences[1,1])){
      calc_dist<-rdist.earth(paleoclimate[[1]][,2:3],occurrences[,2:3])
      min_dist<-1:length(occurrences[,1])
      for(i in 1:length(occurrences[,1])){
        min_dist[i]<-which.min(calc_dist[,i])
      }
      occurrences<-cbind(occurrences,paleoclimate[[1]][min_dist,which.biovars+3])
      return(occurrences)
    }
    if(is.numeric(occurrences[1,1])){
      temp<-array(NA, dim=c(length(occurrences[,1]),length(which.biovars)))
      for(i in 1:length(occurrences[,1])){
        anothertemp <- as.matrix(paleoclimate[[as.integer(occurrences[i,1])]][,2:3])
        calc_dist<-rdist.earth(anothertemp,t(as.matrix(occurrences[i,2:3])))  #'#'ERROR subscript out of bounds, hapopens when climate is extracted from occurrences
        temp[i,]<-unlist(paleoclimate[[as.integer(occurrences[i,1])]][which.min(calc_dist),which.biovars+3])
      }
      occurrences<-cbind(occurrences,temp)
      colnames(occurrences)<-c("Age", "Longitude", "Latitude", paste("bio",which.biovars,sep=""))
      return(occurrences)
    }
  }
}

#' @title node.estimate
#' @description This function estimates the ancestral character states for continuous characters given a model of evolution or using the best fit model of evolution from the fitContinuous function in the geiger pacakage. The ancestral states are estimated using GLS described in Martins and Hansen (1997).
#' @usage node.estimate(treedata.obj, traitnum, model = "BM", plot.est = TRUE)
#' @param treedata.obj an object of the class "treedata".
#' @param traitnum the column number of the trait within the treedata object to be reconstructed.
#' @param model the model of evolution to use in the ancestral state reconstruction. Options are c("estimate","BM","OU","EB","trend","lambda","kappa","delta","drift","white").
#' @param plot.est whether or not to plot the traitgram of the estimated ancester states.
#' @param ... arguments to pass to fitContinuous
#' @details See the fitContinuous details for descriptions of the models of evolution and their parameters being estimated.
#' @return an object of the class "node.estimate".
#' @return model    if model = "estimate", the best fit model of evolution. If the model was specified, then model is the specified model.
#' @return est      the ancestral node estimates of the continuous character.
#' @return phy      the phylogeny used for the estimates, which might be transformed depending on the evolutionary model.
#' @return BM       if model = "BM", returned values from the fitContinuous function where the model is "BM"
#' @return OU       if model = "OU", returned values from the fitContinuous function where the model is "OU"
#' @return EB       if model = "EB", returned values from the fitContinuous function where the model is "EB"
#' @return trend    if model = "trend", returned values from the fitContinuous function where the model is "trend"
#' @return lambda   if model = "lambda", returned values from the fitContinuous function where the model is "lambda"
#' @return kappa    if model = "kappa", returned values from the fitContinuous function where the model is "kappa"
#' @return delta    if model = "delta", returned values from the fitContinuous function where the model is "delta"
#' @return drift    if model = "drift", returned values from the fitContinuous function where the model is "drift"
#' @return white    if model = "white", returned values from the fitContinuous function where the model is "white"
#' @return fitted   if model = "estimate", returned values from the best fit model of evolution.
#' @seealso fitContinuous
#' @references Martins, E. P. and Hansen, T. F. (1997) Phylogenies and the comparative method: a general approach to incorporating phylogenetic information into the analysis of interspecific data. American Naturalist, 149, 646–667.
#' @author A. Michelle Lawing
# @examples
#ex <- treedata(trialestWF$model_min[[1]][[1]][[1]]$phy, trialestWF$cem[,c(1,5)])
#node.estimate(ex, 1, model = 'OU')

node.estimate <- function(treedata.obj,traitnum,model="BM",bounds=list(),control=list(),plot.est=FALSE,...) {
  require(geiger)
  x <- treedata.obj$data[,traitnum]
  phy <- treedata.obj$phy
  was.estimated <- FALSE
  fitted <- model
  if(model=="estimate"){
    models=c("BM","OU","EB","trend","lambda","kappa","delta","drift","white")
    BM=try(fitContinuous(phy,x,model="BM",bounds=bounds,control=control),silent=T)
    OU=try(fitContinuous(phy,x,model="OU",bounds=bounds,control=control),silent=T)
    EB=try(fitContinuous(phy,x,model="EB",bounds=bounds,control=control),silent=T)
    if(!is.ultrametric(phy)){trend=try(fitContinuous(phy,x,model="trend",bounds=bounds,control=control),silent=T)} else {trend=NA}
    lambda=try(fitContinuous(phy,x,model="lambda",bounds=bounds,control=control),silent=T)
    kappa=try(fitContinuous(phy,x,model="kappa",bounds=bounds,control=control),silent=T)
    delta=try(fitContinuous(phy,x,model="delta",bounds=bounds,control=control),silent=T)
    if(!is.ultrametric(phy)){drift=try(fitContinuous(phy,x,model="drift",bounds=bounds,control=control),silent=T)} else {drift=BM}
    white=try(fitContinuous(phy,x,model="white",bounds=bounds,control=control),silent=T)
    trait.macroevo <- list()
    for (mod in 1:length(models)){
      if(is(eval(parse(text=models[mod])),"try-error") | is(eval(parse(text=models[mod])),"NA")) {
        } else {
        trait.macroevo[[mod]] <- eval(parse(text=models[mod]))$opt$aicc
      }
    }
    model <- models[which.min(unlist(trait.macroevo))]
    fitted <- list()
    for (mod in 1:length(models)){
      if(is(eval(parse(text=models[mod])),"try-error") | is(eval(parse(text=models[mod])),"NA")) {
        fitted[[mod]] <- NA
      } else {
        fitted[[models[mod]]] <- eval(parse(text=models[mod]))$opt
      }
    }
    if (model=="OU" & !is(OU,"try-error")) {if(!is.na(OU$opt$alpha)){phy=rescale(phy,model="OU",OU$opt$alpha)}}
    if (model =="EB" & !is(EB,"try-error")){if(!is.na(EB$opt$a)){phy=rescale(phy,model="EB",EB$opt$a)}}
    #'if model is trend, tree should NOT be ultrametric
    if (model=="trend" & !is(trend,"try-error")) {if(!is.na(trend$opt$slope)){ phy=rescale(phy,model="trend",trend$opt$slope)}}
    if (model=="lambda" & !is(lambda,"try-error")) {if(!is.na(lambda$opt$lambda)){phy=rescale(phy,model="lambda",lambda$opt$lambda)}}
    if(model=="kappa" & !is(kappa,"try-error")) {if(!is.na(kappa$opt$kappa)){phy=rescale(phy,model="kappa",kappa$opt$kappa)}}
    if(model=="delta" & !is(delta,"try-error")) {if(!is.na(delta$opt$delta)){phy=rescale(phy,model="delta",delta$opt$delta)}}
    was.estimated <- TRUE
  }
  if(!was.estimated){
    if (model=="BM") {fitted<-BM<-try(fitContinuous(phy,x,model="BM",bounds=bounds,control=control),silent=T)}
    if (model=="OU") {fitted<-OU<-try(fitContinuous(phy,x,model="OU",bounds=bounds,control=control),silent=T)
       if(!is(OU,"try-error")){ if(!is.na(OU$opt$alpha)){phy=rescale(phy,model="OU",OU$opt$alpha)}}}
    if (model =="EB"){fitted<-EB<-try(fitContinuous(phy,x,model="EB",bounds=bounds,control=control),silent=T)
       if(!is(EB,"try-error")){ if(!is.na(EB$opt$a)){ phy=rescale(phy,model="EB",EB$opt$a)}}}
    #'if model is trend, tree should NOT be ultrametric
    if (model=="trend") {fitted<-trend<-try(fitContinuous(phy,x,model="trend",bounds=bounds,control=control),silent=T)
       if(!is(trend,"try-error")){ if(!is.na(trend$opt$slope)){ phy=rescale(phy,model="trend",trend$opt$slope)}}}
    if (model=="lambda") {fitted<-lambda<-try(fitContinuous(phy,x,model="lambda",bounds=bounds,control=control),silent=T)
       if(!is(lambda,"try-error")){if(!is.na(lambda$opt$lambda)){phy=rescale(phy,model="lambda",lambda$opt$lambda)}}}
    if(model=="kappa") {fitted<-kappa<-try(fitContinuous(phy,x,model="kappa",bounds=bounds,control=control),silent=T)
       if(!is(kappa,"try-error")){ if(!is.na(kappa$opt$kappa)){ phy=rescale(phy,model="kappa",kappa$opt$kappa)}}}
    if(model=="delta") {fitted<-delta<-try(fitContinuous(phy,x,model="delta",bounds=bounds,control=control),silent=T)
      if(!is(delta,"try-error")){ if(!is.na(delta$opt$delta)){phy=rescale(phy,model="delta",delta$opt$delta)}}}
    #'if model is drift, tree should NOT be ultrametric
    if (model=="drift") {fitted<-drift<-try(fitContinuous(phy,x,model="drift",bounds=bounds,control=control),silent=T)}
    if (model=="white") {fitted<-white<-try(fitContinuous(phy,x,model="white",bounds=bounds,control=control),silent=T)}
  }  
  M <- dist.nodes(phy) #uses GLS on [rescaled] phylo
  nb.tip <- length(phy$tip.label)
  varAY <- M[-(1:nb.tip), 1:nb.tip]
  varY <- M[1:nb.tip,1:nb.tip]      
  J<-array(1,dim=nb.tip)
  if(try(is.matrix(solve(varY)),silent=T)==TRUE){
    GrandMean<-J%*%solve(varY)%*%x / J%*%solve(varY)%*%J
    node.est<- varAY%*%solve(varY)%*%(x-GrandMean) + GrandMean[1,1]
  } else {
    warning("In node.est(): singular matrix: using dist matrix without the rescale, revert to BM")
    M <- dist.nodes(treedata.obj$phy)
    varAY <- M[-(1:nb.tip), 1:nb.tip]
    varY <- M[1:nb.tip,1:nb.tip]
    GrandMean<-J%*%solve(varY)%*%x / J%*%solve(varY)%*%J
    node.est<- varAY%*%solve(varY)%*%(x-GrandMean) + GrandMean[1,1]
    phy <- treedata.obj$phy
  }
  if(plot.est && !is.na(node.est)) {
    plot(dist.nodes(treedata.obj$phy)[,treedata.obj$phy$edge[1,1]],c(treedata.obj$data[,traitnum],node.est),xlab="Time",ylab="Trait",type="n")
    for(i in 1:length(treedata.obj$phy$edge[,1])){
      lines(dist.nodes(treedata.obj$phy)[,treedata.obj$phy$edge[1,1]][treedata.obj$phy$edge[i,]],c(treedata.obj$data[,traitnum],node.est)[treedata.obj$phy$edge[i,]])
    }
  }
  if(was.estimated){
    return (list(model=model,est=node.est,phy=phy,fitted=fitted))}
  else {return (list(model=model,est=node.est,phy=phy,fitted=fitted$opt))}
}  

#' @title node.estimate.fossils
#' @description To estimate nodes with the placement of fossils on randomly assigned or specifed edges on a tree.
#' @usage node.estimate.fossils()
#' @param treedata_min tree data object with min estimate of the climate envelope
#' @param treedata_max tree data object with max estimate of the climate envelope
#' @param fossils the estimate of the climate envelope of the fossil occurrences
#' @param fossils.edges the edge number that the fossil occurs on
# @examples
#ex_min <- treedata(ex_mytree[[1]],trialest$cem[,c(1,2)])
#colnames(ex_min$data) <- c("bio1","bio4")
#ex_max <- treedata(ex_mytree[[1]],trialest$cem[,c(7,8)])
#colnames(ex_max$data) <- c("bio1","bio4")
#fossils <- getBioclimVars(manipulatedFosssils, which.biovars=c(1,4,15))
#fossils.edges<-fossilsedges
#model <- "BM"
#bounds <- list()
#node.estimate.fossils(ex_min,ex_max,fossils=biovarFossils,fossils.edges=F)

node.estimate.fossils<-function(treedata_min,treedata_max,fossils=FALSE,fossils.edges=FALSE,model="BM", bounds=list(),control=list(),...){
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
        treedata_min <- treedata(addfossil(treedata_min$phy, mintime = fossils[i,1], maxtime = fossils[i,1] + 1, name = rownames(fossils)[i],edge = f.edge), rbind(treedata_min$data,fossils[i,colnames(treedata_min$data),drop=F] - addRangeFossils),sort=TRUE,warnings=F)
        rownames(treedata_min$data)<-treedata_min$phy$tip.label
        colnames(treedata_min$data)<-c.names
        treedata_max <- treedata(addfossil(treedata_max$phy, mintime = fossils[i,1], maxtime = fossils[i,1] + 1, name = rownames(fossils)[i],edge = f.edge), rbind(treedata_max$data,fossils[i,colnames(treedata_max$data),drop=F] + addRangeFossils),sort=TRUE,warnings=F)
        rownames(treedata_max$data)<-treedata_max$phy$tip.label
        colnames(treedata_max$data)<-c.names
      }
  }
  treedata_min$phy$node.label<-treedata_max$phy$node.label<-((length(treedata_min$phy$tip)+1):((length(treedata_min$phy$tip)*2)-1))
  #Here is the trait loop start
  for(trai in 1:num_traits){ 
    min_model[[trai]]<-node.estimate(treedata_min,trai,model=model,bounds=bounds,control=control,plot.est=FALSE)
    tmin<-min_model[[trai]]$est
    max_model[[trai]]<-node.estimate(treedata_max,trai,model=model,bounds=bounds,control=control,plot.est=FALSE)
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

#' @title getTimeSlice
#' @description This function extracts estimated ancestral reconstructions for continuous characters any time specified along a phylogeny for all lineages present at the specified time.
#' @usage getTimeSlice(time, phy, trait, plot.est = TRUE)
#' @param timeSlice single numeric or a vector with the time (or times) to extract the estimated ancestor reconstructions.
#' @param tree an object of the class "phylo" that should be dated
#' @param trait a vector of both tip values and node estimates that correspond to tree
#' @param plot.est a conditional stating whether or not to plot the results
#' @details The estimated reconstruction relies on an interpolation between node or between tip and node estimates of the trait. This method assumes a constant rate of evolution along the lineage where the interpolation is taking place.
#' @return edge for each time specified, a vector of edges that are present during that time are returned
#' @return est for each time specified, a vector of estimates of the ancestral reconstruction along each edge
#' @author A. Michelle Lawing
# @examples
#add example of how to get treedata_min
#ex_est <- node.estimate(treedata_min ,traitnum = 1)
#ex_timeSlice <- getTimeSlice(10,treedata_min$phy,c(treedata_min$data[,1],ex_est$est))

getTimeSlice<-function(timeSlice, tree, trait, model = "BM", plot.est=FALSE){
  M<-dist.nodes(tree)
  treeage<-max(M)/2
  maxedge<-(as.numeric(treeage - M[tree$edge[,1],tree$edge[1,1]]))
  minedge<-(as.numeric(treeage - M[tree$edge[,2],tree$edge[1,1]]))
  edgesample<-list()
  for(i in 1:length(timeSlice)){
    edgesample[[i]]<-which(maxedge>=timeSlice[i] & minedge<=timeSlice[i]+0.001)
  }
  hldtime<-list() 
  for(i in 1:length(timeSlice)){
    hld<-array(NA, dim=length(edgesample[[i]]))
    for(j in 1:length(edgesample[[i]])){
      percentchange<-(maxedge[edgesample[[i]][j]]-timeSlice[i])/(maxedge[edgesample[[i]][j]]-minedge[edgesample[[i]][j]])
      hld[j]<-trait[tree$edge[edgesample[[i]][j],1]] - percentchange*(trait[tree$edge[edgesample[[i]][j],1]]-trait[tree$edge[edgesample[[i]][j],2]])
    }
    hldtime[[i]]<-hld
  }
  if(plot.est){
    xx<-node.estimate(treedata(tree, trait[1:(tree$edge[1,1]-1)]),1, model = model, plot.est = plot.est)
    for(i in 1:length(timeSlice)){
      points(array(max(M)/2-timeSlice[i],dim=length(hldtime[[i]])),hldtime[[i]])
    }
  }
  return(list(edge=edgesample,est=hldtime))
}

#' @title getEnvelopes
#' @description This function gets the bioclimate envelopes of species and nodes.
#' @usage getEnvelopes(treedata_min,treedata_max,node_est)

getEnvelopes<-function(treedata_min,treedata_max,node_est){
  require(ape)
  num_traits<-length(treedata_min$data[1,])
  num_species<-length(treedata_min$data[,1])
  envelope<-array(NA, dim=c(2*num_species-1,5,num_traits))
  for(i in 1:num_traits){
    M<-dist.nodes(treedata_min$phy)[1,num_species+1]-dist.nodes(treedata_min$phy)[,num_species+1]
    temp <- array(unlist(node_est),dim=c(2,length(node_est[[1]][1,,1]),num_traits,length(node_est)))
    traitgram_min_min<-cbind(c(treedata_min$data[,i],sapply(1:(num_species-1),function(x) min(temp[1,x,i,]))),M)
    traitgram_min_max<-cbind(c(treedata_min$data[,i],sapply(1:(num_species-1),function(x) max(temp[1,x,i,]))),M)
    traitgram_max_min<-cbind(c(treedata_max$data[,i],sapply(1:(num_species-1),function(x) min(temp[2,x,i,]))),M)
    traitgram_max_max<-cbind(c(treedata_max$data[,i],sapply(1:(num_species-1),function(x) max(temp[2,x,i,]))),M)
    traitgram_min<-cbind(c(treedata_min$data[,i],colMeans(t(temp[1,1:(num_species-1),i,]))),M)
    traitgram_max<-cbind(c(treedata_max$data[,i],colMeans(t(temp[2,1:(num_species-1),i,]))),M)
    envelope[,,i]<-cbind(traitgram_min_min[,c(2,1)],traitgram_min_max[,1],traitgram_max_min[,1],traitgram_max_max[,1])
  }
  return(envelope)
}

#' @title plotTraitGram
#' @description Combine the node estimates based on random or specified fossil placement and plot them on a phylotrait gram in a specified directory.
#' @usage plotTraitGram()
#' @param treedata_min a tree data object with the min estimate of the climate envelope
#' @param treedata_max a tree data object with the max estimate of the climate envelope
#' @param node_est the estimate of all the nodes, both min and max

plotTraitGram<-function(treedata_min,treedata_max,node_est,fossils=FALSE,which.biovars,path){
  require(ape)
  num_traits<-length(treedata_min$data[1,])
  num_species<-length(treedata_min$data[,1])
  data(paleoclimate) #isotopically scaled paleoclimate bioclimate variables for North America
  for(i in 1:num_traits){
    M<-dist.nodes(treedata_min$phy)[1,num_species+1]-dist.nodes(treedata_min$phy)[,num_species+1]
    temp <- array(unlist(node_est),dim=c(2,length(node_est[[1]][1,,1]),num_traits,length(node_est)))
    traitgram_min_min<-cbind(c(treedata_min$data[,i],sapply(1:(num_species-1),function(x) min(temp[1,x,i,]))),M)
    traitgram_min_max<-cbind(c(treedata_min$data[,i],sapply(1:(num_species-1),function(x) max(temp[1,x,i,]))),M)
    traitgram_max_min<-cbind(c(treedata_max$data[,i],sapply(1:(num_species-1),function(x) min(temp[2,x,i,]))),M)
    traitgram_max_max<-cbind(c(treedata_max$data[,i],sapply(1:(num_species-1),function(x) max(temp[2,x,i,]))),M)
    traitgram_min<-cbind(c(treedata_min$data[,i],colMeans(t(temp[1,1:(num_species-1),i,]))),M)
    traitgram_max<-cbind(c(treedata_max$data[,i],colMeans(t(temp[2,1:(num_species-1),i,]))),M)
    #check on the size of graph and size of text
    pdf(paste(path,"bio",i,".pdf",sep=""),width=960,height=960,pointsize=20,useDingbats=F)
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

#' @title plotTraitGramMultiPhylo
#' @description Combine the node estimates based on random or specified fossil placement and plot them on a phylotrait gram in a specified directory.
#' @usage plotTraitGramMultiPhylo()
#' @param treedata_min tree data object with min estimate of the climate envelope
#' @param treedata_max tree data object with max estimate of the climate envelope

plotTraitGramMultiPhylo<-function(treedata_min,treedata_max,node_est,fossils=FALSE,which.biovars,path,alpha.trans=as.integer(255/(1+log(length(node_est))))){
  require(ape)
  num_traits<-length(treedata_min[[1]]$data[1,])
  num_species<-length(treedata_min[[1]]$data[,1])
  data(paleoclimate) #isotopically scaled paleoclimate bioclimate variables for North America
  traitgram_min_min<-as.list(array(NA,dim=length(node_est)))
  traitgram_min_max<-as.list(array(NA,dim=length(node_est)))
  traitgram_max_min<-as.list(array(NA,dim=length(node_est)))
  traitgram_max_max<-as.list(array(NA,dim=length(node_est)))
  traitgram_min<-as.list(array(NA,dim=length(node_est)))
  traitgram_max<-as.list(array(NA,dim=length(node_est)))
  for(i in 1:num_traits){
    for(j in 1:length(node_est)){
      M<-dist.nodes(treedata_min[[j]]$phy)[1,num_species+1]-dist.nodes(treedata_min[[j]]$phy)[,num_species+1]
      temp <- array(unlist(node_est[[j]]),dim=c(2,length(node_est[[j]][[1]][1,,1]),num_traits,length(node_est[[j]])))
      traitgram_min_min[[j]]<-cbind(c(treedata_min[[j]]$data[,i],sapply(1:(num_species-1),function(x) min(temp[1,x,i,]))),M)
      traitgram_min_max[[j]]<-cbind(c(treedata_min[[j]]$data[,i],sapply(1:(num_species-1),function(x) max(temp[1,x,i,]))),M)
      traitgram_max_min[[j]]<-cbind(c(treedata_max[[j]]$data[,i],sapply(1:(num_species-1),function(x) min(temp[2,x,i,]))),M)
      traitgram_max_max[[j]]<-cbind(c(treedata_max[[j]]$data[,i],sapply(1:(num_species-1),function(x) max(temp[2,x,i,]))),M)
      traitgram_min[[j]]<-cbind(c(treedata_min[[j]]$data[,i],colMeans(t(temp[1,1:(num_species-1),i,]))),M)
      traitgram_max[[j]]<-cbind(c(treedata_max[[j]]$data[,i],colMeans(t(temp[2,1:(num_species-1),i,]))),M)
    }
    #check on the size of graph and size of text
    pdf(paste(path,colnames(treedata_min[[1]]$data)[i],".pdf",sep=""),width=960,height=960,pointsize=20,useDingbats=F)
    plot(traitgram_min[[1]],ylim=rev(range(c(traitgram_min_min[[1]][,2]))),xlim=range(c(paleoclimate[[1]][,which.biovars[i]+3], paleoclimate[[11]][,which.biovars[i]+3], paleoclimate[[16]][,which.biovars[i]+3], traitgram_min_min,traitgram_max_max)), type="n",xlab=colnames(treedata_min[[1]]$data)[i],ylab="Time (mya)")    
    lines(c(min(paleoclimate[[1]][,which.biovars[i]+3],traitgram_min[[j]]),max(paleoclimate[[1]][,which.biovars[i]+3],traitgram_max[[j]])),
          c(0,0),col="antiquewhite",lwd=10)
    lines(c(min(paleoclimate[[6]][,which.biovars[i]+3]),max(paleoclimate[[6]][,which.biovars[i]+3])),c(5,5),col="antiquewhite",lwd=10)
    lines(c(min(paleoclimate[[11]][,which.biovars[i]+3]),max(paleoclimate[[11]][,which.biovars[i]+3])),c(10,10),col="antiquewhite",lwd=10)
    lines(c(min(paleoclimate[[16]][,which.biovars[i]+3]),max(paleoclimate[[16]][,which.biovars[i]+3])),c(15,15),col="antiquewhite",lwd=10)        
    for(j in 1:length(node_est)){
      for(k in 1:104) {lines(c(traitgram_max_min[[j]][k,1],traitgram_max_max[[j]][k,1]),c(traitgram_max_min[[j]][k,2],traitgram_max_min[[j]][k,2]),col=rgb(135,206,235,alpha=alpha.trans,maxColorValue=255),lwd=6)}
      for(k in 1:104) {lines(c(traitgram_min_min[[j]][k,1],traitgram_min_max[[j]][k,1]),c(traitgram_min_min[[j]][k,2],traitgram_min_min[[j]][k,2]),col=rgb(128,128,128,alpha=alpha.trans,maxColorValue=255),lwd=4)}
      for(k in 1:104) {lines(traitgram_max[[j]][treedata_max[[j]]$phy$edge[k,],1],traitgram_max[[j]][treedata_max[[j]]$phy$edge[k,],2],col=rgb(135,206,235,alpha=alpha.trans,maxColorValue=255))}
      for(k in 1:104) {lines(traitgram_min[[j]][treedata_min[[j]]$phy$edge[k,],1],traitgram_min[[j]][treedata_min[[j]]$phy$edge[k,],2],col=rgb(128,128,128,alpha=alpha.trans,maxColorValue=255))}
    }
    if(length(fossils)!=1){for(k in 1:length(fossils[,1])) {points(fossils[k,i+3],fossils[k,1],col="black",pch=16)}}
    dev.off()
  }
}

#' @title plotAnimatedPPGM
#' @description This function creates an animated gif showing the change in modeled suitable habitat through time in geographic space. It require ImageMagick or GraphicsMagick to be previously installed in the operating system. This is easy to do if you have macports. Just type sudo port install ImageMagick into terminal.
#' @usage plotAnimatedPPGM()
#' @param envelope the min and max envelope of each lineage for each time slice
#' @param tree the phylogeny or multiple phylogenies that show the relationship between species

plotAnimatedPPGM<-function(envelope,tree,filename="ppgm.gif",which.biovars,path=""){
  require(animation)
  data(paleoclimate) #isotopically scaled paleoclimate bioclimate variables for North America
  temp_min<-lapply(1:length(paleoclimate),function(i){
    temp<-lapply(1:length(which.biovars),function(j){getTimeSlice(i-1,tree,envelope[,2,j])})
    temp<-t(array(unlist(temp),dim=c(length(unlist(temp[[1]]$edge)),2*length(which.biovars))))
    return(temp)})
  temp_max<-lapply(1:length(paleoclimate),function(i){
    temp<-lapply(1:length(which.biovars),function(j){getTimeSlice(i-1,tree,envelope[,5,j])})
    temp<-t(array(unlist(temp),dim=c(length(unlist(temp[[1]]$edge)),2*length(which.biovars))))
    return(temp)})
  richnesscount<-lapply(1:length(paleoclimate), function(j){
    hld<-array(0,dim=length(paleoclimate[[j]][,1]))
    for(i in 1:length(temp_min[[j]][1,])){
      matching<-
        sapply(1:length(which.biovars),function(x){paleoclimate[[j]][,which.biovars[x]+3]>temp_min[[j]][1:length(which.biovars)*2,i][x] & paleoclimate[[j]][
          ,which.biovars[x]+3]<temp_max[[j]][1:length(which.biovars)*2,i][x]})
      matching<-which(rowSums(matching)==length(which.biovars),arr.ind=TRUE)    
      hld[matching]<-hld[matching]+1
    }
    hld[which(hld==0,arr.ind=TRUE)]=NA
    return(hld)
  })
  saveGIF(for(i in 1:length(paleoclimate)){
    plot(paleoclimate[[i]][,2:3],cex=0.5,pch=16,col="lightgray",xlim=c(-180,0),ylim=c(0,90))
    points(paleoclimate[[i]][,2:3],cex=0.5,pch=16,col=colorRampPalette(c("#'FFE5CC", "#'FF8000", "#'990000"))(length(temp_min[[1]][1,]))[richnesscount[[i]]],xlim=c(-180,0),ylim=c(0,90))
  },movie.name=filename,outdir=getwd()) 
}

#' @title plotAnimatedPPGMMultiPhylo
#' @description This function creates an animated gif showing the change in modeled suitable habitat through time in geographic space. It requires ImageMagick or GraphicsMagick to be previously installed in the operating system. This is easy to do if you have macports. Just type sudo port install ImageMagick into terminal.
#' @usage plotAnimatedPPGMMultiPhylo()
#' @param envelope the min and max envelope of each lineage for each time slice
#' @param tree the phylogeny or multiple phylogenies that show the relationship between species

#make multiphylo and make path work
plotAnimatedPPGMMultiPhylo<-function(envelope,tree,filename="ppgm.gif",which.biovars,path=""){
  require(animation)
  require(Hmisc)
  if(path==""){out=getwd()}
  if(path!=""){out=paste(getwd(),"/",substr(path,1,nchar(path)-1),sep="")}
  data(paleoclimate) #isotopically scaled paleoclimate bioclimate variables for North America
  temp_min<-as.list(array(NA,dim=length(tree)))
  temp_max<-as.list(array(NA,dim=length(tree)))
  richnesscount<-as.list(array(NA,dim=length(tree)))
  for(tr in 1:length(tree)){
    temp_min[[tr]]<-lapply(1:length(paleoclimate),function(i){
    temp<-lapply(1:length(which.biovars),function(j){getTimeSlice(i-1,tree[[tr]],envelope[[tr]][,2,j])})
    temp<-t(array(unlist(temp),dim=c(length(unlist(temp[[1]]$edge)),2*length(which.biovars))))
    return(temp)})
    temp_max[[tr]]<-lapply(1:length(paleoclimate),function(i){
    temp<-lapply(1:length(which.biovars),function(j){getTimeSlice(i-1,tree[[tr]],envelope[[tr]][,5,j])})
    temp<-t(array(unlist(temp),dim=c(length(unlist(temp[[1]]$edge)),2*length(which.biovars))))
    return(temp)})
    richnesscount[[tr]]<-lapply(1:length(paleoclimate), function(j){
    hld<-array(0,dim=length(paleoclimate[[j]][,1]))
    for(i in 1:length(temp_min[[tr]][[j]][1,])){
      matching<-
        sapply(1:length(which.biovars),function(x){paleoclimate[[j]][,which.biovars[x]+3]>temp_min[[tr]][[j]][1:length(which.biovars)*2,i][x] & paleoclimate[[j]][
          ,which.biovars[x]+3]<temp_max[[tr]][[j]][1:length(which.biovars)*2,i][x]})
      matching<-which(rowSums(matching)==length(which.biovars),arr.ind=TRUE)    
      hld[matching]<-hld[matching]+1
    }
    hld[which(hld==0,arr.ind=TRUE)]=NA
    return(hld)
  })
  }
  saveGIF(for(i in 1:length(paleoclimate)){
    richnesscountMean<-rowMeans(array(unlist(lapply(1:length(tree),function(tr) richnesscount[[tr]][[i]])),dim=c(length(richnesscount[[1]][[i]]),length(tree))),na.rm=TRUE)
    lq<-array(unlist(lapply(1:length(tree),function(tr) richnesscount[[tr]][[i]])),dim=c(length(richnesscount[[1]][[i]]),length(tree)))
    richnesscountMIN<-suppressWarnings(unlist(lapply(1:length(lq[,1]),function(r) min(lq[r,],na.rm=TRUE))))
    richnesscountMAX<-suppressWarnings(unlist(lapply(1:length(lq[,1]),function(r) max(lq[r,],na.rm=TRUE))))                         
    par(mar=c(3,0,3,1))
    plot(paleoclimate[[i]][,2:3],cex=0.5,xlab="",ylab="",axes=FALSE,pch=16,col="lightgray",xlim=c(-200,0),ylim=c(0,90))
    points(paleoclimate[[i]][,2:3],cex=0.5,pch=16,col=colorRampPalette(c("#'FFE5CC", "#'FF8000", "#'990000"))(length(temp_min[[tr]][[1]][1,]))[richnesscountMean],xlim=c(-200,0),ylim=c(0,90))
    center<-colMeans(paleoclimate[[i]][,2:3])
    scale<-(cbind(paleoclimate[[i]][,2]-center[1],paleoclimate[[i]][,3]-center[2]))/3
    translate<-cbind(scale[,1]-160,scale[,2]+25)
    points(translate,cex=0.5,pch=16,col="lightgray")
    points(translate,cex=0.5,pch=16,col=colorRampPalette(c("#'FFE5CC", "#'FF8000", "#'990000"))(length(temp_min[[tr]][[1]][1,]))[richnesscountMIN])
    translate<-cbind(scale[,1]-35,scale[,2]+25)
    points(translate,cex=0.5,pch=16,col="lightgray")
    points(translate,cex=0.5,pch=16,col=colorRampPalette(c("#'FFE5CC", "#'FF8000", "#'990000"))(length(temp_min[[tr]][[1]][1,]))[richnesscountMAX])
    text(-160,40,"Min")
    text(-35,40,"Max")
    text(-90,85,"Mean")
    text(-90,90,"Modeled Richness over a distribution of trees")
    },movie.name=filename,outdir=out) 
}

#' @title getGeoRate
#' @description This function calculates the change in suitable habitat through time in geographic space.
#' @usage getGeoRate()
#' @param envelope the min and max envelope of each lineage for each time slice
#' @param tree the phylogeny or multiple phylogenies that show the relationship between species

getGeoRate<-function(envelope,tree,which.biovars){
  require(phangorn)
  require(fields)
  data(paleoclimate) #isotopically scaled paleoclimate bioclimate variables for North America
  temp_min<-lapply(1:length(paleoclimate),function(i){
    temp<-lapply(1:length(which.biovars),function(j){getTimeSlice(i-1,tree,envelope[,2,j])})
    temp<-t(array(unlist(temp),dim=c(length(unlist(temp[[1]]$edge)),2*length(which.biovars))))
    return(temp)})
  temp_max<-lapply(1:length(paleoclimate),function(i){
    temp<-lapply(1:length(which.biovars),function(j){getTimeSlice(i-1,tree,envelope[,5,j])})
    temp<-t(array(unlist(temp),dim=c(length(unlist(temp[[1]]$edge)),2*length(which.biovars))))
    return(temp)})
  lineage<-lapply(1:length(temp_min[[1]][1,]),function(i) {c(i,Ancestors(tree,i))})
  geo_center<-array(NA, dim=c(length(paleoclimate),length(lineage),2))
  geo_size<-array(NA, dim=c(length(paleoclimate),length(lineage)))
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
  #calculate the geo dist from centers and expansions/contractions from areas
  C=array(NA, dim=c(length(lineage),length(paleoclimate),length(paleoclimate)))
  A=array(NA, dim=c(length(lineage),length(paleoclimate),length(paleoclimate)))
  for(i in 1:length(lineage)){
    C[i,,]=rdist.earth(geo_center[,i,],geo_center[,i,],miles=FALSE)
    A[i,,]=sapply(1:length(paleoclimate),function(x){geo_size[x,i]-geo_size[,i]})
  }
  T=as.matrix(dist(1:length(paleoclimate),diag=TRUE))
  return(list(geo_center=C, geo_size=A, time_int=T))
}

#' @title plotBumpChart
#' @usage plotBumpChart()
#' @return bump chart of mean species-climate associations

plotBumpChart<-function(sp_data_mean, path = ""){
  require(ggplot2)
  #require(gridExtra) #' package does not exist anymore... check out functionality...
  require(reshape2)
  jpeg(paste(path,"bump_chart.jpg",sep=""),width=960,height=960,quality=100,pointsize=20)
  trans <- apply(sp_data_mean,2,rank)
  melted <- melt(trans)
  which.biovars<-as.numeric(colnames(sp_data_mean))
  theme_set(theme_bw())
  zmargin<-theme(panel.margin=unit(0,"lines"))
  b1 <- ggplot(melted, aes(X2, value, group = X1, color = X1, label = X1)) + geom_line() + 
    geom_text(aes(label = as.factor(c(paste("S.",as.factor(unique(melted$X1))),rep("",length(melted$X2) - length(as.character(unique(melted$X1)))))), 
                x = melted$X2 - 0.01, size = 3, hjust = 1, fontface = "italic")) + theme_bw() +
    theme(legend.position = "none", panel.border = element_blank(), axis.ticks = element_blank()) +
    scale_x_discrete(breaks = levels(as.factor(melted$X2)), labels = paste("Bio",which.biovars, sep = ""), expand = c(0.5,0.9)) +
    xlab(NULL) + scale_y_continuous(breaks = NULL) + ylab(NULL) +
    annotate("text", x = which.biovars, y = rep(max(trans) + 2, length(which.biovars)),
           label = as.integer(apply(sp_data_mean,2,max)), size = 3) +
    annotate("text", x = which.biovars, y = rep(min(trans) - 2, length(which.biovars)),
           label = as.integer(apply(sp_data_mean,2,min)), size = 3)
  b1
  dev.off()
}

#' @title plotGeoRates
#' @usage plotGeoRates()
#' @return Gray background of all pairwise comparisons of change in geo center and area through time. Blue points on top show the sequential change in geo center and expansion/contraction for all lineages.

plotGeoRates<-function(geo_center,geo_size,time_int,trees,path=""){
  jpeg(paste(path,"geo_rates.jpg",sep=""),width=960,height=960,quality=100,pointsize=20)
  par(mfrow=c(2,1),mar=c(4,4,1,1),mgp=c(2.25,1,0))
  if(is.finite(max(geo_center))){ymax=max(geo_center)}
  if(!is.finite(max(geo_center))){ymax=4000}
  plot(1,1,type="l",xlim=c(0,20),ylim=c(0,ymax),xlab="",ylab="Change in Geo Center")
  for(t in 1:length(trees)){
    for(l in 1:length(geo_center[t,,1,1])){
      points(time_int[,],geo_center[t,l,,],col="lightgray")
    }
  }
  for(t in 1:length(trees)){
    for(l in 1:length(geo_center[1,,1,1])){
      for(i in 1:20){
        points(time_int[1,(i+1)],geo_center[t,l,i,(i+1)],col="blue")
      }
    }
  }
  plot(1,1,type="l",xlim=c(0,20),ylim=c(min(geo_size),max(geo_size)),xlab="Time",ylab="Change in Geo Area")
  for(t in 1:length(trees)){
    for(l in 1:length(geo_size[t,,1,1])){
      points(time_int[,],geo_size[t,l,,],col="lightgray")
    }
  }
  for(t in 1:length(trees)){
    for(l in 1:length(geo_size[1,,1,1])){
      for(i in 1:20){
        points(time_int[1,(i+1)],geo_size[t,l,i,(i+1)],col="blue")
      }
    }
  }
  #reset par
  par(mfrow=c(1,1),mar=c(5,4,4,2),mgp=c(3,1,0))
  dev.off()
}

#' @title ppgm
#' @description ppgm makes a paleophylogeographic species distribution model using the bioclimate envelope method for a specified time period. Currently, models are only available for North America.
#' @usage ppgm(occurrences, fossils, tree, fossils.edges, model = "BM", plot.est=TRUE)
#' @param occurrences a matrix with three columns of species name, longitude, and latitude, in that order, and rows that are entries for species occurrences. The bioclimate variables can be included for each occurrence in following columns. They must be in order 1 through 19.
#' @param fossils a matrix with three columns of age to the closest million year integer, longitude, and latitude, in that order, and rows that are entries for fossil occurrences. The bioclimate variables can be included for each occurrence in following columns. They must be in order 1 through 19.
#' @param tree phylogeny of species from first column of occurrences argument.
#' @param fossils.edges a vector of edges that the fossils belong to. Must be in the same order of the fossils argument. If fossils.edges is false, the the function randomly assigns the location of the fossils depending on the age (see details for more information).
#' @param model the model of evolution to use to estimate ancestor nodes. Argument is passed onto to function node.estimate.
#' @param permut the number of times to randomly place fossils in phylogeny and estimate ancestor states.
#' @param which.biovars a vector with the biovars to include in model (see www.worldclim.org for a list of biovars). If "ALL", then all 19 biovars are included in analysis.
#' @param path path to the directory where the results should be saved.
#' @details If the 19 bioclimate variables are not supplied with the occurrences or with the fossils, they will be extracted from the closest 50km point location in the modern or paleoclimate maps that are loaded in with this function. The paleoclimate maps are isotopically scaled between general circulation models (see Lawing and Polly 2011; Rodder et al. 2013) and modern climate (see Hijmans et al. 2005). The fossils paleoclimate data is extracted to the closest million year paleoclimate map. Paleoclimate maps are derived at one million year intervals for the past 20 Ma. The tree (phylogeny) should be dichotomous and the species names should match the names in the first column of the occurrences argument. fix to work with the nice plots on multiple trees, and check and see if it works when no suitable habitat is identified. also, make option to print any animated lineage through time?
# @examples

ppgm<-function(occurrences, fossils = FALSE, trees, fossils.edges = FALSE, model = "BM", permut = 2, only.biovars = TRUE, 
               which.biovars = c(1,12), path = "", plot.TraitGram = TRUE, plot.AnimatedMaps = TRUE, plot.GeoRates = TRUE, 
               plot.BumpChart = FALSE, bounds = list(), control = list(), use.paleoclimate = TRUE, verbose = TRUE){
  #calculate the alpha.trans, which is the transparancy for all the trees plotted on top of each other
  alpha.trans <- as.integer(255 / (1 + log(length(trees))))
  #assign rownames to fossils
  if(length(fossils)!=1){rownames(fossils)<-paste("fossil",1:length(fossils[,1]),sep="")}
  #Extract bioclim variables for species occurrences, if not supplied by user
  if(only.biovars) {
    occurrences<-getBioclimVars(occurrences,which.biovars)
    if(length(fossils)!=1){
      fossils<-getBioclimVars(fossils,which.biovars)
    }
  }
  #check the occurrences variable conforms with requirments
  if(length(occurrences[1, ]) < 4) {
    return ("ERROR: There are less than four columns for the occurrences input. Check that occurrences input conforms to requirements.")
  }
  #load paleocliamte data
  if(is.logical(use.paleoclimate)) {
    data(paleoclimate)
  } else {
    paleoclimate <- use.paleoclimate
  }
  geo_center<-array(NA,dim=c(length(trees),length(unique(occurrences$Species)),length(paleoclimate),length(paleoclimate)))
  geo_size<-array(NA,dim=c(length(trees),length(unique(occurrences$Species)),length(paleoclimate),length(paleoclimate)))
  time_int<-array(NA,dim=c(length(paleoclimate),length(paleoclimate)))
  treedata_min<-as.list(array(NA,dim=length(trees)))
  treedata_max<-as.list(array(NA,dim=length(trees)))
  node_est<-as.list(array(NA,dim=length(trees)))
  model_min<-as.list(array(NA,dim=length(trees)))
  model_max<-as.list(array(NA,dim=length(trees)))
  envelope<-as.list(array(NA,dim=length(trees)))
  #get bioclimate envelopes
  sp_data_min<-sapply(4:(length(which.biovars)+3),function(x) tapply(occurrences[,x],occurrences$Species,min))
  sp_data_mean<-sapply(4:(length(which.biovars)+3),function(x) tapply(occurrences[,x],occurrences$Species,mean))
  sp_data_max<-sapply(4:(length(which.biovars)+3),function(x) tapply(occurrences[,x],occurrences$Species,max))
  colnames(sp_data_mean)<-which.biovars
  if(plot.BumpChart){
    plotBumpChart(sp_data_mean,path)  #broken?? Error: breaks and labels have unequal lengths
  }
  for(tr in 1:length(trees)){
    #Check if phylogeny is dichotomous, if not, make it dichotomous
    if(!is.binary.tree(trees[[tr]])){trees[[tr]]<-multi2di(trees[[tr]])}
    #make treedata object for bioclimate envelopes and phylogeny
    treedata_min[[tr]]<-treedata(trees[[tr]],sp_data_min,sort=TRUE,warnings=F)
    treedata_max[[tr]]<-treedata(trees[[tr]],sp_data_max,sort=TRUE,warnings=F) 
    colnames(treedata_min[[tr]]$data)<-colnames(treedata_max[[tr]]$data)<-paste("bio",which.biovars,sep="")
    #to estimate nodes, place fossils randomly or as specified on edges from fossils.edges argument    
    full_est <- list()
    for(pr in 1:permut){
      full_est[[pr]] <- node.estimate.fossils(treedata_min=treedata_min[[tr]],treedata_max=treedata_max[[tr]],fossils=fossils,fossils.edges=fossils.edges,model=model,bounds=bounds,control=control)
    }      
#####################
    node_est[[tr]]<-lapply(1:permut, function(p) full_est[[p]]$est)
    model_min[[tr]]<-lapply(1:permut, function(p) full_est[[p]]$min_model)
    model_max[[tr]]<-lapply(1:permut, function(p) full_est[[p]]$max_model)
    #get bioclimate envelopes for species and nodes
    envelope[[tr]]<-getEnvelopes(treedata_min[[tr]],treedata_max[[tr]],node_est[[tr]])
    #get data from geo displacement
    temp<-getGeoRate(envelope[[tr]],tree=trees[[tr]],which.biovars=which.biovars)
    geo_center[tr,,,]<-temp$geo_center
    geo_size[tr,,,]<-temp$geo_size
  }
  time_int<-temp$time_int
  #plot permutations
  if(plot.TraitGram){
    plotTraitGramMultiPhylo(treedata_min,treedata_max,node_est,fossils=fossils,which.biovars=which.biovars,path=path,alpha.trans=alpha.trans)
  }
  #plot animated maps
  if(plot.AnimatedMaps){
    plotAnimatedPPGMMultiPhylo(envelope,tree=trees,filename="ppgm.gif",which.biovars=which.biovars,path=path)
  }
  #plot geo rates
  if(plot.GeoRates){
    plotGeoRates(geo_center,geo_size,time_int,trees,path=path)
  }
  lineage_geo_center=array(NA,dim=c(length(trees),length(geo_center[1,,1,1])))
  for(tr in 1:length(trees)){
    for(l in 1:length(geo_center[tr,,1,1])){
      #rate in km/year
      #need to consider rates and time scale issue
      lineage_geo_center[tr,l]=(mean((geo_center[tr,l,,]/time_int)[lower.tri(geo_center[tr,l,,]/time_int)],na.rm=T))
    }
  }
  lineage_geo_size=array(NA,dim=c(length(trees),length(geo_size[1,,1,1])))
  for(tr in 1:length(trees)){
    for(l in 1:length(geo_size[tr,,1,1])){
      #rate in km/year
      lineage_geo_size[tr,l]=(mean((geo_size[tr,l,,]/time_int)[lower.tri(geo_size[tr,l,,]/time_int)],na.rm=T))
    }
  }    
  cem<-data.frame(sp_data_min,sp_data_mean,sp_data_max)
  names(cem)<-c(paste(colnames(sp_data_mean),"Min",sep=""),paste(colnames(sp_data_mean),"Mean",sep=""),paste(colnames(sp_data_mean),"Max",sep=""))
  geo_move<-data.frame(colMeans(lineage_geo_center,na.rm=T),colMeans(lineage_geo_size,na.rm=T))
  names(geo_move)<-c("RateGeoCenter","RateGeoSize")
  if(verbose){
    return(list(cem=cem,
              geo_move=geo_move,
              change_geo_center=geo_center, 
              change_geo_size=geo_size,
              time_int=time_int,
              treedata_min=treedata_min,
              treedata_max=treedata_max,
              model_min=model_min,
              model_max=model_max,
              node_est=node_est))    
  }
  else{
    return(list(cem=cem,
             geo_move=geo_move,
             change_geo_center=geo_center, 
             change_geo_size=geo_size,
             time_int=time_int))
  }
}


#' @title ppgmMESS
#' @description This creates a MESS map for given time slices, climate envelopes, and paleoclimate models.
#' @usage plotMESS()
#' @param cem_min the cem min output from the ppgm function. cbind() if there are multiple variables.
#' @param cem_max the cem max output from the ppgm function. cbind() if there are multiple variables.
#' @param est the node_est output from the ppgm function.
#' @param timeslice the time in million of years ago to project MESS maps (0 to 20). can handle single timeslice or vector of times.
#' @param tree the phylogeny or multiple phylogenies that show the relationship between species
#' @param which.biovars the biovariable number(s) between 1 and 19. handles vectors. 

#cem_min = cbind(trialest$cem[,1],trialest1$cem[, 1])
#cem_max = cbind(trialest$cem[,7],trialest1$cem[, 3])
#rownames(cem_min) <- rownames(cem_max) <- rownames(trialest$cem)
#extract one variable from the results of multivar run
#relist1 <- lapply(lapply(1:length(tree),function(x) array(unlist(trialest$node_est),dim=c(2,52,3,100))[,,1,x]),list)
#relist2 <- lapply(lapply(1:length(tree),function(x) array(unlist(trialest1$node_est),dim=c(2,52,1,100))[,,1,x]),list)
#est = list(relist1,relist2)
#mess <- ppgmMESS(cem_min, cem_max, est, tree = ex_mytree, timeslice = c(2,5,13,20), which.biovars = c(1,6))

ppgmMESS <- function(cem_min, cem_max, est, tree, fossils, timeslice, which.biovars, path = "", which.plot = "all"){
  
  #isotopically scaled paleoclimate bioclimate variables for North America
  data(paleoclimate)
  
  colorscheme <- colorRampPalette(c("blue","cyan","greenyellow","yellow","darkorange","red"))(250)
  
  MESS_score<-as.list(array(NA,dim=length(timeslice)))
  for(p in 1:length(timeslice)){
    MESS_score[[p]] <- array(NA,dim=c(length(paleoclimate[[(timeslice[p] + 1)]][,1]),length(which.biovars)))
  }
  
  #getenvelope
  for(b in 1:length(which.biovars)){
    treedata_min <- as.list(array(NA,dim = length(tree)))
    treedata_max <- as.list(array(NA,dim = length(tree)))
    yikes <- as.list(array(NA,dim = length(tree)))
    envelope<-as.list(array(NA,dim = length(tree)))
    for(tr in 1:length(tree)){
      if(length(which.biovars) > 1){
        data_min <- cem_min[,b]
        data_max <- cem_max[,b]
        names(data_min) <- names(data_max) <- rownames(cem_min)
      } else {
        data_min <- cem_min
        data_max <- cem_max
        names(data_min) <- names(data_max) <- names(cem_min)
      }
      treedata_min[[tr]] <- treedata(tree[[tr]], data = data_min, sort = TRUE, warnings = F)
      treedata_max[[tr]] <- treedata(tree[[tr]], data = data_max, sort = TRUE, warnings = F)
      yikes[[tr]] <- list(array(est[[b]][[tr]][[1]],dim=c(2,length(est[[b]][[tr]][[1]][1,]),1)))
      #trickery to get yikes in the correct format for getEnvelopes function. would be wise to rethink this.
    }
    for(tr in 1:length(tree)){
      envelope[[tr]] <- getEnvelopes(treedata_min[[tr]], treedata_max[[tr]], yikes[[tr]])
    }
    
    #get dims to alter data format
    env <- dim(envelope[[1]])
    env[3] <- env[3]*length(tree)
    branch_time <- round(envelope[[1]][,1,1], digits = 2)
    adj_data <- array(unlist(envelope), dim = env)
    adj_data <- adj_data[,-1,]
    adj_data <- array(c(adj_data),dim=c(env[1],(env[2]-1)*env[3]))
    
    #calculate the min and max variable that all the species and hyp nodes find suitable
    temp_min <- apply(adj_data, 1, min)
    temp_max <- apply(adj_data, 1, max)
    
    for(p in 1:length(timeslice)){
      spdata <- SpatialPoints(paleoclimate[[(timeslice[p] + 1)]][, 2:3])
      proj4string(spdata)  <- CRS("+init=epsg:4326")
      spdata <- spTransform(spdata, CRS("+init=epsg:26978"))
      if(sum(fossils[, 1] == (timeslice[p] + 1)) != 0){
        spfossils <- SpatialPoints(fossils[, 2:3])
        proj4string(spfossils)  <- CRS("+init=epsg:4326")
        spfossils <- spTransform(spfossils, CRS("+init=epsg:26978"))
        spfossils <- spfossils[fossils[, 1] == (timeslice[p] + 1), ]
      }
      getlineages <- which(branch_time < (timeslice[p] + 1) & branch_time >= ((timeslice[p] + 1) - 1))
      if(!length(getlineages) > 0) {getlineages <- which.max(branch_time)}
      reference_set <- c(adj_data[getlineages,])
      min_val <- min(temp_min[getlineages])
      max_val <- max(temp_max[getlineages])
      for(g in 1:length(paleoclimate[[(timeslice[p] + 1)]][,1])) {
        fi <- (sum(reference_set < paleoclimate[[(timeslice[p] + 1)]][g, (which.biovars[b] + 3)]) / length(reference_set) ) * 100          
        if(fi == 0){
          MESS_score[[p]][g, b] <- 100 * (paleoclimate[[(timeslice[p] + 1)]][g, (which.biovars[b] + 3)] - min_val)/(max_val - min_val)
        } else if (fi < 50) {
          MESS_score[[p]][g, b] <- 0
          #MESS_score[[p]][g, b] <- 2 * fi
        } else if (fi < 100) {
          MESS_score[[p]][g, b] <- 0
          #MESS_score[[p]][g, b] <- 2 * (100 - fi)
        } else if (fi == 100) {
          MESS_score[[p]][g, b] <- (100 * (max_val - paleoclimate[[(timeslice[p] + 1)]][g, (which.biovars[b] + 3)]))/(max_val - min_val)
        }
      }
      MESS_score[[p]] <- (apply(MESS_score[[p]], 2, scale, center = F) * 100) - 1
      MESS_score[[p]][MESS_score[[p]][,b] < -250] <- -250
      MESS_score[[p]][MESS_score[[p]][,b] > 0] <- 0
      #check on the size of graph and size of text
      if(which.plot == "all"){
        pdf(paste("MESS",timeslice[p],"Bio",which.biovars[b],".pdf",sep=""),width=80,height=80,pointsize=100,useDingbats = F)
        plot(spdata,cex=1,xlab="",ylab="",axes=FALSE,pch=16,col="lightgray")
        points(spdata,cex=1,pch=16,col=colorscheme[round(MESS_score[[p]][,b] - min(MESS_score[[p]][,b]) + 1)],xlim=c(-200,0),ylim=c(0,90))
        if(sum(ex_fossils[,1]==(timeslice[p] + 1))!=0){
          points(spfossils,cex=2,pch=16,col="black")
        }
        dev.off()
      }
    }
  }
  if(length(which.biovars) > 1){
    if(which.plot == "all" | which.plot == "mess"){
      for(p in 1:length(timeslice)){
        #transform spatial data
        spdata <- SpatialPoints(paleoclimate[[(timeslice[p] + 1)]][, 2:3])
        proj4string(spdata)  <- CRS("+init=epsg:4326")
        spdata <- spTransform(spdata, CRS("+init=epsg:26978"))
        if(sum(fossils[, 1] == (timeslice[p] + 1)) != 0){
          spfossils <- SpatialPoints(fossils[, 2:3])
          proj4string(spfossils)  <- CRS("+init=epsg:4326")
          spfossils <- spTransform(spfossils, CRS("+init=epsg:26978"))
          spfossils <- spfossils[fossils[, 1] == (timeslice[p] + 1), ]
        }
        #print plots
        pdf(paste("MESS",timeslice[p],"Multi.pdf",sep=""),width=80,height=80,pointsize=100,useDingbats = F)
        plot(spdata,cex=1,xlab="",ylab="",axes=FALSE,pch=16,col="lightgray")
        points(spdata,cex=1,pch=16,col=colorscheme[round(apply(MESS_score[[p]], 1, min) - min(MESS_score[[p]]) + 1)],xlim=c(-200,0),ylim=c(0,90))
        if(sum(ex_fossils[,1]==(timeslice[p] + 1))!=0){
          points(spfossils,cex=2,pch=16,col="black")
        }
        dev.off()
      }
    }
  }
  return(mess = MESS_score)
}
