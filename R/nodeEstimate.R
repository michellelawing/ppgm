#' @title nodeEstimate
#' @description This function estimates the ancestral character states for continuous characters given a model of evolution or using the best fit model of evolution from the fitContinuous function in the geiger package. The ancestral states are estimated using GLS described in Martins and Hansen (1997).
#' @usage nodeEstimate(treedata.obj, traitnum, model = "BM", bounds = list(), control = list(), plot.est = FALSE)
#' @param treedata.obj an object of the class "treedata".
#' @param traitnum the column number of the trait within the treedata object to be reconstructed.
#' @param model the model of evolution to use in the ancestral state reconstruction. Options are "estimate", "BM", "OU", "EB", "lambda", "kappa", "delta".
#' @param plot.est logical. whether or not to plot the traitgram of the estimated ancestor states.
#' @param bounds bounds used for the model, passes to \code{fitContinuous()}, uses default if none specified.
#' @param control setting used for optimization of the model likelihood. Passes to \code{fitContinuous()}.
#' @details See the \code{fitContinuous()} details for descriptions of the models of evolution and parameter estimation. \code{nodeEstimate()} currently supports the following models of evolution: Brownian motion (Felsenstein, 1973), Ornstein-Uhlenbeck (Butler and King, 2004), early-burst (Harmon et al., 2010), lambda (Pagel, 1999), kappa (Pagel, 1999), and delta (Pagel, 1999).
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
#' @seealso \code{fitContinuous()}
#' @references Butler, M. A. and King, A. A. (2004) Phylogenetic comparative analysis: a modeling approach for adaptive evolution. American Naturalist, 164:683-695.
#' @references Felsenstein, J. (1973) Maximum likelihood estimation of evolutionary trees from continuous characters. American Journal of Human Genetics, 25:471-492
#' @references Harmon, L. J. et al. (2010) Early bursts of body size and shape evolution are rare in comparative data. Evolution, 64:2385-2396
#' @references Martins, E. P. and Hansen, T. F. (1997) Phylogenies and the comparative method: a general approach to incorporating phylogenetic information into the analysis of interspecific data. American Naturalist, 149, 646–667.
#' @references Pagel M. (1999) Inferring the historical patterns of biological evolution. Nature, 401:877-884
#' @author A. Michelle Lawing, Alexandra F. C. Howard
#' @export
#' @examples
#' data(beastLeache)
#' data(occurrences)
#' sp_data_min<- tapply(occurrences[,4],occurrences$Species,min)
#' ex <- treedata(beastLeache[[1]], sp_data_min)
#' nodeEstimate(ex, 1, model = 'OU') #runs OU model


nodeEstimate <- function(treedata.obj, traitnum, model="BM", bounds=list(), control=list(), plot.est=FALSE) {
  require(geiger)
  x <- treedata.obj$data[,traitnum]
  phy <- treedata.obj$phy
  was.estimated <- FALSE
  fitted <- model
  if(model=="estimate"){
    models=c("BM","OU","EB","lambda","kappa","delta")
    BM=try(fitContinuous(phy,x,model="BM",bounds=bounds,control=control),silent=T)
    OU=try(fitContinuous(phy,x,model="OU",bounds=bounds,control=control),silent=T)
    EB=try(fitContinuous(phy,x,model="EB",bounds=bounds, control=control),silent=T)
    lambda=try(fitContinuous(phy,x,model="lambda",bounds=bounds,control=control),silent=T)
    kappa=try(fitContinuous(phy,x,model="kappa",bounds=bounds,control=control),silent=T)
    delta=try(fitContinuous(phy,x,model="delta",bounds=bounds,control=control),silent=T)
    trait.macroevo <- list()
    for (mod in 1:length(models)){
      if(is(eval(parse(text=models[mod])),"gfit")) {
        trait.macroevo[[mod]] <- eval(parse(text=models[mod]))$opt$aicc
      } else {
        trait.macroevo[[mod]] <- NA
      }
    }
    model <- models[which.min(unlist(trait.macroevo))]
    fitted <- list()
    for (mod in 1:length(models)){
      if(is(eval(parse(text=models[mod])),"gfit")) {
        fitted[[models[mod]]] <- eval(parse(text=models[mod]))$opt
      } else {
        fitted[[mod]] <- NA
        names(fitted[mod]) <-paste(models[mod])
      }
    }
    if (model=="OU" & !is(OU,"try-error")) {if(!is.na(OU$opt$alpha)){phy=rescale(phy,model="OU",OU$opt$alpha)}}
    if (model =="EB" & !is(EB,"try-error")){if(!is.na(EB$opt$a)){phy=rescale(phy,model="EB",EB$opt$a)}}
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
    if (model=="lambda") {fitted<-lambda<-try(fitContinuous(phy,x,model="lambda",bounds=bounds,control=control),silent=T)
    if(!is(lambda,"try-error")){if(!is.na(lambda$opt$lambda)){phy=rescale(phy,model="lambda",lambda$opt$lambda)}}}
    if(model=="kappa") {fitted<-kappa<-try(fitContinuous(phy,x,model="kappa",bounds=bounds,control=control),silent=T)
    if(!is(kappa,"try-error")){ if(!is.na(kappa$opt$kappa)){ phy=rescale(phy,model="kappa",kappa$opt$kappa)}}}
    if(model=="delta") {fitted<-delta<-try(fitContinuous(phy,x,model="delta",bounds=bounds,control=control),silent=T)
    if(!is(delta,"try-error")){ if(!is.na(delta$opt$delta)){phy=rescale(phy,model="delta",delta$opt$delta)}}}
  }
  #make sure phylo is phylo
  if(class(phy)=="phylo"){
    M <- dist.nodes(phy)
  }else{
    phy <- as.phylo(phy$phy)
    M <- dist.nodes(phy) #uses GLS on [rescaled] phylo
  }
  #Martins & Hansen (1997)
  nb.tip <- length(phy$tip.label)
  varAY <- M[-(1:nb.tip), 1:nb.tip]
  varY <- M[1:nb.tip,1:nb.tip]
  J<-array(1,dim=nb.tip)
  if(try(is.matrix(solve(varY)),silent=T)==TRUE){
    GrandMean<-J%*%solve(varY)%*%x / J%*%solve(varY)%*%J
    node.est<- varAY%*%solve(varY)%*%(x-c(GrandMean)) + GrandMean[1,1]
  } else {
    warning("In node.est(): singular matrix: using dist matrix without the rescale, revert to BM")
    M <- dist.nodes(treedata.obj$phy)
    varAY <- M[-(1:nb.tip), 1:nb.tip]
    varY <- M[1:nb.tip,1:nb.tip]
    GrandMean<-J%*%solve(varY)%*%x / J%*%solve(varY)%*%J
    node.est<- varAY%*%solve(varY)%*%(x-c(GrandMean)) + GrandMean[1,1]
    phy <- treedata.obj$phy
  }
  if(plot.est && !is.na(node.est)) {
    plot(dist.nodes(treedata.obj$phy)[,treedata.obj$phy$edge[1,1]],c(treedata.obj$data[,traitnum],node.est),xlab="Time",ylab="Trait",type="n")
    for(i in 1:length(treedata.obj$phy$edge[,1])){
      lines(dist.nodes(treedata.obj$phy)[,treedata.obj$phy$edge[1,1]][treedata.obj$phy$edge[i,]],c(treedata.obj$data[,traitnum],node.est)[treedata.obj$phy$edge[i,]])
    }
  }
  if(was.estimated){ return (list(model=model,est=node.est,phy=phy,fitted=fitted))} else {return (list(model=model,est=node.est,phy=phy,fitted=fitted$opt))}
}
