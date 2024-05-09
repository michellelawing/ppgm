#' @title ppgmConsensus
#' @description ppgm makes a paleophylogeographic species distribution model using the bioclimate envelope method for a specified time period. consensus version
#' @usage ppgmConsensus(occurrences, fossils = FALSE, trees, 
#' fossils.edges = FALSE, model = "BM", permut = 1, only.biovars = TRUE, 
#' which.biovars = c(1:19), path = "", plot.TraitGram = FALSE, 
#' plot.AnimatedMaps = FALSE, plot.GeoRates = FALSE, bounds = list(), 
#' control = list(), use.paleoclimate = TRUE, paleoclimateUser = NULL, 
#' layerAge = c(0:20), verbose = TRUE)
#' @param occurrences a matrix with three columns of species name, longitude, and latitude, in that order, and rows that are entries for species occurrences. The bioclimate variables can be included for each occurrence in following columns. They must be in order 1 through 19.
#' @param fossils a matrix with four columns of min age, max age, longitude, and latitude, in that order, and rows that are entries for fossil occurrences. The bioclimate variables can be included for each occurrence in following columns. They must be in order 1 through 19. All 19 variables must be included at this stage, variable selection is done with the argument: "which.biovars".
#' @param trees phylogeny of species from first column of occurrences argument. Object of class phylo.
#' @param fossils.edges a vector of edges that the fossils belong to. Must be in the same order of the fossils argument. If fossils.edges is false, the the function randomly assigns the location of the fossils depending on the age (see details for more information).
#' @param model the model of evolution to use to estimate ancestor nodes. Argument is passed onto to function nodeEstimate.
#' @param permut the number of times to randomly place fossils in phylogeny and estimate ancestor states.
#' @param only.biovars logical. If FALSE, user must include biovariables in occurrence object.
#' @param which.biovars a vector with the biovars to include in model (see www.worldclim.org for a list of biovars). If "ALL", then all 19 biovars are included in analysis.
#' @param path path to the directory where the results should be saved.
#' @param plot.TraitGram logical. Whether to plot a TraitGram
#' @param plot.AnimatedMaps Logical. Whether to plot AnimatedMaps. Requires ImageMagick to be installed on the system.
#' @param plot.GeoRates logical. Whether to plot GeoRates
#' @param bounds parameters for the evolutionary model selected. If none are supplied the default is used
#' @param control settings used for optimisation of model likelihood. Passes to \code{geiger::fitContinuous}
#' @param use.paleoclimate if left blank, default North America paleoclimate data is used. If FALSE, user submitted paleoclimate must be provided
#' @param paleoclimateUser list of data frames with paleoclimates, must be dataframes with columns: GlobalID, Longitude, Latitude, bio1, bio2,...,bio19.
#' @param layerAge vector with the ages of the paleoclimate dataframes, if using user submitted paleoclimate data
#' @param verbose default true, returns all outputs. If FALSE then returns only climate envelopes and geographic data
#' @details If the 19 bioclimate variables are not supplied with the occurrences or with the fossils, they will be extracted from the closest 50km point location in the modern or paleoclimate maps that are loaded in with this function. The paleoclimate maps are isotopically scaled between general circulation models (see Lawing and Polly 2011; Rodder et al. 2013) and modern climate (see Hijmans et al. 2005). The fossils paleoclimate data is extracted to the closest million year paleoclimate map. Paleoclimate maps are derived at one million year intervals for the past 20 Ma. The tree (phylogeny) should be dichotomous and the species names should match the names in the first column of the occurrences argument.
#' @return \code{cem} Estimate of climate envelope for each species in present time. A data frame containing species and min mean and max of biovars specified with \code{which.biovars}.
#' @return \code{geo_move} data frame of RateGeoCenter and RateGeoSize
#' @return \code{change_geo_center} array of change in geographic center of suitable climate for each lineage
#' @return \code{change_geo_size} array of change in geographic size of suitable climate for each lineage
#' @return \code{time_int} matrix array of time intervals
#' @return \code{treedata_min} list of trees with minimum bioclimatic variables
#' @return \code{treedata_max} list of trees with maximum bioclimatic variables
#' @return \code{node_est} list of traits at each node for all trees, min and max for each species. As estimated by nodeEstimate and nodeEstimateEnvelopes
#' @return \code{aicmin} if model is estimated, table of aic values for minimum trait values
#' @return \code{aicmax} if model is estimated, table of aic values for maximum trait values
#' @author A. Michelle Lawing, Alexandra F. C. Howard, Maria A. Hurtado-Materon
#' @importFrom utils data
#' @importFrom ape is.binary
#' @importFrom ape multi2di
#' @importFrom geiger treedata
#' @export
#' @examples
#' data(sampletrees)
#' data(occurrences)
#' data(scel_fossils)
#' #set bounds for analysis
#' bounds <- list(sigsq = c(min = 0, max = 1000000))
#' ex_mytree <- sampletrees[[3]] #single tree
#' \donttest{test_fossil_con <- ppgmConsensus(occurrences = occurrences, 
#' fossils = scel_fossils, trees = ex_mytree, fossils.edges = FALSE, model = "BM", 
#' permut = 5, which.biovars = 1, bounds = bounds, control = list(niter = 20))}

ppgmConsensus <- function(occurrences, fossils = FALSE, trees, fossils.edges = FALSE, model = "BM", permut = 1, only.biovars = TRUE,
                          which.biovars = c(1:19), path = "", plot.TraitGram = FALSE, plot.AnimatedMaps = FALSE, plot.GeoRates = FALSE,
                          bounds = list(), control = list(), use.paleoclimate = TRUE, paleoclimateUser = NULL, layerAge=c(0:20), verbose = TRUE){
  #assign rownames to fossils
  if(length(fossils)!=1){rownames(fossils)<-paste("fossil",1:length(fossils[,1]),sep="")}
  #Extract bioclim variables for species occurrences, if not supplied by user
  if(only.biovars) {
    occurrences<-getBioclimVars(occurrences, which.biovars, use.paleoclimate, paleoclimateUser, layerAge)
  }
  #check the occurrences variable conforms with requirements
  if(length(occurrences[1, ]) < 4) {
    return ("ERROR: There are less than four columns for the occurrences input. Check that occurrences input conforms to requirements.")
  }
  #load paleoclimate data
  if(use.paleoclimate) {
    utils::data(paleoclimate) #uses paleoclimate data from package
  } else {
    if(is.null(paleoclimateUser)) {
      stop("paleoclimateUser argument must be provided when use.paleoclimate is FALSE.") #uses user inputted paleoclimate
    } else {
      paleoclimate <- paleoclimateUser
    }
  }
  #create empty lists for outputs
  geo_center<-array(NA,dim=c(length(unique(occurrences$Species)),length(paleoclimate),length(paleoclimate)))
  geo_size<-array(NA,dim=c(length(unique(occurrences$Species)),length(paleoclimate),length(paleoclimate)))
  time_int<-array(NA,dim=c(length(paleoclimate),length(paleoclimate)))
  treedata_min<-list()
  treedata_max<-list()
  node_est<-list()
  model_min<-list()
  model_max<-list()
  envelope<-list()
  #get bioclimate envelopes
  sp_data_min<-sapply(4:(length(which.biovars)+3),function(x) tapply(occurrences[,x],occurrences$Species,min))  #min of biovariable per species
  sp_data_mean<-sapply(4:(length(which.biovars)+3),function(x) tapply(occurrences[,x],occurrences$Species,mean))
  sp_data_max<-sapply(4:(length(which.biovars)+3),function(x) tapply(occurrences[,x],occurrences$Species,max))
  colnames(sp_data_mean)<-which.biovars
  #Check if phylogeny is dichotomous, if not, make it dichotomous
  if(!ape::is.binary(trees)){trees<-ape::multi2di(trees)}
  #make treedata object for bioclimate envelopes and phylogeny
  treedata_min<-geiger::treedata(trees,sp_data_min,sort=TRUE,warnings=F)  #matches species in tree with species in data
  treedata_max<-geiger::treedata(trees,sp_data_max,sort=TRUE,warnings=F)
  colnames(treedata_min$data)<-colnames(treedata_max$data)<-paste("bio",which.biovars,sep="")  #labels biovars
  #to estimate nodes, place fossils randomly or as specified on edges from fossils.edges argument
  full_est <- list()
  for(pr in 1:permut){
    full_est[[pr]] <- nodeEstimateEnvelopes(treedata_min=treedata_min, treedata_max=treedata_max, fossils=fossils, fossils.edges=fossils.edges, model=model, bounds=bounds, control=control, which.biovars=which.biovars, paleoclimateUser=paleoclimateUser, use.paleoclimate=use.paleoclimate, layerAge=layerAge)
  }
  #####################
  node_est <- lapply(1:permut, function(p) full_est[[p]]$est)
  model_min <- lapply(1:permut, function(p) full_est[[p]]$min_model)
  model_max <- lapply(1:permut, function(p) full_est[[p]]$max_model)
  #get bioclimate envelopes for species and nodes
  envelope <- getEnvelopes(treedata_min, treedata_max, node_est)
  #get data from geo displacement
  temp <- getGeoRate(envelope, tree=trees, which.biovars=which.biovars, layerAge=layerAge)
  geo_center <- temp$geo_center
  geo_size <- temp$geo_size
  time_int<-temp$time_int
  #plot permutations
  if(plot.TraitGram){
    plotTraitGram(treedata_min, treedata_max, node_est, fossils=fossils, which.biovars=which.biovars, path=path, paleoclimateUser=paleoclimateUser, use.paleoclimate=use.paleoclimate, layerAge=layerAge)
  }
  #plot animated maps
  if(plot.AnimatedMaps){
    plotAnimatedPPGM(envelope, tree=trees, filename="ppgm.gif", which.biovars=which.biovars, path=path)
  }
  #plot geo rates
  if(plot.GeoRates){
    plotGeoRatesCon(geo_center, geo_size, time_int, trees, path=path)
  }
  lineage_geo_center=array(NA,dim=c(1,length(geo_center[,1,1])))
  for(l in 1:length(geo_center[,1,1])){
    #rate in km/year
    #need to consider rates and time scale issue
    lineage_geo_center[l]=(mean((geo_center[l,,]/time_int)[lower.tri(geo_center[l,,]/time_int)],na.rm=T))
  }

  lineage_geo_size=array(NA,dim=c(length(trees),length(geo_size[,1,1])))
  for(l in 1:length(geo_size[,1,1])){
    #rate in km/year
    lineage_geo_size[l]=(mean((geo_size[l,,]/time_int)[lower.tri(geo_size[l,,]/time_int)],na.rm=T))
  }
  cem<-data.frame(sp_data_min,sp_data_mean,sp_data_max)
  names(cem)<-c(paste(colnames(sp_data_mean),"Min",sep=""),paste(colnames(sp_data_mean),"Mean",sep=""),paste(colnames(sp_data_mean),"Max",sep=""))
  geo_move<-data.frame(colMeans(lineage_geo_center,na.rm=T),colMeans(lineage_geo_size,na.rm=T))
  names(geo_move)<-c("RateGeoCenter","RateGeoSize")
  print_table_min <- list()
  print_table_max <- list()
  AICcmin <- AICcmax <- list()
  if(model=="estimate"){
    models <- c("BM", "OU", "EB", "lambda", "kappa", "delta")
    for(traits in 1:length(model_min[[1]])){
      clean<-list()
      for(trees in 1:length(model_min)){
        temp_min <- cbind(unlist(sapply(models, function(z) model_min[[1]][[traits]]$fitted[[z]]['aicc'])))
        temp_max <- cbind(unlist(sapply(models, function(z) model_min[[1]][[traits]]$fitted[[z]]['aicc'])))
        colnames(temp_min) <- colnames(temp_max) <- paste("aicc", trees, 1, traits, sep = "")
        print_table_min <- as.vector(temp_min)
        print_table_max <- as.vector(temp_max)
      }
      clean <- as.data.frame(print_table_min) 
      clean2 <- as.data.frame(print_table_max)
      colnames(clean) <- colnames(clean2) <- paste("tree1",sep="")
      rownames(clean) <- rownames(clean2) <- models
      AICcmin[[traits]] <- clean
      AICcmax[[traits]] <- clean2
    }
    names(AICcmin) <- names(AICcmax) <- paste("bio",which.biovars,sep="")
  }
  if(verbose){
    return(list(cem=cem,
                geo_move=geo_move,
                change_geo_center=geo_center,
                change_geo_size=geo_size,
                time_int=time_int,
                treedata_min=treedata_min,
                treedata_max=treedata_max,
                node_est=node_est,
                AICcmin=AICcmin,
                AICcmax=AICcmax))
  }  else{
    return(list(cem=cem,
                geo_move=geo_move,
                change_geo_center=geo_center,
                change_geo_size=geo_size,
                time_int=time_int))
  }
}

