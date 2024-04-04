#' @title ppgm
#' @description ppgm makes a paleophylogeographic species distribution model using the bioclimate envelope method for a specified time period. Currently, models are only available for North America.
#' @usage ppgm(occurrences, fossils = FALSE, trees, fossils.edges = FALSE, model = "BM", permut = 1, only.biovars = TRUE, which.biovars = c(1:19), path = "", plot.TraitGram = F, plot.AnimatedMaps = F, plot.GeoRates = F, bounds = list(), control = list(), use.paleoclimate = TRUE, paleoclimateUser = NULL, verbose = TRUE)
#' @param occurrences a matrix with three columns of species name, longitude, and latitude, in that order, and rows that are entries for species occurrences. The bioclimate variables can be included for each occurrence in following columns. They must be in order 1 through 19.
#' @param fossils a matrix with four columns of age to the closest million year integer, longitude, and latitude, in that order, and rows that are entries for fossil occurrences. The bioclimate variables can be included for each occurrence in following columns. They must be in order 1 through 19. All 19 variables must be included at this stage, variable selection is done with the argument: "which.biovars".
#' @param trees phylogenies of species from first column of occurrences argument. Object of class multiphylo.
#' @param fossils.edges a vector of edges that the fossils belong to. Must be in the same order of the fossils argument. If fossils.edges is false, the the function randomly assigns the location of the fossils depending on the age (see details for more information).
#' @param model the model of evolution to use to estimate ancestor nodes. Argument is passed onto to function nodeEstimate.
#' @param permut the number of times to randomly place fossils in phylogeny and estimate ancestor states.
#' @param only.biovars logical. If FALSE, user must include biovariables in occurrence object.
#' @param which.biovars a vector with the biovars to include in model (see www.worldclim.org for a list of biovars). If "ALL", then all 19 biovars are included in analysis.
#' @param path path to the directory where the results should be saved.
#' @param plot.TraitGram logical. Whether to plot a TraitGram
#' @param plot.AnimatedMaps logical. Whether to plot AnimatedMaps. Requires ImageMagick to be installed on the system.
#' @param plot.GeoRates logical. Whether to plot GeoRates
#' @param bounds parameters for the evolutionary model selected. If none are supplied the default is used
#' @param control settings used for optimisation of model likelihood. Passes to \code{geiger::fitContinuous}
#' @param use.paleoclimate if left blank, default North America paleoclimate data is used. If FALSE, user submitted paleoclimate must be provided
#' @param paleoclimateUser list of data frames with paleoclimates, must be dataframes with columns: GlobalID, Longitude, Latitude, bio1, bio2,...,bio19.
#' @param verbose default true, returns all outputs. If FALSE then returns only climate envelopes and geographic data
#' @details If the 19 bioclimate variables are not supplied with the occurrences or with the fossils, they will be extracted from the closest 50km point location in the modern or paleoclimate maps that are loaded in with this function. The paleoclimate maps are isotopically scaled between general circulation models (see Lawing and Polly 2011; Rodder et al. 2013) and modern climate (see Hijmans et al. 2005). The fossils paleoclimate data is extracted to the closest million year paleoclimate map. Paleoclimate maps are derived at one million year intervals for the past 20 Ma. The tree (phylogeny) should be dichotomous and the species names should match the names in the first column of the occurrences argument.
#' For Michelle: fix to work with the nice plots on multiple trees, and check and see if it works when no suitable habitat is identified. also, make option to print any animated lineage through time?
#' @return \code{cem} Estimate of climate envelope for each species in present time. A data frame containing species and min mean and max of biovars specified with \code{which.biovars}.
#' @return geo_move data frame of RateGeoCenter and RateGeoSize
#' @return change_geo_center array of change in geographic center of suitable climate for each lineage
#' @return change_geo_size array of change in geographic size of suitable climate for each lineage
#' @return time_int matrix array of time intervals
#' @return treedata_min list of trees with minimum bioclimatic variables
#' @return treedata_max list of trees with maximum bioclimatic variables
#' @return model_min list of trees with minimum fitted model as specified in \code{model}
#' @return model_max list of trees with maximum fitted model as specified in \code{model}
#' @return node_est list of traits at each node for all trees, min and max for each species. As estimated by nodeEstimate and nodeEstimateFossils
#' @author A. Michelle Lawing, Alexandra F. C. Howard, Maria A. Hurtado-Materon
#' @importFrom utils data
#' @importFrom ape is.binary
#' @importFrom ape multi2di
#' @importFrom geiger treedata
#' @export
#' @examples
#' utils::data(beastLeache)
#' utils::data(occurrences)
#' bounds <- list(sigsq = c(min = 0, max = 1000000))
#' ex_mytree <- sample(beastLeache,size=10) #small number of trees used for example only to save processing time, we recommend at least 100 trees
#' test_ppgm <- ppgm(occurrences = occurrences,trees = ex_mytree, model = "BM", which.biovars = c(1), bounds = bounds, control = list(niter = 20), plot.TraitGram=T)


ppgm <- function(occurrences, fossils = FALSE, trees, fossils.edges = FALSE, model = "BM", permut = 1, only.biovars = TRUE,
                which.biovars = c(1:19), path = "", plot.TraitGram = F, plot.AnimatedMaps = F, plot.GeoRates = F,
                bounds = list(), control = list(), use.paleoclimate = TRUE, paleoclimateUser = NULL, verbose = TRUE){
  if(class(trees)=="phylo"){
    stop("ERROR: only one tree supplied. Please use ppgmConsensus")
  }
  #calculate the alpha.trans, which is the transparency for all the trees plotted on top of each other
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
  sp_data_min<-sapply(4:(length(which.biovars)+3),function(x) tapply(occurrences[,x],occurrences$Species,min))  #min of biovariable per species
  sp_data_mean<-sapply(4:(length(which.biovars)+3),function(x) tapply(occurrences[,x],occurrences$Species,mean))
  sp_data_max<-sapply(4:(length(which.biovars)+3),function(x) tapply(occurrences[,x],occurrences$Species,max))
  colnames(sp_data_mean)<-which.biovars
#  if(plot.BumpChart){
#    plotBumpChart(sp_data_mean,path)  #broken?? Error: breaks and labels have unequal lengths
#  }
  for(tr in 1:length(trees)){
    #Check if phylogeny is dichotomous, if not, make it dichotomous
    if(length(trees)==1){if(!ape::is.binary(trees[[tr]])){trees[[tr]]<-ape::multi2di(trees[[tr]])}}
    #make treedata object for bioclimate envelopes and phylogeny
    treedata_min[[tr]]<-geiger::treedata(trees[[tr]],sp_data_min,sort=TRUE,warnings=F)  #matches species in tree with species in data
    treedata_max[[tr]]<-geiger::treedata(trees[[tr]],sp_data_max,sort=TRUE,warnings=F)
    colnames(treedata_min[[tr]]$data)<-colnames(treedata_max[[tr]]$data)<-paste("bio",which.biovars,sep="")  #labels biovars
    #to estimate nodes, place fossils randomly or as specified on edges from fossils.edges argument
    full_est <- list()
    for(pr in 1:permut){
      full_est[[pr]] <- nodeEstimateFossils(treedata_min=treedata_min[[tr]],treedata_max=treedata_max[[tr]],fossils=fossils,fossils.edges=fossils.edges,model=model,bounds=bounds,control=control)
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
    plotTraitGramMultiPhylo(treedata_min,treedata_max,node_est,fossils=fossils,which.biovars=which.biovars,path=path)
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
  print_table_min <- NA
  print_table_max <- NA
  if(model=="estimate"){
    models <- c("BM", "OU", "EB", "lambda", "kappa", "delta")
    for(trees in 1:length(model_min)){
      for(traits in 1:length(model_min[[1]][[1]])){
        temp_min <- cbind(unlist(sapply(models, function(z) model_min[[trees]][[1]][[traits]]$fitted[[z]]['aicc'])))
        temp_max <- cbind(unlist(sapply(models, function(z) model_min[[trees]][[1]][[traits]]$fitted[[z]]['aicc'])))
        colnames(temp_min) <- colnames(temp_max) <- paste("aicc", trees, 1, traits, sep = "")
        if(trees == 1){
          print_table_min <- temp_min
          print_table_max <- temp_max
        }
        else {
          print_table_min <- cbind(print_table_min, temp_min)
          print_table_max <- cbind(print_table_max, temp_max)
        }
        rownames(print_table_min) <- rownames(print_table_max) <- models
      }
    }
    print_table_min <- print_table_min[-4,]
    print_table_max <- print_table_max[-4,]
  }
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
                node_est=node_est,
                print_table_min=print_table_min,
                print_table_max=print_table_max))
  }  else{
    return(list(cem=cem,
                geo_move=geo_move,
                change_geo_center=geo_center,
                change_geo_size=geo_size,
                time_int=time_int))
  }
}

