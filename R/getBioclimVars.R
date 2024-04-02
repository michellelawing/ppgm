#' @title getBioclimVars
#' @description This function retrieves the bioclimatic variables described in Nix & Busby (1986) for the specified variables and the specified time period.
#' @usage getBioclimVars(occurrences, which.biovars=c(1,12), use.paleoclimate=TRUE, paleoclimateUser=NULL)
#' @param occurrences a matrix or data.frame with three columns and rows to represent individuals. The first column must be species name for extant occurrences or the age in closest Ma for fossil occurrences. Second and third column must be Longitude and Latitude.
#' @param which.biovars a vector of the numbers of the bioclimatic variables that should be returned. The bioclimatic variables number correspond to the table at (https://www.worldclim.org/data/bioclim.html).
#' @param use.paleoclimate if left blank, default North America paleoclimate data is used. If FALSE, user submitted paleoclimate must be provided
#' @param paleoclimateUser list of data frames with paleoclimates, must be dataframes with columns: GlobalID, Longitude, Latitude, bio1, bio2,...,bio19.
#' @details The occurrences argument should contain all extant or all fossils. Columns should be in the format: Species, Longitude, Latitude for extant data.
#' @details If using the provided paleoclimate data:
#' @details Modern time period uses the Hijmans et al. (2005) high resolution climate interpolations.
#' @details The time period 10 Ma uses the GCM by Micheels et al (2011) for the Tortonian.
#' @details The time period 15 Ma uses the GCM by Krapp & Jungclaus (2011) for the Middle Miocene.
#' @details For the one million year intervals outside the modern and past GCMs, the climate was interpolated based on the benthic marine foram stable oxygen isotope ratio curve from Ruddiman et al 1989. The scale of these variables is at a 50 km equidistant point grain size corresponding to Polly XX.
#' @return Returns a data frame with the original occurrences input appended with columns of bioclimate variables as specified. If fossils are included, the returned bioclimate variables are from the closest 1 Ma interval of isotopically scaled climate.
#' @author A. Michelle Lawing, Alexandra F. C. Howard, Maria-Aleja Hurtado-Materon
#' @references Polly;
#' @references Hijmans, R. J. et al. (2005) Very high resolution interpolated climate surfaces for global land areas
#' @references Krapp, M. and Jungclaus, J. H. (2011) The Middle Miocene climate as modeled in an atmosphere-ocean-biosphere model. Climate of the Past 7(4):1169-1188
#' @references Micheels, A. et al. (2011) Analysis of heat transport mechanisms from a Late Miocene model experiment with a fully-coupled atmosphere-ocean general circulation model. Palaeogeography, Palaeoclimatology, Palaeocology 304: 337-350
#' @references Nix, H. and Busby, J. (1986) BIOCLIM, a bioclimatic analysis and prediction system. CSIRO annual report. CSIRO Division of Water and Land Resources, Canberra.
#' @references Ruddiman, W. F. et al. (1989) Pleistocene evolution: Northern hemisphere ice sheets and North Atlantic Ocean. Paleoceanography 4: 353-412
#' @export
#' @examples
#' require(fields)
#' data(occurrences)
#' subsetoccur <- head(occurrences[,c(1:3)], n=100L)
#' biooccur <- getBioclimVars(subsetoccur,which.biovars=c(3,5))
#' #returns data frame with bioclimate variables 3 and 5 for occurrence data


getBioclimVars <- function(occurrences, which.biovars=c(1,12), use.paleoclimate=TRUE, paleoclimateUser=NULL){
  require (fields)
  if(use.paleoclimate) {
    data(paleoclimate) #uses paleoclimate data from package
  } else {
    if(is.null(paleoclimateUser)) {
      stop("paleoclimateUser argument must be provided when use.paleoclimate is FALSE.") #uses user inputted paleoclimate
    } else {
      paleoclimate <- paleoclimateUser
    }
  }
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
      temp <- array(NA, dim=c(length(occurrences[,1]), length(which.biovars)))
      for(i in 1:length(occurrences[,1])){
        anothertemp <- as.matrix(paleoclimate[[as.integer(occurrences[i,1])]][,2:3])
        calc_dist <- rdist.earth(anothertemp,t(as.matrix(occurrences[i,2:3])))
        temp[i,] <- unlist(paleoclimate[[as.integer(occurrences[i,1])]][which.min(calc_dist),which.biovars+3])
      }
      occurrences <- cbind(occurrences,temp)
      colnames(occurrences) <- c("Age", "Longitude", "Latitude", paste("bio",which.biovars,sep=""))
      return(occurrences)
    }
  }
}
