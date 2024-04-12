#' @title ppgmMESS
#' @description This creates a MESS map for given time slices, climate envelopes, and paleoclimate models.
#' @usage ppgmMESS(cem_min, cem_max, est, tree, fossils, timeslice, 
#' which.biovars, path = "", use.paleoclimate=TRUE, paleoclimateUser = NULL, 
#' which.plot = c("all","mess","none"))
#' @param cem_min the cem min output from the ppgm function. cbind() if there are multiple variables.
#' @param cem_max the cem max output from the ppgm function. cbind() if there are multiple variables.
#' @param est the node_est output from the ppgm function, in list format. [tree][1][min and max][no.of species]
#' @param fossils a matrix with four columns of age to the closest million year integer, longitude, and latitude, in that order, and rows that are entries for fossil occurrences.
#' @param timeslice the time in million of years ago to project MESS maps (0 to 20). can handle single timeslice or vector of times.
#' @param tree the phylogeny or multiple phylogenies that show the relationship between species
#' @param which.biovars the biovariable number(s) between 1 and 19. handles vectors.
#' @param path directory where plots should be stored
#' @param use.paleoclimate if left blank, default North America paleoclimate data is used. If FALSE, user submitted paleoclimate must be provided
#' @param paleoclimateUser list of data frames with paleoclimates, must be dataframes with columns: GlobalID, Longitude, Latitude, bio1, bio2,...,bio19.
#' @param which.plot "all" plots trait maps and MESS, "mess" plots MESS map, "none" does not plot
#' @details plots MESS maps of climate envelope model for specific time slices. Can either plot individual biovariables, or combined.
#' @importFrom utils data
#' @importFrom grDevices colorRamp
#' @importFrom geiger treedata
#' @importFrom grDevices pdf
#' @importFrom graphics points
#' @importFrom grDevices dev.off
#' @import sp
#' @export
#' @seealso \code{ppgm()}
#' @author A. Michelle Lawing, Alexandra F. C. Howard, Maria-Aleja Hurtado-Materon
#' @examples
#' data(sampletrees)
#' data(occurrences)
#' data(scel_fossils)
#' biooccu <- getBioclimVars(occurrences, which.biovars=1)
#' sp_data_min<- tapply(biooccu[,4],biooccu$Species,min)
#' sp_data_max<- tapply(biooccu[,4],biooccu$Species,max)
#' treedata_min <- treedata_max <- node_est <- list()
#' \dontrun{biofossils <- getBioclimVars(scel_fossils,which.biovars=1)
#' rownames(biofossils)<-paste("fossil",1:length(biofossils[,1]),sep="")
#' for (tr in 1:length(sampletrees)){
#'   treedata_min[[tr]] <- geiger::treedata(sampletrees[[tr]],sp_data_min,sort=TRUE,warnings=F)
#'   treedata_max[[tr]] <- geiger::treedata(sampletrees[[tr]],sp_data_max,sort=TRUE,warnings=F)
#'   full_est <- nodeEstimateFossils(treedata_min[[tr]],treedata_max[[tr]])
#'   node_est[[tr]] <- full_est$est
#' }
#' relist1 <- lapply(lapply(1:length(sampletrees), function(x) array(unlist(node_est), 
#' dim = c(2, 52, 3, 100))[,,1,x]), list)
#' ppgmMESS(sp_data_min, sp_data_max, est=list(relist1) ,sampletrees, fossils=biofossils,
#' timeslice=10, which.biovars=1, which.plot="all")
#' }

ppgmMESS <- function(cem_min, cem_max, est, tree, fossils, timeslice, which.biovars, path = "", use.paleoclimate=TRUE, paleoclimateUser = NULL, which.plot = c("all","mess","none")){
  #load paleoclimate data
  if(use.paleoclimate) {
    paleoclimate <- paleoclimate #uses paleoclimate data from package
  } else {
    if(is.null(paleoclimateUser)) {
      stop("paleoclimateUser argument must be provided when use.paleoclimate is FALSE.") #uses user inputted paleoclimate
    } else {
      paleoclimate <- paleoclimateUser
    }
  }
  colorscheme <- grDevices::colorRampPalette(c("blue","cyan","greenyellow","yellow","darkorange","red"))(250)
  MESS_score <- as.list(array(NA,dim=length(timeslice)))
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
      treedata_min[[tr]] <- geiger::treedata(tree[[tr]], data = data_min, sort = TRUE, warnings = F)
      treedata_max[[tr]] <- geiger::treedata(tree[[tr]], data = data_max, sort = TRUE, warnings = F)
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
      spdata <- sp::SpatialPoints(paleoclimate[[(timeslice[p] + 1)]][, 2:3])
      sp::proj4string(spdata)  <- sp::CRS("+init=epsg:4326")
      spdata <- sp::spTransform(spdata, CRS("+init=epsg:26978"))
      if(sum(fossils[, 1] == (timeslice[p] + 1)) != 0){
        spfossils <- sp::SpatialPoints(fossils[, 2:3])
        sp::proj4string(spfossils)  <- sp::CRS("+init=epsg:4326")
        spfossils <- sp::spTransform(spfossils, sp::CRS("+init=epsg:26978"))
        spfossils@data[,1] <- fossils[, 1]
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
        grDevices::pdf(paste("MESS",timeslice[p],"Bio",which.biovars[b],".pdf",sep=""),width=80,height=80,pointsize=100,useDingbats = F)
        plot(spdata,cex=1,xlab="",ylab="",axes=FALSE,pch=16,col="red")
        graphics::points(spdata,cex=1,pch=16,col=colorscheme[round(MESS_score[[p]][,b] - min(MESS_score[[p]][,b]) + 1)],xlim=c(-200,0),ylim=c(0,90))
        if(sum(fossils[,1]==(timeslice[p] + 1))!=0){
          points(spfossils,cex=2,pch=16,col="black")
        }
        grDevices::dev.off()
      }
    }
  }
  if(length(which.biovars) > 1){
    if(which.plot == "all" | which.plot == "mess"){
      for(p in 1:length(timeslice)){
        #transform spatial data
        spdata <- sp::SpatialPoints(paleoclimate[[(timeslice[p] + 1)]][, 2:3])
        sp::proj4string(spdata)  <- sp::CRS("+init=epsg:4326")
        spdata <- sp::spTransform(spdata, sp::CRS("+init=epsg:26978"))
        if(sum(fossils[, 1] == (timeslice[p] + 1)) != 0){
          spfossils <- sp::SpatialPoints(fossils[, 2:3])
          sp::proj4string(spfossils)  <- sp::CRS("+init=epsg:4326")
          spfossils <- sp::spTransform(spfossils, sp::CRS("+init=epsg:26978"))
          spfossils@data[,1] <- fossils[, 1]
          spfossils <- spfossils[fossils[, 1] == (timeslice[p] + 1), ]
        }
        #print plots
        grDevices::pdf(paste("MESS",timeslice[p],"Multi.pdf",sep=""),width=80,height=80,pointsize=100,useDingbats = F)
        plot(spdata,cex=1,xlab="",ylab="",axes=FALSE,pch=16,col="red")
        graphics::points(spdata,cex=1,pch=16,col=colorscheme[round(apply(MESS_score[[p]], 1, min) - min(MESS_score[[p]]) + 1)],xlim=c(-200,0),ylim=c(0,90))
        if(sum(fossils[,1]==(timeslice[p] + 1))!=0){
          points(spfossils,cex=2,pch=16,col="black")
        }
        grDevices::dev.off()
      }
    }
  }
  return(mess = MESS_score)
}
