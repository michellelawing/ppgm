library(knitr)
library(ape)
library(phytools)
library(geiger)
library(DBI)
library(RMySQL)
library(raster)
library(rgdal)
library(phangorn)
source("R/source_functions.R")

permut<-10000

load("data/trans_env_data.Rdata")
load("data/env_data.Rdata")
load("data/fossil_data.Rdata")

####################begin giant random fossil place - trait loop###################
#macroevolution, specify function to get trait of interest...e.g., mean, min, max...
sp_data_min<-sapply(4:22,function(x) tapply(env_data[,x],env_data$Species,min))
sp_data_max<-sapply(4:22,function(x) tapply(env_data[,x],env_data$Species,max))
colnames(sp_data_min)<-colnames(sp_data_max)<-paste("Bio",1:19,sep="")
#make names match
rownames(sp_data_min)[c(8,10,18,30,41)]<-rownames(sp_data_max)[c(8,10,18,30,41)]<-c("edwardtaylori90","gadovi","jarrovi","ochoterrane","smithii")

###read in one or a bunch of trees
#mytree<-read.tree(file="beastLeache.nex")
load("data/beastLeache.RData")
mytree <- beastLeache
which_run<-sample(1:length(mytree),1)
sample_mytree<-mytree[which_run]
#clean up
rm(mytree)

####here is the random fossil placement loop start
node_min<-array(NA,dim=c(permut,52,19))
node_max<-array(NA,dim=c(permut,52,19))
for (p in 1:permut) {
  treedata_min<-treedata(sample_mytree[[1]],sp_data_min,sort=TRUE,warnings=TRUE)
  treedata_max<-treedata(sample_mytree[[1]],sp_data_max,sort=TRUE,warnings=TRUE) 
  for(i in 14:length(fossil_data[,1])){
    treedata_min<-treedata(addfossil(treedata_min$phy,fossil_data[i,1],fossil_data[i,1]+1,paste("fossil",i,sep=" ")),rbind(treedata_min$data,fossil_data[i,4:22]),sort=TRUE,warnings=TRUE)
    treedata_max<-treedata(addfossil(treedata_max$phy,fossil_data[i,1],fossil_data[i,1]+1,paste("fossil",i,sep=" ")),rbind(treedata_max$data,fossil_data[i,4:22]),sort=TRUE,warnings=TRUE)
  }
  treedata_min$phy$node.label<-treedata_max$phy$node.label<-((length(treedata_min$phy$tip)+1):((length(treedata_min$phy$tip)*2)-1))
  ####here is the trait loop start
  for(tr in 1:19){ 
    traitnum<-tr
    tmin<-node.estimate(treedata_min,traitnum,model="estimate",plot.est=FALSE)$est
    tmax<-node.estimate(treedata_max,traitnum,model="estimate",plot.est=FALSE)$est
    rm_fossils<-sapply(54:57, function(x) Ancestors(treedata_min$phy,x)[1]-57)
    cntr<-array(NA,dim=4)
    cntr<-sapply(1:length(rm_fossils),function(x) if(duplicated(rm_fossils)[x]){cntr[x]<-2} else {cntr[x]<-1})
    node_min[p,,tr]<-tmin[-sapply(1:4, function(x) Ancestors(treedata_min$phy,(53+x))[cntr[x]]-57)]
    rm_fossils<-sapply(54:57, function(x) Ancestors(treedata_max$phy,x)[1]-57)
    cntr<-array(NA,dim=4)
    cntr<-sapply(1:length(rm_fossils),function(x) if(duplicated(rm_fossils)[x]){cntr[x]<-2} else {cntr[x]<-1})
    node_max[p,,tr]<-tmax[-sapply(1:4, function(x) Ancestors(treedata_max$phy,(53+x))[cntr[x]]-57)]
  }
}
treedata_min<-treedata(drop.tip(treedata_min$phy,54:57),treedata_min$data,warnings=FALSE)
treedata_max<-treedata(drop.tip(treedata_max$phy,54:57),treedata_max$data,warnings=FALSE)

##back transform data
for(i in 1:19){
  node_min[,,i]<-exp(node_min[,,i])-1
  treedata_min$data[,i]<-exp(treedata_min$data[,i])-1
  node_max[,,i]<-exp(node_max[,,i])-1
  treedata_max$data[,i]<-exp(treedata_max$data[,i])-1
  if(!is.na(trans_env_data[i])){
    node_min[,,i]<-node_min[,,i]-abs(trans_env_data[i])
    treedata_min$data[,i]<-treedata_min$data[,i]-abs(trans_env_data[i])
    node_max[,,i]<-node_max[,,i]-abs(trans_env_data[i])
    treedata_max$data[,i]<-treedata_max$data[,i]-abs(trans_env_data[i])
  }
}

save(node_min,file=paste("node_min",which_run,".Rdata",sep=""))
save(node_max,file=paste("node_max",which_run,".Rdata",sep=""))
save(treedata_min,file=paste("treedata_min",which_run,".Rdata",sep=""))
save(treedata_max,file=paste("treedata_max",which_run,".Rdata",sep=""))
