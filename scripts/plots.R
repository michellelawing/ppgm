library(maps)
library(RColorBrewer)
library(BiodiversityR)

par(mfrow=c(1,3),mar=c(5,4,1,0.5))

##############################################################################
#Species Occurrences and Richness
##############################################################################
load("data/occurrences.RData")
data.df<-occurrences
data.df <- subset(data.df, subset=(Longitude >= -180 & Latitude >= -90))
data.df$X <- floor(data.df$Longitude)
data.df$Y <- floor(data.df$Latitude)
data.df$Cell <- data.df$X + 360 * data.df$Y
counts <- by(data.df, data.df$Cell, function(d) c(d$X[1], d$Y[1], length(unique(d$Species))))
counts.m <- matrix(unlist(counts), nrow=3)
rownames(counts.m) <- c("X", "Y", "Count")
count.max <- max(counts.m["Count",])
colorscheme <- colorRampPalette(c("lightblue","yellow","orange","red","darkred"))(count.max)
plot(counts.m["X",], counts.m["Y",], type="n", xlab="Longitude", ylab="Latitude")
map("world",fill=T, col="lightgray",add=T)
map("state",add=T)
points(counts.m["X",] + 1/2, counts.m["Y",] + 1/2, cex=1, pch = 19, col=colorscheme[counts.m["Count",]])
legend("bottomleft",legend=1:count.max,col=colorscheme,bty="n",pch=19,cex = 0.75)
box()


##############################################################################
#Fossil Occurrences
##############################################################################
#fossils occurrence data
load("data/fossils.RData")
ex_fossils<-fossils
timebin <- .bincode(ex_fossils$early_age, c(0,1,2,4,5,15,20,21),TRUE,TRUE)
colorscheme <- colorRampPalette(c("lightblue","yellow","orange","red","darkred"))(7)
plot(counts.m["X",] + 1/2, counts.m["Y",] + 1/2, type="n", xlab="Longitude",ylab="")
map("world",fill=T, col="lightgray",add=T)
map("state",add=T)
points(ex_fossils$Longitude,ex_fossils$Latitude, col=colorscheme[timebin],pch=16,cex=1.5)
legend("bottomleft",legend=rev(c("Early Miocene", "Mid Miocene", "Late Miocene", "Early Pliocene","Late Pliocene","Early Pleistocene","Late Pleistocene")),col=colorscheme,pch=16,bty="n",cex = 0.65)
box()
#display.brewer.all(colorblindFriendly=TRUE)

#fossils taxon data
colorscheme <- colorRampPalette(c("lightblue","orange","red"))(3)
plot(counts.m["X",] + 1/2, counts.m["Y",] + 1/2, type="n", xlab="Longitude",ylab="")
map("world",fill=T, col="lightgray",add=T)
map("state",add=T)
points(ex_fossils$Longitude,ex_fossils$Latitude, col=colorscheme[ex_fossils$Taxon],pch=16,cex=1.5)
legend("bottomleft",legend=c("Sceloporus sp.","Sceloporus robustus","Sceloprous undulatus"),col=colorscheme,pch=16,bty="n",cex=0.75)
box()

