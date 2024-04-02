#' @title plotBumpChart
#' @description
#' @usage plotBumpChart(sp_data_mean, which.biovars, path="")
#' @param sp_data_mean
#' @param which.biovars a vector with the biovars (see www.worldclim.org for a list of biovars).
#' @param path path to the directory where the results should be saved
#' @return bump chart of mean species-climate associations
#' @details
#' @return
#' @author A. Michelle Lawing, Alexandra F. C. Howard
#' @examples

#biovars=c(1,4,15)
#sp_data_mean<-sapply(4:(length(which.biovars)+3),function(x) tapply(occurrences[,x],occurrences$Species,mean))
#plotBumpChart(sp_data_mean,which.biovars=biovars)

plotBumpChart<-function(sp_data_mean, which.biovars, path = ""){
  require(ggplot2)
  require(reshape2)
  jpeg(paste(path,"bump_chart.jpg",sep=""),width=960,height=960,quality=100,pointsize=20)
  trans <- apply(sp_data_mean,2,rank)
  melted <- melt(trans)
  theme_set(theme_bw())
  zmargin<-theme(panel.spacing=unit(0,"lines"))
  b1 <- ggplot(melted, aes(Var2, value, group = Var1, color = Var1, label = Var1)) + geom_line() +
    geom_text(aes(label = as.factor(c(paste("S.",as.factor(unique(melted$Var1))),rep("",length(melted$Var2) - length(as.character(unique(melted$Var1)))))),
                  x = melted$Var2 - 0.01, size = 3, hjust = 1, fontface = "italic")) + theme_bw() +
    theme(legend.position = "none", panel.border = element_blank(), axis.ticks = element_blank()) +
    scale_x_discrete(breaks = levels(as.factor(melted$Var2)), labels = paste("Bio",which.biovars, sep = ""), expand = c(0.5,0.9)) +
    xlab(NULL) + scale_y_continuous(breaks = NULL) + ylab(NULL) +
    annotate("text", x = which.biovars, y = rep(max(trans) + 2, length(which.biovars)),
             label = as.integer(apply(sp_data_mean,2,max)), size = 3) +
    annotate("text", x = which.biovars, y = rep(min(trans) - 2, length(which.biovars)),
             label = as.integer(apply(sp_data_mean,2,min)), size = 3)
  b1
  dev.off()
}



