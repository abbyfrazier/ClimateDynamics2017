#EOF analysis on raster data
#Updated 9/9/2017 for 2 seasons: Nov-Apr & May-Oct

rm(list=ls())

island_list<-c("ka_","oa_","ma_","bi_")
seas_list<-c("na_","mo_")


for (i in 1:4) {
  isl<-island_list[i]
  
  setwd("C:/AbbyF/Oliver/raster_txtfiles/")
  coords<-read.table(paste(isl,"coords.txt",sep=""), header=T)
  
  for (j in 1:2) {
    seas<-seas_list[j]
    
    setwd("C:/AbbyF/Oliver/raster_txtfiles/")
    rasters<-read.table(paste(isl,seas,"rasters_all.txt",sep=""),header=T)
    rasterT<-t(rasters)
    raster.pr<-prcomp(rasterT,center=T,scale=T)
    
    setwd("C:/AbbyF/Oliver/PC_outputs_R_T/new_ts/")
    prmatrix <- rasterT%*%raster.pr$rotation  #NOT USING THIS ANYMORE (as of 9/12/2017)
    spatial<-data.frame(cbind(coords,raster.pr$rotation))
    #rotation is now spatial series - tack on coordinates so this file can be mapped
    
    #grab time series for top 4 PCs:
    PC1 <- raster.pr$x[,'PC1']
    PC2 <- raster.pr$x[,'PC2']
    PC3 <- raster.pr$x[,'PC3']
    PC4 <- raster.pr$x[,'PC4']
    tsPC14<-data.frame(cbind(PC1,PC2,PC3,PC4))
    
    write.table(tsPC14, file = paste(isl,seas,"timeseries.csv",sep=""), sep = ",")
    write.table(spatial[,1:12], file = paste(isl,seas,"spatial_pr.csv",sep=""), sep = ",", col.names = T,row.names=F)
    write.table(raster.pr$sdev, file = paste(isl,seas,"sdev_pr.csv",sep=""), sep = ",")

  }
  
}





#time series:
library(xts)
#grab directly from raster.pr$x instead of using time series format:
#set end year based on season (Nov-Apr ends in 2011) - for ts command only
#if(j==1){
#  endyr=2011
#}
#if(j==2){
#  endyr=2012
#}

PC1 <- ts(raster.pr$x[,'PC1'],start=c(1920,1),end=c(endyr,1),frequency = 1)
PC2 <- ts(raster.pr$x[,'PC2'],start=c(1920,1),end=c(endyr,1),frequency = 1)
PC3 <- ts(raster.pr$x[,'PC3'],start=c(1920,1),end=c(endyr,1),frequency = 1)
PC4 <- ts(raster.pr$x[,'PC4'],start=c(1920,1),end=c(endyr,1),frequency = 1)
tsPC14<-data.frame(cbind(PC1,PC2,PC3,PC4))



##EXTRA 9/12/2017
#Screeplots:
std_dev <- raster.pr$sdev
pr_var <- std_dev^2
pr_var[1:10]
prop_varex <- pr_var/sum(pr_var)
prop_varex[1:10]

plot(prop_varex, xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     type = "b")
plot(cumsum(prop_varex), xlab = "Principal Component",
     ylab = "Cumulative Proportion of Variance Explained",
     type = "b")


##
#Maybe I'm not calculating time series right?
#prmatrix <- rasterT%*%raster.pr$rotation
#write.table(prmatrix, file = paste(isl,seas,"timeseries.csv",sep=""), sep = ",")

#https://stackoverflow.com/questions/41022927/principal-component-analysis-pca-of-time-series-data-spatial-and-temporal-pat
library(xts)
#for Nov-Apr stop in 2011
PC1 <- ts(raster.pr$x[,'PC1'],start=c(1920,1),end=c(2011,1),frequency = 1)
PC2 <- ts(raster.pr$x[,'PC2'],start=c(1920,1),end=c(2011,1),frequency = 1)
PC3 <- ts(raster.pr$x[,'PC3'],start=c(1920,1),end=c(2011,1),frequency = 1)
PC4 <- ts(raster.pr$x[,'PC4'],start=c(1920,1),end=c(2011,1),frequency = 1)
tsPC14<-data.frame(cbind(PC1,PC2,PC3,PC4))

plot(as.xts(PC1),major.format = "%Y-%b", type = 'l', main = "PC")
lines(as.xts(PC2),col='blue',type="l") # the blue one is PC2
cor(PC1,PC2)

##
library(factoextra) #http://www.sthda.com/english/wiki/principal-component-analysis-in-r-prcomp-vs-princomp-r-software-and-data-mining
res.pca<-prcomp(rasterT,center=T,scale=T)
# Eigenvalues
eig <- (res.pca$sdev)^2
# Variances in percentage
variance <- eig*100/sum(eig)
# Cumulative variances
cumvar <- cumsum(variance)
eig.decathlon2.active <- data.frame(eig = eig, variance = variance, cumvariance = cumvar)
head(eig.decathlon2.active)
summary(res.pca)
eig.val <- get_eigenvalue(res.pca)
head(eig.val)

fviz_screeplot(res.pca, ncp=10)
fviz_screeplot(res.pca, ncp=10, choice="eigenvalue")

var <- get_pca_var(res.pca)
var

head(var$coord[, 1:4])
# Helper function : 
# Correlation between variables and principal components
var_cor_func <- function(var.loadings, comp.sdev){
  var.loadings*comp.sdev
}
# Variable correlation/coordinates
loadings <- res.pca$rotation
sdev <- res.pca$sdev
var.coord <- var.cor <- t(apply(loadings, 1, var_cor_func, sdev))
head(var.coord[, 1:4])


# Plot the correlation circle
a <- seq(0, 2*pi, length = 100)
plot( cos(a), sin(a), type = 'l', col="gray",
      xlab = "PC1",  ylab = "PC2")
abline(h = 0, v = 0, lty = 2)
# Add active variables
arrows(0, 0, var.coord[, 1], var.coord[, 2], 
       length = 0.1, angle = 15, code = 2)
# Add labels
text(var.coord, labels=rownames(var.coord), cex = 1, adj=1)

#another way to plot the circle of variables:
fviz_pca_var(res.pca)


