#get correlations between each pixel and the indices, tack on coordinates, then map them in GIS
#Updated 9/9/2017 for 2 seasons: Nov-Apr & May-Oct

#Added p-value calculation
#takes 12-hours for one big island season (with 4 indices), less time with 3 indices


rm(list=ls()) 
library(psych) # for corr.test and corr.p

island_list<-c("ka_","oa_","ma_","bi_")
seas_list<-c("na_","mo_")

setwd("C:/AbbyF/Oliver/Indices/")

indices_na<-read.table("Indices_NA.txt",header=T)
indices_mo<-read.table("Indices_MO.txt",header=T)


for (i in 1:4) {
  isl<-island_list[i]
  
  setwd("C:/AbbyF/Oliver/raster_txtfiles/")
  coords<-read.table(paste(isl,"coords.txt",sep=""), header=T)
  
  seas<-seas_list[1]
  #open file
  setwd("C:/AbbyF/Oliver/raster_txtfiles/")  #where input data is stored
  rasters<-read.table(paste(isl,seas,"rasters_all.txt",sep=""), header=T)
  #transpose:
  rastersT<-t(rasters)
  rastersT.df<-as.data.frame(rastersT)  #NEED THIS STEP FOR corr.test TO WORK (p-values)
  #correlation (UPDATE INDEX FOR SEASON):
  cor.rastersT<-as.matrix(round(cor(rastersT,indices_na,use="pairwise.complete.obs"),digits=3))
  cor.rastersT.p<-corr.test(rastersT.df,indices_na,adjust="none")

  mapcor<-data.frame(cbind(coords,cor.rastersT))
  mapcor.p<-data.frame(cbind(coords,cor.rastersT.p$p))
  setwd("C:/AbbyF/Oliver/Correlations/")
  write.table(mapcor, file = paste(isl,seas,"RasterIndexCorrels.csv",sep=""), sep = ",", col.names = T,row.names=F)
  write.table(round(mapcor.p,8), file = paste(isl,seas,"RasterIndexCorrels_pval.csv",sep=""), sep = ",", col.names = T)
  
  seas<-seas_list[2]
  setwd("C:/AbbyF/Oliver/raster_txtfiles/")  
  rasters<-read.table(paste(isl,seas,"rasters_all.txt",sep=""), header=T)
  rastersT<-t(rasters)
  rastersT.df<-as.data.frame(rastersT)  #NEED THIS STEP FOR corr.test TO WORK (p-values)
  cor.rastersT<-as.matrix(round(cor(rastersT,indices_mo,use="pairwise.complete.obs"),digits=3))
  cor.rastersT.p<-corr.test(rastersT.df,indices_mo,adjust="none")
  mapcor<-data.frame(cbind(coords,cor.rastersT))
  mapcor.p<-data.frame(cbind(coords,cor.rastersT.p$p))
  setwd("C:/AbbyF/Oliver/Correlations/")
  write.table(mapcor, file = paste(isl,seas,"RasterIndexCorrels.csv",sep=""), sep = ",", col.names = T,row.names=F)
  write.table(round(mapcor.p,8), file = paste(isl,seas,"RasterIndexCorrels_pval.csv",sep=""), sep = ",", col.names = T)
  
  
}