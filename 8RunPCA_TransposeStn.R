#EOF analysis on STATION data, calculate % variability from top 10 PCs

#Updated 9/9/2017 for 2 seasons: Nov-Apr & May-Oct



#**Nov-Apr
#Removed stations if they had >5 -9999s

rm(list=ls())

setwd("C:/AbbyF/Oliver/StationPCA/")  #where input data is stored

#make sure input files don't have negative/empty cells
#open 6-month seasonal rainfall files (inches)
bi_na<-read.table("bi_na.txt", header=T)
ma_na<-read.table("ma_na.txt", header=T)
ka_na<-read.table("ka_na.txt", header=T)
oa_na<-read.table("oa_na.txt", header=T)

bi_na[bi_na==-9999]<-NA
ma_na[ma_na==-9999]<-NA
ka_na[ka_na==-9999]<-NA
oa_na[oa_na==-9999]<-NA


#transpose:
bi_na_T<-t(bi_na)
ma_na_T<-t(ma_na)
ka_na_T<-t(ka_na)
oa_na_T<-t(oa_na)


#run PCA (on z-scores, using the "center" and "scale" functions), save outputs
bi_na.pr<-prcomp(na.omit(bi_na_T),center=T,scale=T)
ma_na.pr<-prcomp(na.omit(ma_na_T),center=T,scale=T)
ka_na.pr<-prcomp(na.omit(ka_na_T),center=T,scale=T)
oa_na.pr<-prcomp(na.omit(oa_na_T),center=T,scale=T)

setwd("C:/AbbyF/Oliver/StationPCA/output/")

#matrix mult to get PC series
bi_na_matrix <- bi_na_T%*%bi_na.pr$rotation
ma_na_matrix <- ma_na_T%*%ma_na.pr$rotation
ka_na_matrix <- ka_na_T%*%ka_na.pr$rotation
oa_na_matrix <- oa_na_T%*%oa_na.pr$rotation

write.table(ma_na_matrix, file = "ma_na_timeseries.csv", sep = ",")
write.table(ma_na.pr$rotation, file = "ma_na_spatial_pr.csv", sep = ",", col.names = T)
write.table(ma_na.pr$sdev, file = "ma_na_sdev_pr.csv", sep = ",")
#rotation is now spatial series

write.table(bi_na_matrix, file = "bi_na_timeseries.csv", sep = ",")
write.table(bi_na.pr$rotation, file = "bi_na_spatial_pr.csv", sep = ",", col.names = T)
write.table(bi_na.pr$sdev, file = "bi_na_sdev_pr.csv", sep = ",")

write.table(ka_na_matrix, file = "ka_na_timeseries.csv", sep = ",")
write.table(ka_na.pr$rotation, file = "ka_na_spatial_pr.csv", sep = ",", col.names = T)
write.table(ka_na.pr$sdev, file = "ka_na_sdev_pr.csv", sep = ",")

write.table(oa_na_matrix, file = "oa_na_timeseries.csv", sep = ",")
write.table(oa_na.pr$rotation, file = "oa_na_spatial_pr.csv", sep = ",", col.names = T)
write.table(oa_na.pr$sdev, file = "oa_na_sdev_pr.csv", sep = ",")




#**May-Oct
#Removed stations if they had >5 -9999s

rm(list=ls())

setwd("C:/AbbyF/Oliver/StationPCA/")  #where input data is stored

#make sure input files don't have negative/empty cells
#open 3-month seasonal rainfall files (inches)
bi_mo<-read.table("bi_mo.txt", header=T)
ma_mo<-read.table("ma_mo.txt", header=T)
ka_mo<-read.table("ka_mo.txt", header=T)
oa_mo<-read.table("oa_mo.txt", header=T)

bi_mo[bi_mo==-9999]<-NA
ma_mo[ma_mo==-9999]<-NA
ka_mo[ka_mo==-9999]<-NA
oa_mo[oa_mo==-9999]<-NA


#transpose:
bi_mo_T<-t(bi_mo)
ma_mo_T<-t(ma_mo)
ka_mo_T<-t(ka_mo)
oa_mo_T<-t(oa_mo)


#run PCA (on z-scores, using the "center" and "scale" functions), save outputs
bi_mo.pr<-prcomp(na.omit(bi_mo_T),center=T,scale=T)
ma_mo.pr<-prcomp(na.omit(ma_mo_T),center=T,scale=T)
ka_mo.pr<-prcomp(na.omit(ka_mo_T),center=T,scale=T)
oa_mo.pr<-prcomp(na.omit(oa_mo_T),center=T,scale=T)

setwd("C:/AbbyF/Oliver/StationPCA/output/")

#matrix mult to get PC series
bi_mo_matrix <- bi_mo_T%*%bi_mo.pr$rotation
ma_mo_matrix <- ma_mo_T%*%ma_mo.pr$rotation
ka_mo_matrix <- ka_mo_T%*%ka_mo.pr$rotation
oa_mo_matrix <- oa_mo_T%*%oa_mo.pr$rotation

write.table(ma_mo_matrix, file = "ma_mo_timeseries.csv", sep = ",")
write.table(ma_mo.pr$rotation, file = "ma_mo_spatial_pr.csv", sep = ",", col.names = T)
write.table(ma_mo.pr$sdev, file = "ma_mo_sdev_pr.csv", sep = ",")
#rotation is now spatial series

write.table(bi_mo_matrix, file = "bi_mo_timeseries.csv", sep = ",")
write.table(bi_mo.pr$rotation, file = "bi_mo_spatial_pr.csv", sep = ",", col.names = T)
write.table(bi_mo.pr$sdev, file = "bi_mo_sdev_pr.csv", sep = ",")

write.table(ka_mo_matrix, file = "ka_mo_timeseries.csv", sep = ",")
write.table(ka_mo.pr$rotation, file = "ka_mo_spatial_pr.csv", sep = ",", col.names = T)
write.table(ka_mo.pr$sdev, file = "ka_mo_sdev_pr.csv", sep = ",")

write.table(oa_mo_matrix, file = "oa_mo_timeseries.csv", sep = ",")
write.table(oa_mo.pr$rotation, file = "oa_mo_spatial_pr.csv", sep = ",", col.names = T)
write.table(oa_mo.pr$sdev, file = "oa_mo_sdev_pr.csv", sep = ",")




#NOW CALC %VAR for top 10 PCs:

rm(list=ls())

setwd("C:/AbbyF/Oliver/StationPCA/output/")

island_list<-c("ka_","oa_","ma_","bi_")
seas_list<-c("na_","mo_")

top10<-data.frame(matrix(NA,nrow=10,ncol=2))

for (i in 1:4) {
  isl<-island_list[i]
  
  for (j in 1:2) {
    seas<-seas_list[j]
    
    sdev<-read.csv(paste(isl,seas,"sdev_pr.csv",sep=""))
    sdev$x2<-sdev$x^2
    sumvar<-sum(sdev$x2)
    sdev$pctvar<-sdev$x2/sumvar
    
    top10[j]<-data.frame(sdev[1:10,3])
    colnames(top10)[j]<-paste(isl,seas,"pctvar",sep="")
    
  }
  write.table(top10, file = paste(isl,"top10PCs.csv",sep=""), sep = ",", col.names = T)
}




