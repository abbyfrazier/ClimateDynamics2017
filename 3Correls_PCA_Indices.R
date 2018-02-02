#calculate correlations between top PCs and indices of teleconnections
#Updated 9/9/2017 for 2 seasons: Nov-Apr & May-Oct


rm(list=ls())

island_list<-c("ka_","oa_","ma_","bi_")
seas_list<-c("na_","mo_")


#location of indices:
setwd("C:/AbbyF/Oliver/Indices/")

indices_na<-read.table("Indices_NA.txt",header=T)
indices_mo<-read.table("Indices_MO.txt",header=T)


#get correlations between indices and output as csv:
cor_na<-as.matrix(round(cor(indices_na,use="pairwise.complete.obs"),digits=3))
cor_mo<-as.matrix(round(cor(indices_mo,use="pairwise.complete.obs"),digits=3))

write.table(cor_na, file = "index_correls_na.csv", sep = ",", col.names = T)
write.table(cor_mo, file = "index_correls_mo.csv", sep = ",", col.names = T)


#calculate significance values for correlations!:
library("psych")  #will compute cor.test as a matrix:
cor_na.p<-corr.test(indices_na,adjust="none")
cor_mo.p<-corr.test(indices_mo,adjust="none")

write.table(round(cor_na.p$p,4), file = "index_correls_na_pval.csv", sep = ",", col.names = T)
write.table(round(cor_mo.p$p,4), file = "index_correls_mo_pval.csv", sep = ",", col.names = T)


library(car)
attach(indices_na)
scatterplotMatrix(~MEI_NA+PDO_NA+PNA_NA,main="Scatter plot matrix NA")
detach(indices_na)

attach(indices_mo)
scatterplotMatrix(~MEI_MO+PDO_MO+PNA_MO,main="Scatter plot matrix MO")
detach(indices_mo)



#now get PCs and get correlations between those and indices:

for (i in 1:4) {
  isl<-island_list[i]
  
  seas<-seas_list[1]
  setwd("C:/AbbyF/Oliver/PC_outputs_R_T/new_ts/")
  pc_ts<-read.csv(paste(isl,seas,"timeseries.csv",sep=""))
  PCnacor<-cor(pc_ts[,1:4],indices_na,use="pairwise.complete.obs")
  PCnacor.p<-corr.test(pc_ts[,1:4],indices_na,adjust="none")
  setwd("C:/AbbyF/Oliver/Correlations/new_ts/")
  write.table(PCnacor, file = paste(isl,seas,"correls.csv",sep=""), sep = ",", col.names = T)
  write.table(round(PCnacor.p$p,4), file = paste(isl,seas,"correls_pval.csv",sep=""), sep = ",", col.names = T)

  seas<-seas_list[2]
  setwd("C:/AbbyF/Oliver/PC_outputs_R_T/new_ts/")
  pc_ts<-read.csv(paste(isl,seas,"timeseries.csv",sep=""))
  PCmocor<-cor(pc_ts[,1:4],indices_mo,use="pairwise.complete.obs")
  PCmocor.p<-corr.test(pc_ts[,1:4],indices_mo,adjust="none")
  setwd("C:/AbbyF/Oliver/Correlations/new_ts/")
  write.table(PCmocor, file = paste(isl,seas,"correls.csv",sep=""), sep = ",", col.names = T)
  write.table(round(PCmocor.p$p,4), file = paste(isl,seas,"correls_pval.csv",sep=""), sep = ",", col.names = T)
  
  
}


#***Slopes of PCs...need to standardize units!
