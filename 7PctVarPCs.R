#this code will go through the "sdev" outputs from the RunPCA code
# and will compute the percent variability from each component.
# Then it will pull out these %s for the top 10 pcs and output them

#Updated 9/9/2017 for 2 seasons: Nov-Apr & May-Oct


rm(list=ls())

setwd("C:/AbbyF/Oliver/PC_outputs_R_T/")

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


