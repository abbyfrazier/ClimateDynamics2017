#MLR with PC data and MEI, PDO & PNA indices AND CO2 DATA (in situ) FROM MLO: 1974-2012
#Using PCA of rainfall maps

#1920-2012
#....it will delete any that don't match, so really it's 1974-2012!


#Updated 9/9/2017 for 2 seasons: Nov-Apr & May-Oct


rm(list=ls())
library(car)

island_list<-c("ka_","oa_","ma_","bi_")
seas_list<-c("na_","mo_")

k=1 #counter of total number of isl x seas combos

#location of indices:
setwd("C:/AbbyF/Oliver/Indices/")

indices_na<-read.table("Indices_NA.txt",header=T)
indices_mo<-read.table("Indices_MO.txt",header=T)

co2_na<-read.table("co2_na.txt",header=T)
co2_mo<-read.table("co2_mo.txt",header=T)



#set up output file for Rsquared values of models, 4x2 rows (4 isl x 2 seas)
rsq<-data.frame(matrix(NA,nrow=8,ncol=4))
colnames(rsq)[1]<-"PC1"
colnames(rsq)[2]<-"PC2"
colnames(rsq)[3]<-"PC3"
colnames(rsq)[4]<-"PC4"

rsqA<-rsq  #make rsq file for additive models
pvalA<-rsq


#Loop through all islands and seasons, run regression on PCs 1-4 (Full - all interactions),
#  save anova tables and R-sq of each model
for (i in 1:4) {
  isl<-island_list[i]
  
  for (j in 1:2) {
    seas<-seas_list[j]
    
    #set which index based on season
    if(j==1){
      ind<-indices_na
      CO2<-co2_na$AvgCO2
      nyr=92
    }
    if(j==2){
      ind<-indices_mo
      CO2<-co2_mo$AvgCO2
      nyr=93
    }

    
    MEI<-ind[55:nyr,1]
    PDO<-ind[55:nyr,2]
    PNA<-ind[55:nyr,4]
    
    
    setwd("C:/AbbyF/Oliver/PC_outputs_R_T/new_ts/")
    pc_ts<-read.csv(paste(isl,seas,"timeseries.csv",sep=""))
    pc_ts<-pc_ts[55:nyr,] #****************************************to line up with CO2
    
    rownames(rsq)[k]<-paste(isl,seas,"rsq",sep="")
    rownames(rsqA)[k]<-paste(isl,seas,"rsq",sep="")
    
     
    #***************Additive models:
    #PC1:
    fit1a<-lm(pc_ts$PC1~MEI+PDO+PNA+CO2)
    rsqA[k,1]<-summary(fit1a)$r.squared
    #coefficients(fit1)
    a1a<-Anova(fit1a,type=2)
    f=summary(fit1a)$fstatistic
    pvalA[k,1]<-pf(f[1],f[2],f[3],lower.tail=FALSE)
    
    #PC2:
    fit2a<-lm(pc_ts$PC2~MEI+PDO+PNA+CO2)
    rsqA[k,2]<-summary(fit2a)$r.squared
    a2a<-Anova(fit2a,type=2)
    f=summary(fit2a)$fstatistic
    pvalA[k,2]<-pf(f[1],f[2],f[3],lower.tail=FALSE)
    
    #PC3:
    fit3a<-lm(pc_ts$PC3~MEI+PDO+PNA+CO2)
    rsqA[k,3]<-summary(fit3a)$r.squared
    a3a<-Anova(fit3a,type=2)
    f=summary(fit3a)$fstatistic
    pvalA[k,3]<-pf(f[1],f[2],f[3],lower.tail=FALSE)
    
    #PC4:
    fit4a<-lm(pc_ts$PC4~MEI+PDO+PNA+CO2)
    rsqA[k,4]<-summary(fit4a)$r.squared
    a4a<-Anova(fit4a,type=2)
    f=summary(fit4a)$fstatistic
    pvalA[k,4]<-pf(f[1],f[2],f[3],lower.tail=FALSE)
    
    #output anova tables:
    setwd("C:/AbbyF/Oliver/MLR/NEW_TS/co2/")
    write.table(a1a, file = paste(isl,seas,"AdditAnovaPC1.csv",sep=""), sep = ",", col.names = T)
    write.table(a2a, file = paste(isl,seas,"AdditAnovaPC2.csv",sep=""), sep = ",", col.names = T)
    write.table(a3a, file = paste(isl,seas,"AdditAnovaPC3.csv",sep=""), sep = ",", col.names = T)
    write.table(a4a, file = paste(isl,seas,"AdditAnovaPC4.csv",sep=""), sep = ",", col.names = T)
    
    k = k + 1
  }
  
}

write.table(rsqA, file = "RSq_AdditModels.csv", sep = ",", col.names = T)
write.table(pvalA, file = "PVals_AdditModels.csv", sep = ",", col.names = T)

