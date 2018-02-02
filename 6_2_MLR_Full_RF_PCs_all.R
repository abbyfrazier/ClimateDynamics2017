#MLR with PC data and MEI, PDO & PNA indices AND TREND TERM
# OUTPUT *FULL* ANOVA MODELS (all interactions)
#Using PCA of rainfall maps (Grid EOF)
#All 3 Trend terms: Global temperature, CO2, and linear trend
#2 seasons: Nov-Apr & May-Oct


#1920-2012
#....PNA makes it 1950-2012


rm(list=ls())
library(car)

island_list<-c("ka_","oa_","ma_","bi_")
seas_list<-c("na_","mo_")

k=1 #counter of total number of isl x seas combos

#location of indices:
setwd("C:/AbbyF/Oliver/Indices/")

indices_na<-read.table("Indices_NA.txt",header=T)
indices_mo<-read.table("Indices_MO.txt",header=T)

#GLOBAL TEMPERATURE as trend:
temp_na<-read.table("GlobTemp_NA.txt",header=T)
temp_mo<-read.table("GlobTemp_MO.txt",header=T)

#CO2 as trend:
co2_na<-read.table("co2_na.txt",header=T)
co2_mo<-read.table("co2_mo.txt",header=T)

#Linear trend as trend term:
lineartrends<-read.table("LinearTrends.txt",header=T)
trnd_na<-lineartrends[1:92,2:5]

#Which trend term?  Un-comment the one you want:
#Update the "if" statement too
#   and for CO2, need to un-comment the pc_ts lines to shorten it (line 88-91, 96)

tr_term<-"temp/"
#tr_term<-"LinearTrend/trnd05/"
#tr_term<-"co2/"



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
      tr<-temp_na$TempC
      #tr<-trnd_na$trnd05
      #tr<-co2_na$AvgCO2
      nyr=92
    }
    if(j==2){
      ind<-indices_mo
      tr<-temp_mo$TempC
      #tr<-lineartrends$trnd05
      #tr<-co2_mo$AvgCO2
      nyr=93
    }
    
    #comment these out for CO2:
    MEI<-ind[,1]
    PDO<-ind[,2]
    PNA<-ind[,4]
    
    #FOR CO2 ONLY:
    #MEI<-ind[55:nyr,1]
    #PDO<-ind[55:nyr,2]
    #PNA<-ind[55:nyr,4]
    
    
    setwd("C:/AbbyF/Oliver/PC_outputs_R_T/new_ts/")
    pc_ts<-read.csv(paste(isl,seas,"timeseries.csv",sep=""))
    #pc_ts<-pc_ts[55:nyr,] #****************************************to line up with CO2
    
    rownames(rsq)[k]<-paste(isl,seas,"rsq",sep="")
    rownames(rsqA)[k]<-paste(isl,seas,"rsq",sep="")
    rownames(pvalA)[k]<-paste(isl,seas,"pval",sep="")
    
    
    #***************Full models - all interactions:
    #PC1:
    fit1a<-lm(pc_ts$PC1~MEI*PDO*PNA*tr)
    rsqA[k,1]<-summary(fit1a)$r.squared
    #coefficients(fit1)
    a1a<-Anova(fit1a,type=2)
    f=summary(fit1a)$fstatistic
    pvalA[k,1]<-pf(f[1],f[2],f[3],lower.tail=FALSE)
    
    #PC2:
    fit2a<-lm(pc_ts$PC2~MEI*PDO*PNA*tr)
    rsqA[k,2]<-summary(fit2a)$r.squared
    a2a<-Anova(fit2a,type=2)
    f=summary(fit2a)$fstatistic
    pvalA[k,2]<-pf(f[1],f[2],f[3],lower.tail=FALSE)
    
    #PC3:
    fit3a<-lm(pc_ts$PC3~MEI*PDO*PNA*tr)
    rsqA[k,3]<-summary(fit3a)$r.squared
    a3a<-Anova(fit3a,type=2)
    f=summary(fit3a)$fstatistic
    pvalA[k,3]<-pf(f[1],f[2],f[3],lower.tail=FALSE)
    
    #PC4:
    fit4a<-lm(pc_ts$PC4~MEI*PDO*PNA*tr)
    rsqA[k,4]<-summary(fit4a)$r.squared
    a4a<-Anova(fit4a,type=2)
    f=summary(fit4a)$fstatistic
    pvalA[k,4]<-pf(f[1],f[2],f[3],lower.tail=FALSE)
    
    #output anova tables:
    setwd(paste("C:/AbbyF/Oliver/MLR/NEW_TS/",tr_term,sep=""))
    write.table(a1a, file = paste(isl,seas,"FullAnovaPC1.csv",sep=""), sep = ",", col.names = T)
    write.table(a2a, file = paste(isl,seas,"FullAnovaPC2.csv",sep=""), sep = ",", col.names = T)
    write.table(a3a, file = paste(isl,seas,"FullAnovaPC3.csv",sep=""), sep = ",", col.names = T)
    write.table(a4a, file = paste(isl,seas,"FullAnovaPC4.csv",sep=""), sep = ",", col.names = T)
    
    k = k + 1
  }
  
}

write.table(rsqA, file = "RSq_FullModels.csv", sep = ",", col.names = T)
write.table(pvalA, file = "PVals_FullModels.csv", sep = ",", col.names = T)
