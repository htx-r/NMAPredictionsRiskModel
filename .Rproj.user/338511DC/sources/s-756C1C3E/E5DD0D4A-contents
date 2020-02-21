##############################################################################################################################
####################################### Script that make tha DATA proper fot JAGS IPD AD NMA model  ###########################
###################################################################################################################################

#### Aggregated data
MS <- read_excel("MS.xlsx")
MS <- as.data.frame(MS)

ADdata<-MS[which(MS$study=="Bornstein" | MS$study=="Johnson"),]



######DATA
## arm column
RiskData$arm<-NA
RiskData$arm[RiskData$STUDYID==1 & RiskData$TRT01A==1]<-1
RiskData$arm[RiskData$STUDYID==1 & RiskData$TRT01A==4]<-2
RiskData$arm[RiskData$STUDYID==2 & RiskData$TRT01A==1]<-1
RiskData$arm[RiskData$STUDYID==2 & RiskData$TRT01A==2]<-2
RiskData$arm[RiskData$STUDYID==2 & RiskData$TRT01A==4]<-3
RiskData$arm[RiskData$STUDYID==3 & RiskData$TRT01A==3]<-1
RiskData$arm[RiskData$STUDYID==3 & RiskData$TRT01A==4]<-2


jagsdataIPDADNMA<- list(
  N.IPD.studies=3,
  Np=sum(as.numeric(table(as.numeric(as.factor(RiskData$STUDYID))))),
  studyid=as.numeric(as.factor(RiskData$STUDYID)),
  outcome=as.numeric(RiskData$RELAPSE2year)-1,
  treat= rbind(c(1,4,NA),c(1,2,4),c(3,4,NA), c(2,4,NA),c(2,4,NA)),
  na=c(2,3,2,2,2),
  arm=RiskData$arm,
  ref=4,
  nt=4,
  N.AD.studies=2,
  outcome.ad=rbind(c(NA,NA,NA,NA),c(NA,NA,NA,NA),c(NA,NA,NA,NA),c(NA,11,NA,19),c(NA,89,NA,97)),
  n=rbind(c(NA,NA,NA,NA),c(NA,NA,NA,NA),c(NA,NA,NA,NA),c(NA,25,NA,25),c(NA,125,NA,126)),
  Risk=RiskData$logitRiskFabio

)
