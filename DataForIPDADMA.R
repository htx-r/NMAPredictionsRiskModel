



#### Aggregated data
MS <- read_excel("MS.xlsx")
MS <- as.data.frame(MS)

ADdata<-MS[which(MS$study=="Bornstein" | MS$study=="Johnson"),]



######DATA


RiskData$treatma<-NA
RiskData$treatma[RiskData$TRT01A==1]<-2
RiskData$treatma[RiskData$TRT01A==2]<-2
RiskData$treatma[RiskData$TRT01A==3]<-2
RiskData$treatma[RiskData$TRT01A==4]<-1


jagsdataIPDADMA<- list(
  N.IPD.studies=3,
  Np=sum(as.numeric(table(as.numeric(as.factor(RiskData$STUDYID))))),
  studyid=as.numeric(as.factor(RiskData$STUDYID)),
  treatma=RiskData$treatma,
  outcome=as.numeric(RiskData$RELAPSE2year)-1,
  N.AD.studies=2,
  outcome.ad=rbind(c(NA,NA),c(NA,NA),c(NA,NA),c(19,11),c(97,89)),
  n=rbind(c(NA,NA),c(NA,NA),c(NA,NA),c(25,25),c(125,126))
)

IPDADMAJAGSmodel <- jags.parallel(data = jagsdataIPDADMA ,inits=NULL,parameters.to.save = c('u','ORref',"mean","d"),model.file =modelIPDADMA,
                                       n.chains=2,n.iter = 10000,n.burnin = 100,DIC=F,n.thin = 1)
