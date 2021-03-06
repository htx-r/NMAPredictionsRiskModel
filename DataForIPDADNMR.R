##############################################################################################################################
#####################      Script that make tha DATA proper for all JAGS IPD AD models    ###########################
###################################################################################################################################



#### 1. LOAD the Aggregated data
MS <- read_excel("MS.xlsx")
MS <- as.data.frame(MS)
### keep studied of interest
ADdata<-MS[which(MS$study=="Bornstein" | MS$study=="Johnson"),]

## DATA for IPD AD NMA model
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
  Risk=RiskData$logitRiskPreSpecified

)

## DATA for IPD AD MA model
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


## DATA for IPD AD NMR model with risk as prognostic factor only
jagsdataIPDADNMRPr <- list(
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
  Risk=RiskData$logitRiskPreSpecified,
  meanRisk=c(tapply(RiskData$logitRiskPreSpecified, RiskData$STUDYID, summary)$`1`[4],tapply(RiskData$logitRiskPreSpecified, RiskData$STUDYID, summary)$`2`[4],tapply(RiskData$logitRiskPreSpecified, RiskData$STUDYID, summary)$`3`[4]) ##here is the mean of logit of risk


)

## DATA for IPD AD NMR model with risk as effect modifier only within
jagsdataIPDADNMREMwithin <- list(
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
  Risk=RiskData$logitRiskPreSpecified,
  meanRisk=c(tapply(RiskData$logitRiskPreSpecified, RiskData$STUDYID, summary)$`1`[4],tapply(RiskData$logitRiskPreSpecified, RiskData$STUDYID, summary)$`2`[4],tapply(RiskData$logitRiskPreSpecified, RiskData$STUDYID, summary)$`3`[4], -0.5, -0.6) ##here is the mean of logit of risk
)

  ## DATA for IPD AD NMR model (all included)
  jagsdataIPDADNMR <- list(
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
    Risk=RiskData$logitRiskPreSpecified,
    meanRisk=c(tapply(RiskData$logitRiskPreSpecified, RiskData$STUDYID, summary)$`1`[4],tapply(RiskData$logitRiskPreSpecified, RiskData$STUDYID, summary)$`2`[4],tapply(RiskData$logitRiskPreSpecified, RiskData$STUDYID, summary)$`3`[4], -0.5, -0.6)##here is the mean of logit of risk


  )




