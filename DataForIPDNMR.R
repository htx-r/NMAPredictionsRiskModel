####################################################################
################## Script to add columns arm and meanRisk ############
################## that are needed for the IPD NMR model ################

RiskData$meanlogitRisk<-NA
RiskData$meanlogitRisk[RiskData$STUDYID==1]<-mean(RiskData$logitRisk[RiskData$STUDYID==1])
RiskData$meanlogitRisk[RiskData$STUDYID==2]<-mean(RiskData$logitRisk[RiskData$STUDYID==2])
RiskData$meanlogitRisk[RiskData$STUDYID==3]<-mean(RiskData$logitRisk[RiskData$STUDYID==3])

RiskData$meanRisk<-NA
RiskData$meanRisk[RiskData$STUDYID==1]<-mean(RiskData$Risk[RiskData$STUDYID==1])
RiskData$meanRisk[RiskData$STUDYID==2]<-mean(RiskData$Risk[RiskData$STUDYID==2])
RiskData$meanRisk[RiskData$STUDYID==3]<-mean(RiskData$Risk[RiskData$STUDYID==3])

RiskData$arm<-NA
RiskData$arm[RiskData$STUDYID==1 & RiskData$TRT01A==1]<-1
RiskData$arm[RiskData$STUDYID==1 & RiskData$TRT01A==4]<-2
RiskData$arm[RiskData$STUDYID==2 & RiskData$TRT01A==1]<-1
RiskData$arm[RiskData$STUDYID==2 & RiskData$TRT01A==2]<-2
RiskData$arm[RiskData$STUDYID==2 & RiskData$TRT01A==4]<-3
RiskData$arm[RiskData$STUDYID==3 & RiskData$TRT01A==3]<-1
RiskData$arm[RiskData$STUDYID==3 & RiskData$TRT01A==4]<-2

###data for jagsmodel with metarigression on logit of Risk
jagsdataIPDNMR <- list(
  Nstudies=3,
  Np=sum(as.numeric(table(as.numeric(as.factor(RiskData$STUDYID))))),
  studyid=as.numeric(RiskData$STUDYID),
  outcome=RiskData$RELAPSE2year,
  treat= rbind(c(1,4,NA),c(1,2,4),c(3,4,NA)),
  na=c(2,3,2),
  arm=RiskData$arm,
  Risk=RiskData$logitRisk,
  meanRisk=c(-0.5360,-0.6501,-0.5110), ##here is the mean of logit of risk
  nt=4,
  ref=4
)
