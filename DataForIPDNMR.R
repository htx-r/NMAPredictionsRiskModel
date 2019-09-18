####################################################################
################## Script to add columns arm and meanRisk ############
################## that are needed for the IPD NMR model ################

#mean of logit
RiskData$meanlogitRisk<-NA
RiskData$meanlogitRisk[RiskData$STUDYID==1]<-mean(RiskData$logitRisk[RiskData$STUDYID==1])
RiskData$meanlogitRisk[RiskData$STUDYID==2]<-mean(RiskData$logitRisk[RiskData$STUDYID==2])
RiskData$meanlogitRisk[RiskData$STUDYID==3]<-mean(RiskData$logitRisk[RiskData$STUDYID==3])
#mean of risk
RiskData$meanRisk<-NA
RiskData$meanRisk[RiskData$STUDYID==1]<-mean(RiskData$Risk[RiskData$STUDYID==1])
RiskData$meanRisk[RiskData$STUDYID==2]<-mean(RiskData$Risk[RiskData$STUDYID==2])
RiskData$meanRisk[RiskData$STUDYID==3]<-mean(RiskData$Risk[RiskData$STUDYID==3])
#arm
RiskData$arm<-NA
RiskData$arm[RiskData$STUDYID==1 & RiskData$TRT01A==1]<-1
RiskData$arm[RiskData$STUDYID==1 & RiskData$TRT01A==4]<-2
RiskData$arm[RiskData$STUDYID==2 & RiskData$TRT01A==1]<-1
RiskData$arm[RiskData$STUDYID==2 & RiskData$TRT01A==2]<-2
RiskData$arm[RiskData$STUDYID==2 & RiskData$TRT01A==4]<-3
RiskData$arm[RiskData$STUDYID==3 & RiskData$TRT01A==3]<-1
RiskData$arm[RiskData$STUDYID==3 & RiskData$TRT01A==4]<-2


##new Risk & logit Risk
Risknew<-seq(0.01,0.99,0.01)
Risknew<-as.data.frame(Risknew)
logit<-function(x) {log(x/(1-x))}
expit<-function(x) {exp(x)/(1+exp(x))}

logitRisknew<-NA
logitRisknew<-as.data.frame(logitRisknew)
for (i in 1:99) {
  logitRisknew[i,1]<-logit(Risknew[i,1])
}


logitRisknew<-as.data.frame(logitRisknew)
logitmeanRisknew<-mean(logitRisknew[,1])

###data for jagsmodel with metarigression on logit of Risk
jagsdataIPDNMR <- list(
  Nstudies=3,
  Np=sum(as.numeric(table(as.numeric(as.factor(RiskData$STUDYID))))),
  studyid=as.numeric(RiskData$STUDYID),
  outcome=RiskData$RELAPSE2year,
  outcomeP=PlaceboArms$Relapse2year,
  NpPlacebo=nrow(PlaceboArms),
  treat= rbind(c(1,4,NA),c(1,2,4),c(3,4,NA)),
  na=c(2,3,2),
  logitRisknew=logitRisknew,
  logitmeanRisknew=logitmeanRisknew,
  arm=RiskData$arm,
  Risk=RiskData$logitRisk,
  meanRisk=c(-0.5360,-0.6501,-0.5110), ##here is the mean of logit of risk
  nt=4,
  ref=4
)
