####################################################################
################## Script to add columns arm and meanRisk ############
################## that are needed for the IPD NMR model ################

#mean of logit
RiskData$meanlogitRisk<-NA
RiskData$meanlogitRisk[RiskData$STUDYID==1]<-mean(RiskData$logitRisk[RiskData$STUDYID==1])
RiskData$meanlogitRisk[RiskData$STUDYID==2]<-mean(RiskData$logitRisk[RiskData$STUDYID==2])
RiskData$meanlogitRisk[RiskData$STUDYID==3]<-mean(RiskData$logitRisk[RiskData$STUDYID==3])
RiskData$meanlogitRisk[RiskData$STUDYID==4]<-mean(RiskData$logitRisk[RiskData$STUDYID==4])
#mean of risk
RiskData$meanRisk<-NA
RiskData$meanRisk[RiskData$STUDYID==1]<-mean(RiskData$Risk[RiskData$STUDYID==1])
RiskData$meanRisk[RiskData$STUDYID==2]<-mean(RiskData$Risk[RiskData$STUDYID==2])
RiskData$meanRisk[RiskData$STUDYID==3]<-mean(RiskData$Risk[RiskData$STUDYID==3])
RiskData$meanRisk[RiskData$STUDYID==4]<-mean(RiskData$Risk[RiskData$STUDYID==4])
#arm
RiskData$arm<-NA
RiskData$arm[RiskData$STUDYID==1 & RiskData$TRT01A==2]<-1
RiskData$arm[RiskData$STUDYID==1 & RiskData$TRT01A==5]<-2
RiskData$arm[RiskData$STUDYID==2 & RiskData$TRT01A==2]<-1
RiskData$arm[RiskData$STUDYID==2 & RiskData$TRT01A==3]<-2
RiskData$arm[RiskData$STUDYID==2 & RiskData$TRT01A==5]<-3
RiskData$arm[RiskData$STUDYID==3 & RiskData$TRT01A==4]<-1
RiskData$arm[RiskData$STUDYID==3 & RiskData$TRT01A==5]<-2
RiskData$arm[RiskData$STUDYID==4 & RiskData$TRT01A==1]<-1
RiskData$arm[RiskData$STUDYID==4 & RiskData$TRT01A==5]<-2
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
  Nstudies=4,
  Np=nrow(RiskData),
  studyid=as.numeric(RiskData$STUDYID),
  outcome=RiskData$RELAPSE2year,
  outcomeP=PlaceboArms$Relapse2year,
  NpPlacebo=nrow(PlaceboArms),
  treat= rbind(c(2,5,NA),c(2,3,5),c(4,5,NA),c(1,5,NA)),
  na=c(2,3,2,2),
  logitRisknew=logitRisknew,
  logitmeanRisknew=logitmeanRisknew,
  arm=RiskData$arm,
  Risk=RiskData$logitRisk,
  meanRisk=c(-0.5257,-0.6467,-0.5047,-0.6143), ##here is the mean of logit of risk
  nt=5,
  ref=5
)

