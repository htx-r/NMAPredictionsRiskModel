####################################################################
################## Script to add columns arm and meanRisk ############
################## that are needed for the IPD NMR model ################


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

###data for jagsmodel with metarigression on logit of Risk for LASSO model
jagsdataIPDNMRLASSO <- list(
  Nstudies=3,
  Np=nrow(RiskData),
  studyid=as.numeric(RiskData$STUDYID),
  outcome=as.numeric(RiskData$RELAPSE2year),
  outcomeP=PlaceboArms$Relapse2year,
  NpPlacebo=nrow(PlaceboArms),
  treat= rbind(c(1,4,NA),c(1,2,4),c(3,4,NA)),
  na=c(2,3,2),
  logitRisknew=logitRisknew,
  logitmeanRisknew=logitmeanRisknew,
  arm=RiskData$arm,
  Risk=RiskData$logitRiskLASSO,
  meanRisk=c(-0.5335,-0.5832,-0.5091), ##here is the mean of logit of risk
  nt=4,
  ref=4
)

jagsdataIPDNMRFabio <- list(
  Nstudies=3,
  Np=nrow(RiskData),
  studyid=as.numeric(RiskData$STUDYID),
  outcome=RiskData$RELAPSE2year,
  outcomeP=PlaceboArms$Relapse2year,
  NpPlacebo=nrow(PlaceboArms),
  treat= rbind(c(1,4,NA),c(1,2,4),c(3,4,NA)),
  na=c(2,3,2),
  logitRisknew=logitRisknew,
  logitmeanRisknew=logitmeanRisknew,
  arm=RiskData$arm,
  Risk=RiskData$logitRiskFabio,
  meanRisk=c(-0.5450,-0.5757,-0.5374), ##here is the mean of logit of risk
  nt=4,
  ref=4
)

