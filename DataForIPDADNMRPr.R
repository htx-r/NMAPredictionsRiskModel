####################################################################################################################################
########################################### DATA FOR IPD AND AD NMR with risk as prognostic factor ONLY ############################
####################################################################################################################################


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
  Risk=RiskData$logitRiskFabio,
  meanRisk=c(tapply(RiskData$logitRiskFabio, RiskData$STUDYID, summary)$`1`[4],tapply(RiskData$logitRiskFabio, RiskData$STUDYID, summary)$`2`[4],tapply(RiskData$logitRiskFabio, RiskData$STUDYID, summary)$`3`[4]) ##here is the mean of logit of risk


)
