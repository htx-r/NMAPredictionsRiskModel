
##################################################################################################################
############################# CHECK RESULTS WITH IPD DATA ONLY for Risk as EM within###############################################
####################################################################################################################



jagsdataIPDNMRFabio <- list(
  Nstudies=3,
  Np=nrow(RiskData),
  studyid=as.numeric(RiskData$STUDYID),
  outcome=as.numeric(RiskData$RELAPSE2year)-1,
  outcomeP=PlaceboArms$Relapse2year,
  NpPlacebo=nrow(PlaceboArms),
  treat= rbind(c(1,4,NA),c(1,2,4),c(3,4,NA)),
  na=c(2,3,2),
  logitRisknew=logitRisknew,
  logitmeanRisknew=-0.545,
  Nnew=99,
  arm=RiskData$arm,
  Risk=RiskData$logitRiskFabio,
  meanRisk=c(tapply(RiskData$logitRiskFabio, RiskData$STUDYID, summary)$`1`[4],tapply(RiskData$logitRiskFabio, RiskData$STUDYID, summary)$`2`[4],tapply(RiskData$logitRiskFabio, RiskData$STUDYID, summary)$`3`[4]), ##here is the mean of logit of risk
  ##here is the mean of logit of risk
  nt=4,
  ref=4
)

IPDNMRJAGSmodelFabio <- jags.parallel(data = jagsdataIPDNMRFabio,inits=NULL,parameters.to.save = c('gamma.w','gamma', 'ORref','delta','u'),model.file = modelIPDNMRWithinonly,
                                      n.chains=2,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 10)

IPDNMRJAGSmodelFabio
