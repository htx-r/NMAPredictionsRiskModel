jagsdataIPDNMRLASSOconcordance <- list(
  Nstudies=3,
  Np=nrow(RiskData),
  studyid=as.numeric(RiskData$STUDYID),
  outcome=as.numeric(RiskData$RELAPSE2year)-1,
  outcomeP=PlaceboArms$Relapse2year,
  NpPlacebo=nrow(PlaceboArms),
  treat= rbind(c(1,4,NA),c(1,2,4),c(3,4,NA)),
  na=c(2,3,2),
  logitRisknew=as.data.frame(RiskData$logitRiskLASSO, ncol=1),
  logitmeanRisknew=mean(RiskData$logitRiskLASSO),
  Nnew=nrow(RiskData),
  arm=RiskData$arm,
  Risk=RiskData$logitRiskLASSO,
  meanRisk=c(tapply(RiskData$logitRiskLASSO, RiskData$STUDYID, summary)$`1`[4],tapply(RiskData$logitRiskLASSO, RiskData$STUDYID, summary)$`2`[4],tapply(RiskData$logitRiskLASSO, RiskData$STUDYID, summary)$`3`[4]), ##here is the mean of logit of risk
  nt=4,
  ref=4
)
#n <- parallel::detectCores()/2 # experiment!
#cl <- parallel::makeCluster(n)
#doParallel::registerDoParallel(cl)
IPDNMRJAGSmodelLASSOconcordance <- jags.parallel(data = jagsdataIPDNMRLASSOconcordance,inits=NULL,parameters.to.save = c('logitp'),model.file = modelIPDNMR,
                                      n.chains=2,n.iter = 100000,n.burnin = 1000,DIC=F,n.thin = 10)

#parallel::stopCluster(cl)



