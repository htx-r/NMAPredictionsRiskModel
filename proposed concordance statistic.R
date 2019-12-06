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
                                      n.chains=2,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 10)

#parallel::stopCluster(cl)

p<-IPDNMRJAGSmodelLASSOconcordance$BUGSoutput$mean$logitp
expit<-function(x) {exp(x)/(1+exp(x))}
p<-expit(IPDNMRJAGSmodelLASSOconcordance$BUGSoutput$mean$logitp)
p<-as.data.frame(p)
Relapse2year<-RiskData$RELAPSE2year
Relapse2year<-as.data.frame(Relapse2year)
Data<-cbind(p,Relapse2year)
benefit<-as.data.frame(Data$V4-Data$V3)
FinalData<-cbind(Data,benefit)
Treatment<-RiskData$TRT01A
FinalData<-cbind(FinalData,Treatment)
FinalData<-FinalData[which(FinalData$Treatment==3  | FinalData$Treatment==4),]
FinalDataNat<-FinalData[which(FinalData$Treatment==3),]
FinalDataPla<-FinalData[which(FinalData$Treatment==4),]
FinalDataPla<-sample_n(FinalDataPla,577)

FinalDataPla<-FinalDataPla[order(FinalDataPla$`Data$V4 - Data$V3`),]

FinalDataNat<-FinalDataNat[order(FinalDataNat$`Data$V4 - Data$V3`),]




pred.ben.avg<-(FinalDataNat$`Data$V4 - Data$V3`+FinalDataPla$`Data$V4 - Data$V3`)/2
obs.ben<-as.numeric(FinalDataNat$Relapse2year)-as.numeric(FinalDataPla$Relapse2year)
# Matched patient pair data
x<-data.frame(FinalDataNat$`Data$V4 - Data$V3`,FinalDataNat$Relapse2year,FinalDataPla$`Data$V4 - Data$V3`,FinalDataPla$Relapse2year,pred.ben.avg,obs.ben)

cindex <- rcorr.cens(pred.ben.avg, obs.ben)
c.benefit <- cindex["C Index"][[1]]
c.benefit.se <- cindex["S.D."][[1]]/2	# The sd of the c-index is half the sd of Dxy

c.benefit
# [1] 0.6363636		# The c-for-benefit
c.benefit - 1.96*c.benefit.se
# [1] 0.3833374		# The 95% lower bound of the c-for-benefit
c.benefit + 1.96*c.benefit.se
# [1] 0.8893899

