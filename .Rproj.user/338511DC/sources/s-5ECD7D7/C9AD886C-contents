jagsdataSMSC <- list(
  Nobservations=nrow(SMSCdataC),
  outcome=SMSCdataC$outcome,
  npf=9,
  npid=length(unique(SMSCdataC$patient.id)),
  #subj=sort(as.integer((SMSCdataC$NumericID))),
  subj=as.integer(factor(SMSCdataC$patient.id)),
  age=SMSCdataC$age,
  disease.duration=SMSCdataC$disease.duration,
  edss=SMSCdataC$edss,
  nr.Gd.enhanced.lesions=SMSCdataC$nr.Gd.enhanced.lesions,
  nr.relapses.2y.prior.study=SMSCdataC$nr.relapses.2y.prior.study,
  months.since.last.relapse=SMSCdataC$months.since.last.relapse,
  treatment.naive.prior.visit=SMSCdataC$treatment.naive.prior.visit ,
  gender=as.factor(SMSCdataC$gender),
  treatment.time.during.cycle.months=SMSCdataC$treatment.time.during.cycle.months

)

SMSCjagsResults <- jags.parallel(data = jagsdataSMSC ,inits=NULL,parameters.to.save = c('b'),model.file = jagsmodelSMSC,
                                        n.chains=2,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 10)

