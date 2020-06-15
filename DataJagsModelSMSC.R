#dummy variables for nr.relapses.2y.prior.study
SMSCdataC$nr.relapses.2y.prior.study1<-SMSCdataC$nr.relapses.2y.prior.study
SMSCdataC$nr.relapses.2y.prior.study2<-SMSCdataC$nr.relapses.2y.prior.study
SMSCdataC$nr.relapses.2y.prior.study1[which(SMSCdataC$nr.relapses.2y.prior.study1==2)]<-0
SMSCdataC$nr.relapses.2y.prior.study2[which(SMSCdataC$nr.relapses.2y.prior.study2==1)]<-0
SMSCdataC$nr.relapses.2y.prior.study2[which(SMSCdataC$nr.relapses.2y.prior.study2==2)]<-1




#give the data
jagsdataSMSC <- list(
  Nobservations=nrow(SMSCdataC),
  outcome=as.numeric(factor(SMSCdataC$relapse.2y.after.study))-1,
  npf=10,
  npid=length(unique(SMSCdataC$patient.id)),
  subj=as.integer(factor(SMSCdataC$patient.id)),
  age=SMSCdataC$age-mean(SMSCdataC$age),
  disease.duration=SMSCdataC$disease.duration-mean(SMSCdataC$disease.duration),
  edss=SMSCdataC$edss-mean(SMSCdataC$edss),
  nr.Gd.enhanced.lesions=SMSCdataC$nr.Gd.enhanced.lesions,
  nr.relapses.2y.prior.study1=SMSCdataC$nr.relapses.2y.prior.study1,
  nr.relapses.2y.prior.study2=SMSCdataC$nr.relapses.2y.prior.study2,
  months.since.last.relapse=SMSCdataC$months.since.last.relapse-mean(SMSCdataC$months.since.last.relapse),
  treatment.naive.prior.visit=SMSCdataC$treatment.naive.prior.visit ,
  gender=as.factor(SMSCdataC$gender),
  treatment.time.during.cycle.months=SMSCdataC$treatment.time.during.cycle.months-mean(SMSCdataC$treatment.time.during.cycle.months)

)
