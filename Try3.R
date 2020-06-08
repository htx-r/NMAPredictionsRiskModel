########################## Only sigma^2 in the diagonals ##############################################


jagsmodelSMSC2 <-function()
{


  for( i in 1:Nobservations)
  {
    outcome[i]~dbern(p[i])
    #likelihood
    logit(p[i])<-b[1]+u[subj[i],1] + (b[2]+u[subj[i],2])*age[i]+(b[3]+u[subj[i],3])*disease.duration[i]+
      (b[4]+u[subj[i],4])*edss[i]+(b[5]+u[subj[i],5])*nr.Gd.enhanced.lesions[i]+(b[6]+u[subj[i],6])*nr.relapses.2y.prior.study[i]+
      (b[7]+u[subj[i],7])*months.since.last.relapse[i]+(b[8]+u[subj[i],8])*treatment.naive.prior.visit[i] +
      (b[9]+u[subj[i],9])*gender[i] + (b[10]+u[subj[i],10])*treatment.time.during.cycle.months[i]


  }


  #priors for fixed effects for b[i] (fixef intercept and slopes)
  for (i in 1:(npf+1)){
    b[i]~dnorm(0,0.001)
  }


  for(i in 1:(npf+1)){
    zero[i]<-0
  }

  for (i in 1:npid){
    u[i,1:(npf+1)]~dmnorm(zero, Omega.u) #precision matrix
  }

  for(i in 1:(npf+1)){
    Omega.u[i,i]<-1/(pow(sigma,2))
  }

  for(i in 1:npf) {
    for(j in (i+1):(npf+1)){
      Omega.u[i,j]<-0
    }
  }

  for (i in 2:(npf+1)){
    for(j in (1:(i-1))){
      Omega.u[i,j]<-0
    }
  }

  sigma~dunif(0,2)

}


#give the data
jagsdataSMSC <- list(
  Nobservations=nrow(SMSCdataC),
  outcome=as.numeric(factor(SMSCdataC$relapse.2y.after.study))-1,
  npf=9,
  npid=length(unique(SMSCdataC$patient.id)),
  subj=as.integer(factor(SMSCdataC$patient.id)),
  age=SMSCdataC$age-mean(SMSCdataC$age),
  disease.duration=SMSCdataC$disease.duration-mean(SMSCdataC$disease.duration),
  edss=SMSCdataC$edss-mean(SMSCdataC$edss),
  nr.Gd.enhanced.lesions=SMSCdataC$nr.Gd.enhanced.lesions,
  nr.relapses.2y.prior.study=SMSCdataC$nr.relapses.2y.prior.study,
  months.since.last.relapse=SMSCdataC$months.since.last.relapse-mean(SMSCdataC$months.since.last.relapse),
  treatment.naive.prior.visit=SMSCdataC$treatment.naive.prior.visit ,
  gender=as.factor(SMSCdataC$gender),
  treatment.time.during.cycle.months=SMSCdataC$treatment.time.during.cycle.months-mean(SMSCdataC$treatment.time.during.cycle.months)

)

# run the jugs model
SMSCjagsResults <- jags.parallel(data = jagsdataSMSC ,inits=NULL,parameters.to.save = c('b','sigma'),model.file = jagsmodelSMSC2,
                                 n.chains=2,n.iter = 100000,n.burnin = 1000,n.thin = 10)

traceplot(SMSCjagsResults)
########################## sigma[i]^2 in the diagonials ##############################################


jagsmodelSMSC2 <-function()
{


  for( i in 1:Nobservations)
  {
    outcome[i]~dbern(p[i])
    #likelihood
    logit(p[i])<-b[1]+u[subj[i],1] + (b[2]+u[subj[i],2])*age[i]+(b[3]+u[subj[i],3])*disease.duration[i]+
      (b[4]+u[subj[i],4])*edss[i]+(b[5]+u[subj[i],5])*nr.Gd.enhanced.lesions[i]+(b[6]+u[subj[i],6])*nr.relapses.2y.prior.study[i]+
      (b[7]+u[subj[i],7])*months.since.last.relapse[i]+(b[8]+u[subj[i],8])*treatment.naive.prior.visit[i] +
      (b[9]+u[subj[i],9])*gender[i] + (b[10]+u[subj[i],10])*treatment.time.during.cycle.months[i]


  }


  #priors for fixed effects for b[i] (fixef intercept and slopes)
  for (i in 1:(npf+1)){
    b[i]~dnorm(0,0.001)
  }


  for(i in 1:(npf+1)){
    zero[i]<-0
  }

  for (i in 1:npid){
    u[i,1:(npf+1)]~dmnorm(zero, Omega.u) #precision matrix
  }

  for(i in 1:(npf+1)){
    Omega.u[i,i]<-1/(pow(sigma[i],2))
  }

  for(i in 1:npf) {
    for(j in (i+1):(npf+1)){
      Omega.u[i,j]<-0
    }
  }

  for (i in 2:(npf+1)){
    for(j in (1:(i-1))){
      Omega.u[i,j]<-0
    }
  }
  for(i in 1:(npf+1)){
  sigma[i]~dunif(0,2)
  }
}


#give the data
jagsdataSMSC <- list(
  Nobservations=nrow(SMSCdataC),
  outcome=as.numeric(factor(SMSCdataC$relapse.2y.after.study))-1,
  npf=9,
  npid=length(unique(SMSCdataC$patient.id)),
  subj=as.integer(factor(SMSCdataC$patient.id)),
  age=SMSCdataC$age-mean(SMSCdataC$age),
  disease.duration=SMSCdataC$disease.duration-mean(SMSCdataC$disease.duration),
  edss=SMSCdataC$edss-mean(SMSCdataC$edss),
  nr.Gd.enhanced.lesions=SMSCdataC$nr.Gd.enhanced.lesions,
  nr.relapses.2y.prior.study=SMSCdataC$nr.relapses.2y.prior.study,
  months.since.last.relapse=SMSCdataC$months.since.last.relapse-mean(SMSCdataC$months.since.last.relapse),
  treatment.naive.prior.visit=SMSCdataC$treatment.naive.prior.visit ,
  gender=as.factor(SMSCdataC$gender),
  treatment.time.during.cycle.months=SMSCdataC$treatment.time.during.cycle.months-mean(SMSCdataC$treatment.time.during.cycle.months)

)

# run the jugs model
SMSCjagsResults <- jags.parallel(data = jagsdataSMSC ,inits=NULL,parameters.to.save = c('b','sigma'),model.file = jagsmodelSMSC2,
                                 n.chains=2,n.iter = 100000,n.burnin = 1000,n.thin = 10)

traceplot(SMSCjagsResults)


