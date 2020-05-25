
#########################################  Independent uis  #######################################################

#  The model

jagsmodelSMSC2 <-function(){
  for( i in 1:Nobservations)
    {
    outcome[i]~dbern(p[i])
    #likelihood
    logit(p[i])<-b[1]+u[subj[i],1] + (b[2]+u[subj[i],2])*age[i]+(b[3]+u[subj[i],3])*disease.duration[i]+
      (b[4]+u[subj[i],4])*edss[i]+(b[5]+u[subj[i],5])*nr.Gd.enhanced.lesions[i]+(b[6]+u[subj[i],6])*nr.relapses.2y.prior.study[i]+(b[7]+u[subj[i],7])*months.since.last.relapse[i]+(b[8]+u[subj[i],8])*treatment.naive.prior.visit[i] + (b[9]+u[subj[i],9])*gender[i] + (b[10]+u[subj[i],10])*treatment.time.during.cycle.months[i]+h[i]
    h[i]~dnorm(0,prec.h)
    }
  #priors for fixed effects b[i] (fixed efffects intercept and slopes)
  for (i in 1:(npf+1))
    {
    b[i]~dnorm(0,0.001)
    }
  for(i in 1:npid)
    {
    for (j in 1:(npf+1))
      {
      u[i,j]~dnorm(0, prec.u[j])
    }
  }
    for(i in 1:(npf+1))
      {
      prec.u[i]<-1/sigma2.u[i]
      sigma2.u[i]<-pow(sigma.u[i], 2)
      sigma.u[i]~dunif(0,2)
      }
    # priors for variance
    prec.h<-1/sigma2.h
    sigma2.h<-pow(sigma.h, 2)
    sigma.h~dunif(0,2)
    }

  # The data
  jagsdataSMSC <- list(
    Nobservations=nrow(SMSCdataC),
    outcome=as.numeric(factor(SMSCdataC$relapse.2y.after.study))-1,
    npf=9,
    npid=length(unique(SMSCdataC$patient.id)),
    subj=as.integer(factor(SMSCdataC$patient.id)),
    age=SMSCdataC$age,
    disease.duration=SMSCdataC$disease.duration,
    edss=SMSCdataC$edss,
    nr.Gd.enhanced.lesions=SMSCdataC$nr.Gd.enhanced.lesions,
    nr.relapses.2y.prior.study=SMSCdataC$nr.relapses.2y.prior.study,
    months.since.last.relapse=SMSCdataC$months.since.last.relapse,
    treatment.naive.prior.visit=SMSCdataC$treatment.naive.prior.visit ,
    gender=as.factor(SMSCdataC$gender),
    treatment.during.cycle =as.factor(SMSCdataC$treatment.during.cycle),
    treatment.time.during.cycle.months=SMSCdataC$treatment.time.during.cycle.months

  )

  # run the jugs model
  SMSCjagsResults <- jags.parallel(data = jagsdataSMSC ,inits=NULL,parameters.to.save = c('b', 'sigma2.h', 'sigma2.u'),model.file = jagsmodelSMSC2,
                                   n.chains=2,n.iter = 100000,n.burnin = 1000,DIC=F,n.thin = 10)
  traceplot(SMSCjagsResults)

  #################################### Unstructured ###############################################################
  jagsmodelSMSC <-function()
  {


    for( i in 1:Nobservations)
    {
      outcome[i]~dbern(p[i])
      #likelihood
      logit(p[i])<-b[1]+u[subj[i],1] + (b[2]+u[subj[i],2])*age[i]+(b[3]+u[subj[i],3])*disease.duration[i]+
        (b[4]+u[subj[i],4])*edss[i]+(b[5]+u[subj[i],5])*nr.Gd.enhanced.lesions[i]+(b[6]+u[subj[i],6])*nr.relapses.2y.prior.study[i]+
        (b[7]+u[subj[i],7])*months.since.last.relapse[i]+(b[8]+u[subj[i],8])*treatment.naive.prior.visit[i] +
        (b[9]+u[subj[i],9])*gender[i] + (b[10]+u[subj[i],10])*treatment.time.during.cycle.months[i]+h[i]

      h[i]~dnorm(0,prec.h)
    }


    #priors for fixed effects for b[i] (fixef intercept and slopes)
    for (i in 1:(npf+1)){
      b[i]~dnorm(0,0.001)
    }


    for(i in 1:(npf+1)){
      zero[i]<-0
    }

    for (i in 1:npid){
      u[i,1:(npf+1)]~dmnorm(zero, Omega.u)
    }
    ##fill in the R.u

    for(i in 1:(npf+1)) {

      R.u[i,i]<-sigma[i]*sigma[i]
    }

    for(i in 1:npf) {
      for(j in (i+1):(npf+1)){
        R.u[i,j]<-rho.u*sigma[i]*sigma[j]
      }
    }

    for (i in 2:npf){
      for(j in (1:(i-1))){
        R.u[i,j]<-rho.u*sigma[i]*sigma[j]
      }
    }



    #priors for varying intercept and slopes
    for(i in (1:(npf+1))){
      sigma[i]~dunif(0,10)
    }

    #prior for correlation
    rho.u~dunif(-1,1)

    #prior on precision Omega.u
    Omega.u~dwish(R.u,npf+1)

    # priors for variance
    prec.h<-1/tau.h
    tau.h<-pow(sigma.h, 2)
    sigma.h~dunif(0,10)

  }

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

  SMSCjagsResults <- jags.parallel(data = jagsdataSMSC ,inits=NULL,parameters.to.save = c('b', 'sigma.h', 'sigma', 'rho'),model.file = jagsmodelSMSC,
                                   n.chains=2,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 10)

  ########################## Only sigma^2 in the diagonials ##############################################


  jagsmodelSMSC2 <-function()
  {


    for( i in 1:Nobservations)
    {
      outcome[i]~dbern(p[i])
      #likelihood
      logit(p[i])<-b[1]+u[subj[i],1] + (b[2]+u[subj[i],2])*age[i]+(b[3]+u[subj[i],3])*disease.duration[i]+
        (b[4]+u[subj[i],4])*edss[i]+(b[5]+u[subj[i],5])*nr.Gd.enhanced.lesions[i]+(b[6]+u[subj[i],6])*nr.relapses.2y.prior.study[i]+
        (b[7]+u[subj[i],7])*months.since.last.relapse[i]+(b[8]+u[subj[i],8])*treatment.naive.prior.visit[i] +
        (b[9]+u[subj[i],9])*gender[i] + (b[10]+u[subj[i],10])*treatment.time.during.cycle.months[i]+h[i]

      h[i]~dnorm(0,prec.h)
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

    sigma~dunif(0,10)


    # priors for variance
    prec.h<-1/tau.h
    tau.h<-pow(sigma.h, 2)
    sigma.h~dunif(0,10)

  }


  #give the data
  jagsdataSMSC <- list(
    Nobservations=nrow(SMSCdataC),
    outcome=as.numeric(factor(SMSCdataC$relapse.2y.after.study))-1,
    npf=9,
    npid=length(unique(SMSCdataC$patient.id)),
    subj=as.integer(factor(SMSCdataC$patient.id)),
    age=SMSCdataC$age,
    disease.duration=SMSCdataC$disease.duration,
    edss=SMSCdataC$edss,
    nr.Gd.enhanced.lesions=SMSCdataC$nr.Gd.enhanced.lesions,
    nr.relapses.2y.prior.study=SMSCdataC$nr.relapses.2y.prior.study,
    months.since.last.relapse=SMSCdataC$months.since.last.relapse,
    treatment.naive.prior.visit=SMSCdataC$treatment.naive.prior.visit ,
    gender=as.factor(SMSCdataC$gender),
    treatment.during.cycle =as.factor(SMSCdataC$treatment.during.cycle),
    treatment.time.during.cycle.months=SMSCdataC$treatment.time.during.cycle.months

  )

  # run the jugs model
  SMSCjagsResults <- jags.parallel(data = jagsdataSMSC ,inits=NULL,parameters.to.save = c('b','sigma', 'sigma.h'),model.file = jagsmodelSMSC2,
                                   n.chains=2,n.iter = 100000,n.burnin = 1000,n.thin = 10)


