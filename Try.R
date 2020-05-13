###############         Script to check what is the problem with the model   #####################



#### 1. Let's try to simplify the variance-covariance matrix
#### We will use here as precision matrix the diagonial matrix 1/sigma^2 with dmnorm


jagsmodelSMSC2 <-function()
{


  for( i in 1:Nobservations)
  {
    outcome[i]~dbern(p[i])
    #likelihood
    logit(p[i])<-b[1]+u[subj[i],1] + (b[2]+u[subj[i],2])*age[i]+(b[3]+u[subj[i],3])*disease.duration[i]+
      (b[4]+u[subj[i],4])*edss+(b[5]+u[subj[i],5])*nr.Gd.enhanced.lesions+(b[6]+u[subj[i],6])*nr.relapses.2y.prior.study+
      (b[7]+u[subj[i],7])*months.since.last.relapse+(b[8]+u[subj[i],8])*treatment.naive.prior.visit +
      (b[9]+u[subj[i],9])*gender + (b[10]+u[subj[i],10])*treatment.time.during.cycle.months+h[i]

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
  SMSCjagsResults <- jags.parallel(data = jagsdataSMSC ,inits=NULL,parameters.to.save = c('b'),model.file = jagsmodelSMSC2,
                                   n.chains=2,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 10)
  #### 2. Let's try to simplify the variance-covariance matrix
  #### We will use here as precision matrix the diagonial matrix 1/sigma^2 with dmnorm


  jagsmodelSMSC2 <-function()
  {


    for( i in 1:Nobservations)
    {
      outcome[i]~dbern(p[i])
      #likelihood
      logit(p[i])<-b[1]+u[subj[i],1] + (b[2]+u[subj[i],2])*age[i]+(b[3]+u[subj[i],3])*disease.duration[i]+
        (b[4]+u[subj[i],4])*edss+(b[5]+u[subj[i],5])*nr.Gd.enhanced.lesions+(b[6]+u[subj[i],6])*nr.relapses.2y.prior.study+
        (b[7]+u[subj[i],7])*months.since.last.relapse+(b[8]+u[subj[i],8])*treatment.naive.prior.visit +
        (b[9]+u[subj[i],9])*gender + (b[10]+u[subj[i],10])*treatment.time.during.cycle.months+h[i]

      h[i]~dnorm(0,prec.h)
    }


    #priors for fixed effects for b[i] (fixef intercept and slopes)
    for (i in 1:(npf+1)){
      b[i]~dnorm(0,0.001)
    }


    for(i in 1:npid){

    for (j in 1:(npf+1)){
      u[i,j]~dnorm(0, prec.u[j]) #precision matrix
    }
    }
for(i in 1:(npf+1)){
    prec.u[i]<-1/tau.u[i]
    tau.u[i]<-pow(sigma.u[i], 2)
    sigma.u[i]~dunif(0,10)
}


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
  SMSCjagsResults <- jags.parallel(data = jagsdataSMSC ,inits=NULL,parameters.to.save = c('b'),model.file = jagsmodelSMSC2,
                                   n.chains=2,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 10)

  #### 3. Let's try to simplify the variance-covariance matrix
  #### We will use here as precision matrix the diagonial matrix 1/sigma^2 with dmnorm.vcov where you do not use the inverse matrix


  jagsmodelSMSC3 <-function()
  {


    for( i in 1:Nobservations)
    {
      outcome[i]~dbern(p[i])
      #likelihood
      logit(p[i])<-b[1]+u[subj[i],1] + (b[2]+u[subj[i],2])*age[i]+(b[3]+u[subj[i],3])*disease.duration[i]+
        (b[4]+u[subj[i],4])*edss+(b[5]+u[subj[i],5])*nr.Gd.enhanced.lesions+(b[6]+u[subj[i],6])*nr.relapses.2y.prior.study+
        (b[7]+u[subj[i],7])*months.since.last.relapse+(b[8]+u[subj[i],8])*treatment.naive.prior.visit +
        (b[9]+u[subj[i],9])*gender + (b[10]+u[subj[i],10])*treatment.time.during.cycle.months+h[i]

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
      u[i,1:(npf+1)]~dmnorm.vcov(zero, Omega.u) #precision matrix
    }

    for(i in 1:(npf+1)){
      Omega.u[i,i]<-(pow(sigma,2))
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
  SMSCjagsResults <- jags.parallel(data = jagsdataSMSC ,inits=NULL,parameters.to.save = c('b'),model.file = jagsmodelSMSC2,
                                   n.chains=2,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 10)

## 4. Let's do an empty model with a random intercept only


  jagsmodelSMSC4 <-function()
  {

    for( i in 1:Nobservations)
    {
      outcome[i]~dbern(p[i])
      #likelihood
      logit(p[i])<-b+u[subj[i]]
      h[i]~dnorm(logit(p[i]),prec.h)
    }


    #priors for fixed effects for b[i] (fixef intercept and slopes)

      b~dnorm(0,0.001)

   for(i in 1:npid){

      u[i]~dnorm(0, prec.u) #precision
   }

      prec.u<-pow(sigma.u,-2)
      sigma.u~dunif(0,10)



    # priors for variance
    prec.h<-pow(sigma.h, -2)
    sigma.h~dunif(0,10)

  }


  #give the data
  jagsdataSMSC <- list(
    Nobservations=nrow(SMSCdataC),
    outcome=as.numeric(factor(SMSCdataC$relapse.2y.after.study))-1,
    npid=length(unique(SMSCdataC$patient.id)),
    subj=as.integer(factor(SMSCdataC$patient.id))

  )

  # run the jugs model
  SMSCjagsResults <- jags.parallel(data = jagsdataSMSC ,inits=NULL,parameters.to.save = c('b'),model.file = jagsmodelSMSC4,
                                   n.chains=2,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 10)

  # 5. Let's try a totally empty model with just a fixed intercept


  jagsmodelSMSC5 <-function()
  {

    for( i in 1:Nobservations)
    {
      outcome[i]~dbern(p[i])
      #likelihood
      logit(p[i])<-b
      h[i]~dnorm(logit(p[i]),prec.h)
    }


    #priors for fixed effects for b[i] (fixef intercept and slopes)

    b~dnorm(0,0.001)

    # priors for variance
    prec.h<-pow(sigma.h, -2)
    sigma.h~dunif(0,10)

  }


  #give the data
  jagsdataSMSC <- list(
    Nobservations=nrow(SMSCdataC),
    outcome=as.numeric(factor(SMSCdataC$relapse.2y.after.study))-1


  )

  # run the jugs model
  SMSCjagsResults <- jags.parallel(data = jagsdataSMSC ,inits=NULL,parameters.to.save = c('b'),model.file = jagsmodelSMSC5,
                                   n.chains=2,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 10)

