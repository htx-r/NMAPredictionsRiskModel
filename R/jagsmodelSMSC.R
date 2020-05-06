jagsmodelSMSC <-function()
{


for( i in 1:Nobservations)
{
  outcome[i]~dbern(p[i])
  #likelihood
  logit(p[i])<-b[1]+u[subj[i],1] + (b[2]+u[subj[i],2])*age[i]+(b[3]+u[subj[i],3])*disease.duration[i]+
    (b[4]+u[subj[i],4])*edss+(b[5]+u[subj[i],5])*nr.Gd.enhanced.lesions+(b[6]+u[subj[i],6])*nr.relapses.2y.prior.study+
    (b[7]+u[subj[i],7])*months.since.last.relapse+(b[8]+u[subj[i],8])*treatment.naive.prior.visit +
    (b[9]+u[subj[i],9])*gender+
    (b[10]+u[subj[i],10])*treatment.during.cycle+ (b[11]+u[subj[i],11])*treatment.time.during.cycle.months+h[i]

    h[i]~dnorm(0,tau.h)
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

    for(i in 1:11) {

        R.u[i,i]<-sigma[i]*sigma[i]
      }

   for(i in 1:10) {
      for(j in (i+1):11){
        R.u[i,j]<-rho.u*sigma[i]*sigma[j]
      }
   }

    for (i in 2:11){
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
  Omega.u~dwish(R.u,(npf+1))

  # priors for variance
  tau.h<-pow(sigma.h, -2)
  sigma.h~dunif(0,10)

}




