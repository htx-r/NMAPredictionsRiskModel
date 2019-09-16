############################ IPD NMR Bayesian Model ###########################3

modelIPDNMR<-function(){
  ###likelihood
  for (i in 1:Np){
    outcome[i]~dbern(p[i])
    ###formula
    logit(p[i])<-u[studyid[i]] + delta[studyid[i], arm[i]] + beta[studyid[i]]*(Risk[i]-meanRisk[studyid[i]]) + beta.eff[studyid[i],arm[i]]*(Risk[i]-meanRisk[studyid[i]])
  }


  ###fixed effects for u[i]

  for (i in 1:Nstudies){
    u[i]~dnorm(0,0.001)
    beta[i]<-Beta
  }
  Beta~dnorm(0,0.001)

  #####treatment effect - random effects
  for(i in 1:Nstudies){
    delta[i,1] <- 0
    w[i,1] <- 0
    beta.eff[i,1]<-0

    for(k in 2:na[i]){
      delta[i,k]<-md[i, k]
      md[i, k] <- mean[i, k] + sw[i, k]
      w[i, k] <- (delta[i, k] - mean[i, k])
      sw[i, k] <- sum(w[i, 1:(k - 1)])/(k - 1)
      mean[i, k] <- d[treat[i, k]] - d[treat[i, 1]]

      beta.eff[i,k]<-be[treat[i,k]]-be[treat[i,1]]

    }
  }


  ###priors
  d[ref] <- 0 # treatment effect is zero for reference treatment = PLACEBO
  be[ref] <- 0

  for (k in 1:(ref-1)){
    d[k] ~ dnorm(0, 0.01)
    be[k] ~ dnorm(0, 0.01)

  }
  for (k in (ref+1):nt){
    d[k] ~ dnorm(0, 0.01)
    be[k] ~ dnorm(0, 0.01)


  }

  for(j in 1:nt){ORref[j]<- exp(d[j] - d[4])}
}
