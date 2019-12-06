modelIPDADNMRPr<-function(){

  #############Part I: model for IPD data
  for (i in 1:Np){
    outcome[i]~dbern(p[i])
    ###formula
    logit(p[i])<-u[studyid[i]] + d[studyid[i], arm[i]] + g[studyid[i]]*(Risk[i]-meanRisk[studyid[i]])
  }

  #####treatment effects - fixed across studies & correction for multi-arm studies
  for(i in 1:N.IPD.studies){
    d[i,1] <- 0
    w[i,1] <- 0

    for(k in 2:na[i]){
      d[i,k]<-md[i, k]
      md[i, k] <- mean[i, k] + sw[i, k]
      w[i, k] <- (d[i, k] - mean[i, k])
      sw[i, k] <- sum(w[i, 1:(k - 1)])/(k - 1)
      mean[i, k] <- delta[treat[i, k]] - delta[treat[i, 1]]

    }

  }

  ##fixed across studies for g
  for (i in 1:N.IPD.studies){
    g[i]<-gamma
  }
  #vague prior for gamma
  gamma~dnorm(0,0.001)

  ##Part II: Model for AD data


  for(i in (N.IPD.studies+1):(N.AD.studies +N.IPD.studies)){
    d[i,treat[i,1] ]<- 0
    w[i,1] <- 0
    ### formula for study's reference
    logit(pa[i, treat[i,1] ])<-u[i]

    for (k in 1:na[i]){
      ##likelihood
      outcome.ad[i,treat[i,k]]~dbin(pa[i,treat[i,k]], n[i,treat[i,k]])
    }
    for(k in 2:na[i]){
      ### formula - fixed across studies
      logit(pa[i, treat[i,k] ])<-u[i] + d[i, treat[i,k]]

      d[i, treat[i,k]]<-md.ad [i, treat[i,k]]
      md.ad[i, treat[i,k]]<- mean[i, k] + sw[i, k]
      w[i, k] <- (d[i, treat[i, k]] - mean[i, k])
      sw[i, k] <- sum(w[i, 1:(k - 1)])/(k - 1)
      mean[i, k] <- delta[treat[i, k]] - delta[treat[i, 1]]


    }
  }

  ####PART III: Model for combining all treatment effect estimates
  #Vague priors

  ##independent ui for each study
  for (i in 1:(N.IPD.studies+N.AD.studies)){
    u[i]~dnorm(0,0.001)
  }

  delta[ref]<-0

  for(k in 1:(ref-1)){
    delta[k]~dnorm(0,0.001)

  }
  for(k in (ref+1):nt){
    delta[k]~dnorm(0,0.001)

  }



  ###Calculation of ORs
  for(j in 1:nt){ORref[j]<- exp(delta[j] - delta[ref])}

}


jagsdataIPDADnetmeta <- list(
  N.IPD.studies=3,
  Np=sum(as.numeric(table(as.numeric(as.factor(RiskData$STUDYID))))),
  studyid=as.numeric(as.factor(RiskData$STUDYID)),
  outcome=as.numeric(RiskData$RELAPSE2year)-1,
  treat= rbind(c(1,4,NA),c(1,2,4),c(3,4,NA), c(2,4,NA),c(2,4,NA)),
  na=c(2,3,2,2,2),
  arm=RiskData$arm,
  ref=4,
  nt=4,
  N.AD.studies=2,
  outcome.ad=rbind(c(NA,NA,NA,NA),c(NA,NA,NA,NA),c(NA,NA,NA,NA),c(NA,11,NA,19),c(NA,89,NA,97)),
  n=rbind(c(NA,NA,NA,NA),c(NA,NA,NA,NA),c(NA,NA,NA,NA),c(NA,25,NA,25),c(NA,125,NA,126)),
  Risk=RiskData$logitRiskFabio,
  meanRisk=c(tapply(RiskData$logitRiskFabio, RiskData$STUDYID, summary)$`1`[4],tapply(RiskData$logitRiskFabio, RiskData$STUDYID, summary)$`2`[4],tapply(RiskData$logitRiskFabio, RiskData$STUDYID, summary)$`3`[4]) ##here is the mean of logit of risk


)

####RUN the model
IPDADnetmetaJAGSmodel <- jags.parallel(data = jagsdataIPDADnetmeta ,inits=NULL,parameters.to.save = c('delta','u','ORref','gamma'),model.file = modelIPDADNMRPr,
                                       n.chains=2,n.iter = 10000,n.burnin = 100,DIC=F,n.thin = 1)


