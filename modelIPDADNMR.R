modelIPDADNMR<-function(){

  #############Part I: model for IPD data
  for (i in 1:Np){
    outcome[i]~dbern(p[i])
    ###formula
    logit(p[i])<-u[studyid[i]] + d[studyid[i], arm[i]]+ g[studyid[i]]*(Risk[i]-meanRisk[studyid[i]])
    + g.w[studyid[i],arm[i]]*(Risk[i]-meanRisk[studyid[i]])
    + g.b[studyid[i],arm[i]]*meanRisk[studyid[i]]
  }

  #####treatment effects - fixed across studies & correction for multi-arm studies
  for(i in 1:N.IPD.studies){
    d[i,1] <- 0
    w[i,1] <- 0
    g.w[i,1]<-0
    g.b[i,1]<-0
    for(k in 2:na[i]){
      d[i,k]<-md[i, k]
      md[i, k] <- mean[i, k] + sw[i, k]
      w[i, k] <- (d[i, k] - mean[i, k])
      sw[i, k] <- sum(w[i, 1:(k - 1)])/(k - 1)
      mean[i, k] <- delta[treat[i, k]] - delta[treat[i, 1]]

      g.w[i,k]<-gamma.w[treat[i,k]]-gamma.w[treat[i,1]]
      g.b[i,k]<-gamma.b[treat[i,k]]-gamma.b[treat[i,1]]
    }

  }

  ###fixed across studies for beta
  for (i in 1:N.IPD.studies){
    g[i]<-gamma
  }
  ###Vague prior for beta
  gamma~dnorm(0,0.001)

  ##Part II: Model for AD data


  for(i in (N.IPD.studies+1):(N.AD.studies +N.IPD.studies)){
    d[i,treat[i,1] ]<- 0
    g.b[i,treat[i,1]]<-0
    w[i,1] <- 0
    ### formula for study's reference
    logit(pa[i, treat[i,1] ])<-u[i]

    for (k in 1:na[i]){
      ##likelihood
      outcome.ad[i,treat[i,k]]~dbin(pa[i,treat[i,k]], n[i,treat[i,k]])
    }
    for(k in 2:na[i]){
      ### formula - fixed across studies
      logit(pa[i, treat[i,k] ])<-u[i] + d[i, treat[i,k]] + g.b[i,treat[i,k]]*meanRisk[i]

      d[i, treat[i,k]]<-md.ad [i, treat[i,k]]
      md.ad[i, treat[i,k]]<- mean[i, k] + sw[i, k]
      w[i, k] <- (d[i, treat[i, k]] - mean[i, k])
      sw[i, k] <- sum(w[i, 1:(k - 1)])/(k - 1)
      mean[i, k] <- delta[treat[i, k]] - delta[treat[i, 1]]

      g.b[i,treat[i,k]]<-gamma.b[treat[i,k]]-gamma.b[treat[i,1]]
    }
  }

  ####PART III: Model for combining all treatment effect estimates
 #Vague priors

  ##independent ui for each study
  for (i in 1:(N.IPD.studies+N.AD.studies)){
    u[i]~dnorm(0,0.001)
  }

  delta[ref]<-0
  gamma.w[ref]<-0
  gamma.b[ref]<-0
  for(k in 1:(ref-1)){
    delta[k]~dnorm(0,0.001)
    gamma.w[k] ~ dnorm(0, 0.001)
    gamma.b[k] ~ dnorm(0, 0.001)
  }
  for(k in (ref+1):nt){
    delta[k]~dnorm(0,0.001)
    gamma.w[k] ~ dnorm(0, 0.001)
    gamma.b[k] ~ dnorm(0, 0.001)
  }



  ###Calculation of ORs
  for(j in 1:nt){ORref[j]<- exp(delta[j] - delta[ref])}

}



#### Aggregated data
MS <- read_excel("C:/Users/kc19o338/articles/TracemereADMS/MS replication/MS.xlsx")
MS <- as.data.frame(MS)

ADdata<-MS[which(MS$study=="Bornstein" | MS$study=="Johnson"),]


######DATA
## arm column
RiskData$arm<-NA
RiskData$arm[RiskData$STUDYID==1 & RiskData$TRT01A==1]<-1
RiskData$arm[RiskData$STUDYID==1 & RiskData$TRT01A==4]<-2
RiskData$arm[RiskData$STUDYID==2 & RiskData$TRT01A==1]<-1
RiskData$arm[RiskData$STUDYID==2 & RiskData$TRT01A==2]<-2
RiskData$arm[RiskData$STUDYID==2 & RiskData$TRT01A==4]<-3
RiskData$arm[RiskData$STUDYID==3 & RiskData$TRT01A==3]<-1
RiskData$arm[RiskData$STUDYID==3 & RiskData$TRT01A==4]<-2

ADdata$treat[which(ADdata$treat=="Placebo")]<-1
ADdata$treat[which(ADdata$treat=="Glatiramer acetate")]<-2
ADdata$id<-c(4,4,5,5)

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
  meanRisk=c(tapply(RiskData$logitRiskFabio, RiskData$STUDYID, summary)$`1`[4],tapply(RiskData$logitRiskFabio, RiskData$STUDYID, summary)$`2`[4],tapply(RiskData$logitRiskFabio, RiskData$STUDYID, summary)$`3`[4], -0.5, -0.6)##here is the mean of logit of risk


)

####RUN the model
IPDADnetmetaJAGSmodel <- jags.parallel(data = jagsdataIPDADnetmeta ,inits=NULL,parameters.to.save = c('delta','u','ORref','gamma','gamma.b','gamma.w'),model.file = modelIPDADNMR,
                                       n.chains=2,n.iter = 10000,n.burnin = 100,DIC=F,n.thin = 1)
