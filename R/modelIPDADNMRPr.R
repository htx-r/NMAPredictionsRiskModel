##############################################################################################################################
##############  Model for IPD and AD NMR, with Risk only as prognostic factor ###############################################
##############################################################################################################################


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
