
modelIPDADMA<-function(){

  #############Part I: model for IPD data
  for (i in 1:Np){
    outcome[i]~dbern(p[i])
    ###formula
    logit(p[i])<-u[studyid[i]] + d[studyid[i], treatma[i]]
  }


  ##Part II: Model for AD data

 for(i in (N.IPD.studies+1):(N.AD.studies +N.IPD.studies)){

    ##likelihood
     outcome.ad[i,1]~dbin(pa[i,1], n[i,1])
      outcome.ad[i,2]~dbin(pa[i,2], n[i,2])

    ### formula for study's reference
    logit(pa[i, 1] )<-u[i]
    logit(pa[i, 2] )<-u[i]+d[i, 2]

 }

  ####PART III: Model for combining all treatment effect estimates
  #Vague priors

  ##independent ui for each study
  for (i in 1:(N.AD.studies +N.IPD.studies)){
    d[i,1]<-0
    u[i]~dnorm(0,0.001)
    d[i,2]<-delta
  }
  delta~dnorm(0,0.001)
 # d~dnorm(0,0.001)



}


