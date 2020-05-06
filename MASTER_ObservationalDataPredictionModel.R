############################################################
#         Master analysis for MS NMA Prediction MODEL
############################################################



##########################################################
############### LIBRARIES #################################
### Load needed variables
library(readxl)
library(dplyr)
library(rms)
#library(gee)
library(lme4)
library(R2jags)

# give/change your datapath (the datapath that includes the data)
mydatapath="C:/Users/kc19o338/Desktop/SMSC Basel"

# function that keeps only the needed variables (selected via pre-existing prognostic models on the literature, for RRMS patients)
# and makses the proper transformations for continues and categorical variables
# (in case you want to see all summary statistics for the SMSC data before and after the transformations you can run the SMSC_Summary.R script)
source("R/FinalDataSMSC.fun.R")
SMSCdata<-FinalDataSMSC.fun(mydatapath)
# the dataset with complete cases only, no missing values at all
SMSCdataC<-na.omit(SMSCdata)
  #create a variable that indicates the cycle of the patient (i.e. 1st cycle, 2nd cycle, 3rd cycle, etc.) for the complete dataset
  SMSCdataC$cycle <- ave(SMSCdataC$unique.visit.id, SMSCdataC$patient.id,  FUN = seq_along)
  #create a column with patient ids as numeric for the complete dataset
  SMSCdataC$NumericID <- as.numeric(factor(SMSCdataC$patient.id,
                                        levels=unique(SMSCdataC$patient.id)))


############################ Frequentist framework ######################################################################



# just a test using generalized linear mixed effects model in a frequentist setting
source("R/Frequentist_glmm.fun.R")
glmmModel<-Frequentist_glmm.fun(SMSCdataC)
    #the results of the mixed-effects model
    summary(glmmModel$glmer_out) #message that the model failed to converge
    # compare them with a model with fixed effects only (glm)
    summary(glmmModel$glm_out)


############################ Bayesian framework ######################################################################

##### Read the data needed for the jags model
    jagsdataSMSC <- list(
      Nobservations=nrow(SMSCdataC),
      outcome=SMSCdataC$outcome,
      npf=10,
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
      treatment.during.cycle =as.factor(SMSCdataC$treatment.during.cycle),
      treatment.time.during.cycle.months=SMSCdataC$treatment.time.during.cycle.months

    )
source("R/jagsmodelSMSC.R")
# run the jugs model
    SMSCjagsResults <- jags.parallel(data = jagsdataSMSC ,inits=NULL,parameters.to.save = c('b'),model.file = jagsmodelSMSC,
                                     n.chains=2,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 10)



