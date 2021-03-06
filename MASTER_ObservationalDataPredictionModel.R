############################################################
#         Master analysis for MS NMA Prediction MODEL
#     with SMSC observational study - repeated measures
#               Generalized mixed effects model
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


############################ Frequentist framework ######################################################################



# just a test using generalized linear mixed effects model in a frequentist setting
source("Frequentist_glmm.R")

    #the results of the mixed-effects model
    summary(glmer_out) #message that the model failed to converge
    # compare them with a model with fixed effects only (glm)
    summary(glm_out)


############################ Bayesian framework ######################################################################

##### Read the data needed for the jags model
source("DataJagsModelSMSC.R")
#source the jags model
source("R/jagsmodelSMSC.R")
# run the jugs model
    SMSCjagsResults <- jags.parallel(data = jagsdataSMSC ,inits=NULL,parameters.to.save = c('b','sigma', 'rho' ),model.file = jagsmodelSMSC,
                                     n.chains=2,n.iter = 100000,n.burnin = 1000,n.thin = 10)
#results
    print(SMSCjagsResults)
    traceplot(SMSCjagsResults)


