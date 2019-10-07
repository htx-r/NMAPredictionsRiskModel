############################################################
#         Master analysis for MS NMA Prediction MODEL
############################################################

##########################################################
############### LIBRARIES #################################

###load the needed libraries
library(devtools)
install_github("htx-r/CleaningData",force=TRUE)
install_github("htx-r/NMAPredictionsRiskModel", force = TRUE)
library(NMAPredictionsRiskModel)
library(CleaningData)
library(R2jags)
library(dplyr)
library(glmpath)
library(readxl)
library(car)
library(glmnet)
library(Hmisc)
library(rms)
library(gridExtra)
library(ggpubr)
library(ggplot2)
library(synthpop)
library(pmsampsize)
library(selectiveInference)
library(plyr)
#######################################################################################
####################################  DATA   ###########################################

###### Give your path of data
mydatapath="C:/Users/kc19o338/Desktop/Real world predictions project/HTx/data/IPD data from 6 Biogen trials"
mydatapath1="C:/Users/kc19o338/Desktop/Real world predictions project/HTx/data cleaning"
######## load data
###cleaning the data from BIOGEN, defined the outcomes in columns: RELAPSE02Year, RELAPSE01Year, names of Treatments and Drugs
cleanBIOGENtrials<-cleanBIOGENtrials.fun(mydatapath)
PlaceboArms<-cleanPLACEBOtrials.fun(mydatapath1)
adsl01<-cleanBIOGENtrials$adsl01
### Select variables that I need- exclude variables with a huge ammount of missing values,
#exclude factors with just one category, exclude factors that are transformations from already existing variables)
#exclude highly correlated variables
###and recode them in numerical values (e.g. Male=1, Female=0)
MSrelapse<-numericalDataRisk.fun(adsl01)  ##final full dataset

#######################################################################################
############################ RISK MODEL ###############################################
######################################################################################
source('EPVandSampleSize.R')
####################### CHECK DIFFERENT RISK MODELS + SHRINKAGE ###############

######## Model 1 - results of Internal risk score
InternalModel<-RiskModelSelection.fun(MSrelapse,"Internal")
######## Model 2 - model Internal with splines to only significant non-linear continues variables
InternalSplinesSignModel<-RiskModelSelection.fun(MSrelapse,"InternalSplinesSign")
######  Model 3 - model Internal with Interactions
InternalInteractionsModel<-RiskModelSelection.fun(MSrelapse,"InternalInteractions")
#########  Model 4 - Fabio's model
FabioModel<-RiskModelSelection.fun(MSrelapse,"Fabio")

###comparison of models to select the best one with respect to discrimination and calibration of each model
source('ComparisonOfModels.R')
ComparisonOfModelsTable

### Chosen model (internal + sign. splines - model 2)
## and graphs for this model


#the final model and the risk data (with extra columns: Risk and logitRisk)
source('FinalRiskModel.R')
## calibration plot of the model
source('Calibration plot.R')
###plots of risk score
source('Plots.R')
RiskDist #distribution of Risk in the whole dataset
RandomizationRisk  #distribution of Risk in each of the studies - randomization of risk
PrognosticRisk #distribution of Risk for those who relapsed and those who did not relapse - Prognostic factor
EffectModRisk #distribution of Risk for those who relapsed and those who did not relapse in each arm and study - Effect modifier


#######################################################################################
####################### NMA PREDICTION MODEL ###############################################
######################################################################################
#add proper columns in the RiskData, like arm, meanRisk, etc.
source('DataForIPDNMR.R')

#run the model & results - it needs some time (around 5 minutes)
IPDNMRJAGSmodel <- jags.parallel(data = jagsdataIPDNMR ,inits=NULL,parameters.to.save = c('be', 'Beta', 'ORref','u','logitp'),model.file = modelIPDNMR,
                                        n.chains=2,n.iter = 100000,n.burnin = 1000,DIC=F,n.thin = 10)
#IPDNMRJAGSmodelFORlogitp <- jags.parallel(data = jagsdataIPDNMR ,inits=NULL,parameters.to.save = c('logitp'),model.file = modelIPDNMR,
#                                 n.chains=2,n.iter = 100000,n.burnin = 1000,DIC=F,n.thin = 10)

print(IPDNMRJAGSmodel,varname=c("be","ORref","u"))
#print(IPDNMRJAGSmodelFORlogitp)
#credible intervals: IPDNMRJAGSmodelFORlogitp$BUGSoutput$summary[,3]
# traceplots
traceplot(IPDNMRJAGSmodel$BUGSoutput,varname=c("be","ORref","u"))

####plot of IPD NMR
source('GraphForPredictedRisk.R')
IPDplot



##remove list
rm(list=ls())
