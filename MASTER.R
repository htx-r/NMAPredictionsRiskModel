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
library(vcd)
library(pROC)
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
###drop SENTTINEL STUDY because of combination of treatments and
### drop ADVANCE study because does not provide information for Relapse in 2 years (only for 1 year)
adsl<-adsl01[adsl01$STUDYID!="SENTINEL" & adsl01$STUDYID!="ADVANCE" ,]
### Select variables that I need- exclude variables with a huge ammount of missing values (more than 50%),
#exclude factors with just one category, exclude factors that are transformations from already existing variables)
#exclude highly correlated variables
###and recode them in numerical values (e.g. Male=1, Female=0)
## transformations of continuous variables to approximate normal distribution
MSrelapse<-numericalDataRisk.fun(adsl)  ##final full dataset

#######################################################################################
############################ STAGE 1 - RISK MODEL ###############################################
######################################################################################
source('EPVandSampleSize.R')
####################### CHECK DIFFERENT RISK MODELS + SHRINKAGE ###############

######## Model 1 - results of LASSO model
LASSOModel<-RiskModels.fun(MSrelapse,"LASSOModel")
#########  Model 2 - Results of Pellegrini's model
FabioModel<-RiskModels.fun(MSrelapse,"FabioModel")

###needs more than 3 hours to run! You can skip it without any problem further
###Performance table of models : discrimination and calibration
data1<-na.omit(MSrelapse)
todrop<-c("STUDYID","USUBJID","TRT01A")
data1<-data1[ , !(names(data1) %in% todrop)]
Internal_validation<-BootstrapValidation.fun(data=data1, samples = 500, alpha = 1, modelElasticNet = LASSOModel$lassomodel, modelSpecific = FabioModel$fabiomodel)
Discrimination_Calibration<-as.data.frame(Internal_validation[[7]])
Discrimination_Calibration### bootstrap optimism corrected discriminatio and calibration of the models

### Calibration Plots
#For LASSO model
Calibrationplots.fun(model="LASSOModel")
# For Pellegrini's model
Calibrationplots.fun(model="FabioModel")

### Create a dataset that includes Risk's and logit Risk's predictions for each individual and for both models
##ALso make treatment and studies numerical
source("RiskData.R")

###plots of risk score
## boxplot of both models
source('Plots.R')

#######################################################################################
####################### STAGE 2 - NMA PREDICTION MODEL ###############################################
######################################################################################
#add proper columns in the RiskData, like arm, meanRisk, etc.
source('DataForIPDNMR.R')

#run the model & results - it needs some time (around 5 minutes)
IPDNMRJAGSmodelLASSO <- jags.parallel(data = jagsdataIPDNMRLASSO,inits=NULL,parameters.to.save = c('gamma.w', 'logitpplacebo','gamma', 'ORref','delta','u','logitp'),model.file = modelIPDNMR,
                                        n.chains=2,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 10)

IPDNMRJAGSmodelFabio <- jags.parallel(data = jagsdataIPDNMRFabio,inits=NULL,parameters.to.save = c('gamma.w', 'logitpplacebo','gamma', 'ORref','delta','u','logitp'),model.file = modelIPDNMR,
                                      n.chains=2,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 10)

# Results using LASSO model
print(IPDNMRJAGSmodelLASSO,varname=c("gamma.w","ORref","u","delta"))
# Results using Pellegrini's model
print(IPDNMRJAGSmodelFabio,varname=c("gamma.w","ORref","u","delta"))

#
#credible intervals: IPDNMRJAGSmodelFORlogitp$BUGSoutput$summary[,3]

# traceplots

traceplot(IPDNMRJAGSmodelLASSO$BUGSoutput,varname=c("be","ORref","u"))
traceplot(IPDNMRJAGSmodelFabio$BUGSoutput,varname=c("be","ORref","u"))


####plot of IPD NMR with both models
source('GraphForPredictedRisk.R')

##remove list
rm(list=ls())
