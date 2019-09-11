############################################################
#         Master analysis for MS NMA Prediction MODEL
############################################################

##########################################################
############### LIBRARIES #################################

###load the github libraries
library(devtools)
install_github("htx-r/CleaningData",force=TRUE)
install_github("htx-r/NMAPredictionsRiskModel", force = TRUE)
library(NMAPredictionsRiskModel)
library(CleaningData)
library(R2jags)
library(dplyr)
library(glmpath)

#######################################################################################
####################################  DATA   ###########################################

###### Give your path of data
#mydatapath="C:/Users/kc19o338/Desktop/Real world predictions project/HTx/data/IPD data from 6 Biogen trials"
mydatapath="~/Google Drive/_mydrive/HTx/Biogen data/IPD data from 6 Biogen trials"
######## load data
###cleaning the data from BIOGEN, defined the outcomes in columns: RELAPSE02Year, RELAPSE01Year, names of Treatments and Drugs
cleanBIOGENtrials<-cleanBIOGENtrials.fun(mydatapath)
adsl01<-cleanBIOGENtrials$adsl01
### Select variables that I need- exclude variables with a huge ammount of missing values,
#exclude factors with just one category, exclude factors that are transformations from already existing variables)
#exclude highly correlated variables
###and recode them in numerical values (e.g. Male=1, Female=0)
MSrelapse<-numericalDataRisk.fun(adsl01)  ##final full dataset

#######################################################################################
############################ RISK MODEL ###############################################
######################################################################################

####################### CHECK DIFFERENT RISK MODELS + SHRINKAGE ###############

#########results of Internal risk score - model 1
InternalModel<-RiskModelSelection.fun(MSrelapse,"Internal")
######## model Internal with splines to all continues variables
#model 2
InternalSplinesModel<-RiskModelSelection.fun(MSrelapse,"InternalSplines")
######## model Internal with splines to only significant non-linear
# continues variables
#model 3
InternalSplinesSignModel<-RiskModelSelection.fun(MSrelapse,"InternalSplinesSign")
###### model Internal with Interactions - model 4
InternalInteractionsModel<-RiskModelSelection.fun(MSrelapse,"InternalInteractions")
### Cross Internal model - model 5
CrossInternalModel<-RiskModelSelection.fun(MSrelapse,"CrossInternal")
######### Fabio's model - model 6
FabioModel<-RiskModelSelection.fun(MSrelapse,"Fabio")##!!!!!!!!!!!!! the putput of the function should be the validate function discrimination and calibration for the 4 shrinkage approaches


#!!!!the RiskModelSelection.fun shoud have output list(discrimination=discrimiation, claibration=claibration)
##!!! discrimination is from the validate (4x1) same for claibration
#!!!!!!here you need a script that creates the table with all the comparison of the models

#creation of table 2
comparisonofmodelstable<-cbind.data.frame(c(InternalSplinesSignModel$discrimination......), c(InternalSplinesSignModel$calibration))
# 12x2 matrix

### Chosen model (internal + sign. splines - model 3)
## and graphs for this model
Risk<-FinalRiskModel.fun(MSrelapse)###!! this should be a script

#data including the risk and logitof risk for each patient (2 extra columns)
RiskData<-Risk[[2]]#!!!include in scirpt
## the lrm final model
RiskModel<-Risk[[1]]#!!!include in scirpt
###plots of risk score
source('Plots.R')
RiskDist
RandomizationRisk
PrognosticRisk
EffectModRISK
#######################################################################################
####################### NMA PREDICTION MODEL ###############################################
######################################################################################
#add proper columns in the RiskData, like arm, meanRisk, etc.
RiskData<-DataForIPDNMR.fun()##!!!! script

###data for jagsmodel with metarigression on logit of Risk
##!!!! script with the precvious
jagsdataIPDNMR <- list(
  Nstudies=3,
  Np=sum(as.numeric(table(as.numeric(as.factor(RiskData$STUDYID))))),
  studyid=as.numeric(RiskData$STUDYID),
  outcome=RiskData$RELAPSE2year,
  treat= rbind(c(1,4,NA),c(1,2,4),c(3,4,NA)),
  na=c(2,3,2),
  arm=RiskData$arm,
  Risk=RiskData$logitRisk,
  meanRisk=c(-0.5360,-0.6501,-0.5110), ##here is the mean of logit of risk
  nt=4,
  ref=4
)

#run the model - it needs some time (around 5 minutes)
IPDNMRJAGSmodel <- jags.parallel(data = jagsdataIPDNMR ,inits=NULL,parameters.to.save = c('d','be', 'Beta', 'ORref','u'),model.file = modelIPDNMR,
                                        n.chains=2,n.iter = 100000,n.burnin = 1000,DIC=F,n.thin = 10)
print(IPDNMRJAGSmodel)
# traceplots
traceplot(IPDNMRJAGSmodel$BUGSoutput)
###plot of IPD NMR
source('graph ipd.R')
IPDplot

##remove list
rm(list=ls())
