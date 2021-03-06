############################################################
#         Master analysis for MS NMA Prediction MODEL
############################################################



##########################################################
############### LIBRARIES #################################

###load the needed libraries
library(devtools)
install_github("htx-r/CleaningData",force=TRUE)
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
library(metafor)
library(meta)
library(netmeta)
library(DiagrammeR)

#######################################################################################
####################################  DATA   ###########################################

###### Give your path of data
mydatapath="C:/Users/kc19o338/Desktop/RealWorldPredictionModel/HTx/data/IPD data from 6 Biogen trials"
mydatapath1="C:/Users/kc19o338/Desktop/RealWorldPredictionModel/HTx/Placebo Arms"
######## load data
###cleaning the data from BIOGEN, defined the outcomes in columns: RELAPSE02Year, RELAPSE01Year, names of Treatments and Drugs
cleanBIOGENtrials<-cleanBIOGENtrials.fun(mydatapath)
PlaceboArms<-cleanPLACEBOtrials.fun(mydatapath1)
RCTs0<-cleanBIOGENtrials$adsl01
###drop SENTTINEL STUDY because of combination of treatments and
### drop ADVANCE study because does not provide information for Relapse in 2 years (only for 1 year)
RCTs<-RCTs0[RCTs0$STUDYID!="SENTINEL" & RCTs0$STUDYID!="ADVANCE" ,]
### Select variables that I need- exclude variables with a huge ammount of missing values (more than 50%),
#exclude factors with just one category, exclude factors that are transformations from already existing variables)
#exclude highly correlated variables
###and recode them in numerical values (e.g. Male=1, Female=0)
## transformations of continuous variables to approximate normal distribution
source("R/numericalDataRisk.fun.R")
MSrelapse<-numericalDataRisk.fun(RCTs)  ##final full dataset


#######################################################################################
############################ STAGE 1 - RISK MODEL ###############################################
######################################################################################

#Step 1. Obtain the EPV and the required sample size for the development of the models

source('EPVandSampleSize.R')


#Step 2. Build two different risk models with shrinkage

######## Model 1 - results of LASSO model with LASSO shrinkage

source("R/RiskModels.fun.R")
LASSOModel<-RiskModels.fun(MSrelapse,"LASSOModel")
#########  Model 2 - Results of Pre-specified model using Penalized Maximum Likelihood Estimation for shrinkage
PreSpecifiedModel<-RiskModels.fun(MSrelapse,"PreSpecifiedModel")

#Step 3. Bootstrap calibration and discrimination by using the same bootstrap sample for two models each time
###needs more than 3 hours to run! You can skip it without any problem further
###Performance table of models : discrimination and calibration
#data1<-na.omit(MSrelapse)
#todrop<-c("STUDYID","USUBJID","TRT01A")
#data1<-data1[ , !(names(data1) %in% todrop)]
#Internal_validation<-BootstrapValidation.fun(data=data1, samples = 500, alpha = 1, modelElasticNet = LASSOModel$lassomodel, modelSpecific = PreSpecifiedModel$PreSpecifiedmodel)
#Discrimination_Calibration<-as.data.frame(Internal_validation[[7]])
#Discrimination_Calibration### bootstrap optimism corrected discriminatio and calibration of the models

#Step 4. Create a dataset that includes Risk's and logit Risk's predictions for each individual and for both models
##Also make treatment and studies numerical values
source("RiskData.R")

#Step 5. Source the script for the plots of risk score
source('Plots.R')

#######################################################################################
####################### STAGE 2 - NMR PREDICTION MODEL ###############################################
######################################################################################

#Step 1.  Add proper columns in the RiskData, like arm, meanRisk, and make data for jags model etc.
source('DataForIPDNMR.R')

source("R/modelIPDNMR.R") #source the model will be used for rjags
#Step 2. Run the model & results & check of traceplots
IPDNMRJAGSresultsLASSO <- jags.parallel(data = jagsdataIPDNMRLASSO,inits=NULL,parameters.to.save = c('gamma.w', 'logitpplacebo','gamma', 'ORref','delta','u','logitp'),model.file = modelIPDNMR,
                                      n.chains=2,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 10)

IPDNMRJAGSresultsPreSpecified <- jags.parallel(data = jagsdataIPDNMRPreSpecified,inits=NULL,parameters.to.save = c('gamma.w', 'logitpplacebo','gamma', 'ORref','delta','u','logitp'),model.file = modelIPDNMR,
                                             n.chains=2,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 10)

# Results using LASSO model
print(IPDNMRJAGSresultsLASSO,varname=c("gamma.w","ORref","u","delta"))
# Results using Pellegrini's model
print(IPDNMRJAGSresultsPreSpecified,varname=c("gamma.w","ORref","u","delta"))

# traceplots

traceplot(IPDNMRJAGSresultsLASSO$BUGSoutput,varname=c("ORref","u", "gamma.w", "gamma"))

traceplot(IPDNMRJAGSresultsPreSpecified$BUGSoutput,varname=c("ORref","u", "gamma.w", "gamma"))

#Step 3. Plot of IPD NMR with both models
source('GraphForPredictedRisk.R')


###################################################################################################################################
####################### STAGE 2 - NMR PREDICTION MODEL COMBINING BOTH AD AND IPD ###############################################
##################################################################################################################################



################################################## A. NMA model  #################################################################

## We check the NMA model for IPD and AD

# Here we load the AD, make proper arms in IPD data and make the jagsdata
source('DataForIPDADNMR.R')

source("R/modelIPDADNMA.R")#source the model will be used for rjags
####RUN the model
IPDADNMAJAGSresults <- jags.parallel(data = jagsdataIPDADNMA ,inits=NULL,parameters.to.save = c('delta','u','ORref'),model.file =modelIPDADNMA,
                                       n.chains=2,n.iter = 10000,n.burnin = 100,DIC=F,n.thin = 1)
#traceplots
traceplot(IPDADNMAJAGSresults)
#results
IPDADNMAJAGSresults

#check the results with netmeta
source('CheckNMA.R')
summary(net1, digits = 2)
forest(net1,ref=4,fontsize=10)

####################################### Meta-regression to see what is the problem with the uis ###############################

source("R/modelIPDADMA.R") #source the model will be used for rjags
IPDADMAJAGSresults <- jags.parallel(data = jagsdataIPDADMA ,inits=NULL,parameters.to.save = c("delta",'u'),model.file =modelIPDADMA,
                                  n.chains=2,n.iter = 10000,n.burnin = 100,DIC=F,n.thin = 1)
#results
IPDADMAJAGSresults

### We can see that u[1:3] are much lower than u[4:5]. In NMA is like all of them drop down so that u[1:3] become less than 0
# and u[4:5] become less than before but still higher than zero


################################################## B. NMR model with only prognostic factor #################################################################
# Here we make the data for jags (for all the models)

source("R/modelIPDADNMRPr.R") #source the model will be used for rjags
####RUN the model
IPDADNMRPrJAGSresults<- jags.parallel(data = jagsdataIPDADNMRPr ,inits=NULL,parameters.to.save = c('delta','u','ORref','gamma'),model.file = modelIPDADNMRPr,
                                       n.chains=2,n.iter = 10000,n.burnin = 100,DIC=F,n.thin = 1)
#traceplots
traceplot(IPDADNMRPrJAGSresults)
#results
IPDADNMRPrJAGSresults
#check the results with IPD NMR model
source('CheckNMRPr.R')
IPDNMRJAGSresultsPreSpecified

################################### C. NMR model with Risk as prognostic factor and only within effect modifier #################################################################


source("R/modelIPDADNMREMwithin.R") #source the model will be used for rjags
####RUN the model
IPDADNMREMwithinJAGSresults <- jags.parallel(data =jagsdataIPDADNMREMwithin ,inits=NULL,parameters.to.save = c('delta','u','ORref','gamma','gamma.w'),model.file = modelIPDADNMREMwithin,
                                       n.chains=2,n.iter = 10000,n.burnin = 100,DIC=F,n.thin = 1)
#traceplot
traceplot(IPDADNMREMwithinJAGSresults)
#results
IPDADNMREMwithinJAGSresults

#check the results with the corresponding IPD NMR model
source('CheckNMREMwithin.R')
IPDNMRJAGSresultsPrespecified
################################### D. NMR model AD and IPD #################################################################

source("R/modelIPDADNMR.R") #source the model will be used for rjags
####RUN the model
IPDADNMRJAGSresults <- jags.parallel(data = jagsdataIPDADNMR ,inits=NULL,parameters.to.save = c('delta','u','ORref','gamma','gamma.b','gamma.w'),model.file = modelIPDADNMR,
                                       n.chains=2,n.iter = 10000,n.burnin = 100,DIC=F,n.thin = 1)
#traceplots
traceplot(IPDADNMRJAGSresults)
#results
IPDADNMRJAGSresults

###################################################################################################################################
#######################                     FIGURES AND TABLES from paper ###############################################
##################################################################################################################################



############################ All tables and figures of presented in the paper #########################
source('Paper - Tables and Figures.R')
#Table 1
AFFIRM
CONFIRM
DEFINE
Placebo


#Table 2
Table_LASSOModel
Table_PreSpecifiedModel

#Table 3

IPDNMR_Table

#Table 4
##predicted probabilities to relapse in 2 years
round(LASSOtable,0) #for LASSO model
round(PreSpecifiedtable,0) #for pre-specified model
###ORs of relapsing in 2 years
round(LASSOtableOR,2) #for LASSO model
round(PreSpecifiedtableOR,2) #for pre-specified model
#Figure 1
ggarrange(PrognosticRiskLASSO,PrognosticRiskPreSpecified,ncol = 1,nrow = 2,labels = c("A  LASSO model","B  Pre-specified model"), hjust=-0.2,font.label = list(size = 11))
#Figure 2
ggarrange(IPDplotLASSO,IPDplotPreSpecified,ncol = 1,nrow = 2,labels = c("A","B"),font.label = list(size = 11))
#Appendix figure 1
ggarrange(IPDplotLASSO_OR,IPDplotPreSpecified_OR,ncol = 1,nrow = 2,labels = c("A","B"),font.label = list(size = 11))
#Appendix figure 3
Flowchart


############################# R-shiny app ##################################################################
source('R-shiny.R') #here you have to change the path of the data
shinyApp(ui = ui, server = server) #also available in https://cinema.ispm.unibe.ch/shinies/koms/

#########################   END    #################################################################
rm(list=ls())
