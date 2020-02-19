#######################################################################################################################
###################################   SCRIPT FOR c-for benefit   #######################################################
########################################################################################################################


# Step 1. I calculate via our NMR model the predicted probabilities to relapse for each patient in our dataset (2000 patients)
# So here we have 2000 patients, their baseline risk of patients prior to treatment
# and their predicted probability to relapse in two years under each treatment based on our NMR model

jagsdataIPDNMR_PreSpecified_concordance <- list(
  Nstudies=3,
  Np=nrow(RiskData),
  studyid=as.numeric(RiskData$STUDYID),
  outcome=as.numeric(RiskData$RELAPSE2year)-1,
  outcomeP=PlaceboArms$Relapse2year,
  NpPlacebo=nrow(PlaceboArms),
  treat= rbind(c(1,4,NA),c(1,2,4),c(3,4,NA)),
  na=c(2,3,2),
  logitRisknew=as.data.frame(RiskData$logitRiskPreSpecified, ncol=1),
  logitmeanRisknew=mean(RiskData$logitRiskPreSpecified),
  Nnew=nrow(RiskData),
  arm=RiskData$arm,
  Risk=RiskData$logitRiskPreSpecified,
  meanRisk=c(tapply(RiskData$logitRiskPreSpecified, RiskData$STUDYID, summary)$`1`[4],tapply(RiskData$logitRiskPreSpecified, RiskData$STUDYID, summary)$`2`[4],tapply(RiskData$logitRiskPreSpecified, RiskData$STUDYID, summary)$`3`[4]), ##here is the mean of logit of risk
  nt=4,
  ref=4
)

#NMR model
IPDNMRJAGSmodel_PreSpecified_concordance <- jags.parallel(data = jagsdataIPDNMR_PreSpecified_concordance,inits=NULL,parameters.to.save = c('logitp'),model.file = modelIPDNMR,
                                      n.chains=2,n.iter = 10000,n.burnin = 1000,DIC=F,n.thin = 10)

##obtain probabilities instead of logit p
p<-IPDNMRJAGSmodel_PreSpecified_concordance$BUGSoutput$mean$logitp
expit<-function(x) {exp(x)/(1+exp(x))}
p<-expit(IPDNMRJAGSmodel_PreSpecified_concordance$BUGSoutput$mean$logitp)
p<-as.data.frame(p) # this is a table with 2000 patients and their probabilities under each one of the treatments
colnames(p)<-c("Dimethyl Fumarate","Glatiramer Acetate","Natalizumab","Placebo")

## Step 2. We add to the previous table the observed outcome for each one of the patients
Relapse2year<-RiskData$RELAPSE2year
Relapse2year<-as.data.frame(Relapse2year)
Data<-cbind(p,Relapse2year)

## Step 3. We add a column with the predicted benefit of Natalizumab vs Placebo
benefit<-as.data.frame(Data$Placebo-Data$Natalizumab) #Predicted benefit of Natalizumab vs Placebo
names(benefit)<-c("Benefit")
FinalData<-cbind(Data,benefit)

### Step 4. We add the column with the actual treatment that each patient was randomized to take
Treatment<-as.data.frame(RiskData$TRT01A)
colnames(Treatment)<-c("Treatment")
FinalData<-cbind(FinalData,Treatment)

### Step 5. Keep only patients under Placebo or Natalizumab (Now we examine only the c-for benefit of Natalizumab)
#FinalData<-FinalData[which(FinalData$Treatment==3  | FinalData$Treatment==4),]
FinalDataNat<-FinalData[which(FinalData$Treatment==3),]
FinalDataPla<-FinalData[which(FinalData$Treatment==4),]

## Step 6. keep equal number of patients in each arm, select randomly from the largest arm as many patients as there are in the
# smallest arm
nrow(FinalDataNat)
nrow(FinalDataPla)
min(nrow(FinalDataNat),nrow(FinalDataNat))
FinalDataPla<-sample_n(FinalDataPla,min(nrow(FinalDataNat),nrow(FinalDataNat)))
# Step 7. Order each arm based on predicted benefits
FinalDataPla<-FinalDataPla[order(FinalDataPla$Benefit),]
FinalDataNat<-FinalDataNat[order(FinalDataNat$Benefit),]

# Step 8. Match each patient in Natalizumab arm with a patient in Placebo arm based on their rank of probabilities
x<-data.frame(FinalDataNat$Benefit,FinalDataNat$Relapse2year,FinalDataPla$Benefit,FinalDataPla$Relapse2year)

#  Step 9. Add two columns: a. The average of their individual benefit predictions
# their observed benefit = outcome of patient under Placebo - Outcome under Natalizumab

pred.ben.avg<-as.data.frame((FinalDataNat$Benefit+FinalDataPla$Benefit)/2)
colnames(pred.ben.avg)<-c("AveragePredictedBenefit")
obs.ben<-as.data.frame(as.numeric(FinalDataPla$Relapse2year)-as.numeric(FinalDataNat$Relapse2year)) # -1 , 0, or 1
colnames(obs.ben)<-c("ObservedBenefit")
x<-as.data.frame(cbind(x,pred.ben.avg,obs.ben))

# Step 10. Calculate the c-for-benefit
cindex <- rcorr.cens(x$AveragePredictedBenefit, x$ObservedBenefit)
c.benefit <- cindex["C Index"][[1]]
c.benefit.se <- cindex["S.D."][[1]]/2	# The sd of the c-index is half the sd of Dxy

c.benefit
# [1] 0.4871		# The c-for-benefit
c.benefit - 1.96*c.benefit.se
# [1] 0.4460425		# The 95% lower bound of the c-for-benefit
c.benefit + 1.96*c.benefit.se
# [1] 0.5282942

