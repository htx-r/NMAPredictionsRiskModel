############################### SAMPLE SIZE CRRITERIA ####################

#######EPV of candidate variables
dataset<-MSrelapse
todrop<-c("STUDYID","USUBJID","TRT01A","RELAPSE1year")
X<-dataset[ , !(names(dataset) %in% todrop)]
X<-na.omit(X)
fullmodel<-lrm(RELAPSE2year~rcs(AGE,4)+SEX+RACE+rcs(HEIGHTBL,4)+rcs(WEIGHTBL,4)+rcs(EDSSBL,4)+rcs(ONSYRS,4)+DOMIHAND+rcs(RLPS3YR,4)+rcs(TRELMOS,4)+MCDBL+PRMSGR+REGION+rcs(T25FWABL,4)+rcs(NHPTMBL,4)+rcs(PASATABL,4)+rcs(MSFCBL,4)+GDLESBL+rcs(T2VOLBL,4)+rcs(BVZBL,4)+rcs(VFT25BL,4)+rcs(SFPCSBL,4)+rcs(SFMCSBL,4)+VISUALBL+BRAINBL+PYRAMIBL+SENSORBL+BOWLBLBL+CEREBRBL+DISTWKBL+rcs(T25FWPC,4)+rcs(NHPTMPC,4)+rcs(PASATPC,4),x=TRUE,y=TRUE,data=X)
anova(fullmodel)
fullmodel<-lrm(RELAPSE2year~AGE+SEX+RACE+HEIGHTBL+WEIGHTBL+EDSSBL+ONSYRS+DOMIHAND+rcs(RLPS3YR,4)+TRELMOS+MCDBL+PRMSGR+REGION+T25FWABL+NHPTMBL+PASATABL+MSFCBL+GDLESBL+T2VOLBL+rcs(BVZBL,4)+VFT25BL+SFPCSBL+SFMCSBL+VISUALBL+BRAINBL+PYRAMIBL+SENSORBL+BOWLBLBL+CEREBRBL+DISTWKBL+T25FWPC+NHPTMPC+PASATPC,x=TRUE,y=TRUE,data=X)
df<-fullmodel[["stats"]][["d.f."]]
events<- nrow(X[which(X$RELAPSE2year==1),])
EPV<-events/df
cat("The EPV of the model is", EPV, fill=TRUE)

#### sample size by Riley et al
###validate command
#### bootstrap validation for the penalized model
set.seed(1)
val.pen<-validate(RiskModel,method="boot",B=500)
R2Nagel=val.pen[2,5]#that gives the
#R2 Nagelkerdes = 0.0842 (corrected for optimism)

#null model
mod0 <- lrm(RELAPSE2year~1,x=TRUE,y=TRUE,data=X)
MaxRcs<-(1-exp(2*as.numeric(logLik(mod0)/1997)))
##R2csadj=R2Nagelkerdes*Max(R2xsadj)
R2csadj=R2Nagel*MaxRcs

sample_size<-pmsampsize(type="b", rsquared=0.11, shrinkage=0.90, parameters = df, prevalence= 741/1997)
cat("The needed sample size is")
print(sample_size$results_table)
