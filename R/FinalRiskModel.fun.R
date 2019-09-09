###################### Function that returns the final selected model ##############

FinalRiskModel.fun<-function(dataset) {
  
  ##loading libraries
  library(glmnet)
  library(Hmisc)
  library(rms)
  library(glmpath)
  
  ####data
  MSrelapse<-dataset[which(MSrelapse$STUDYID!="ADVANCE"),]
  dataset<-MSrelapse
  #drop STUDYID, USIBJID (no needed for the matrix), drop TRT01A (blinded to treatment)
  ##### drop RELAPSE1year - we need only the analysis for 2 years
  todrop<-c("RELAPSE1year")
  X<-dataset[ , !(names(dataset) %in% todrop)]
  X<-na.omit(X)
  ###logistic regression
  finalmodel<-lrm(formula = RELAPSE2year ~ AGE + WEIGHTBL + EDSSBL + rcs(RLPS3YR,4) + 
                    TRELMOS + PRMSGR + REGION + NHPTMBL + GDLESBL + SFPCSBL + 
                    SENSORBL + DISTWKBL, data = X, x = TRUE, y = TRUE, linear.predictors = T)
  
  # Make a penalized model
  p8	<- pentrace(finalmodel, c(0,1,2,3,4,5,6,7,8,10,12,14,20, 24, 32,40))
  finalmodel.pen <- update (finalmodel, penalty=p8$penalty)
  
  cat("The final model after the shrinkage of coefficients (penalized) is:", fill=TRUE)
  print(finalmodel.pen)
  X$logitRisk<-finalmodel.pen$linear.predictors
  X$Risk<-exp(X$logitRisk)/(1+exp(X$logitRisk))
  
  ##make treatment and study id numeric for the prediction model
  X$TRT01A<-as.numeric(X$TRT01A)-1
  X$TRT01A<-recode(X$TRT01A, "1='1'; 3='2';4='3'; 6='4';")
  X$TRT01A<-as.factor(X$TRT01A)
  X$STUDYID<-as.numeric(X$STUDYID)-1
  X$STUDYID<-as.factor(X$STUDYID)
  
  return(list(finalmodel.pen,X))
  
  
  
}