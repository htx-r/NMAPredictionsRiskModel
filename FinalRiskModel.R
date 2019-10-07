###################### Script that returns the final selected model ##############



  ##loading libraries
  #library(glmnet)
  #library(Hmisc)
  #library(rms)
  #library(glmpath)

  ####data
  dataset<-MSrelapse[which(MSrelapse$STUDYID!="ADVANCE"),]
  #drop STUDYID, USIBJID (no needed for the matrix), drop TRT01A (blinded to treatment)
  ##### drop RELAPSE1year - we need only the analysis for 2 years
  todrop<-c("RELAPSE1year")
  X<-dataset[ , !(names(dataset) %in% todrop)]
  X<-na.omit(X)

  ###logistic regression
  finalmodel<-lrm(formula = RELAPSE2year ~ AGE + WEIGHTBL + EDSSBL + RLPS3YR +
                    TRELMOS + PRMSGR + REGION + NHPTMBL + GDLESBL + SFPCSBL +
                    DISTWKBL, data = X, x = TRUE, y = TRUE, linear.predictors = T)


  set.seed(1)
  penalized	<- pentrace(finalmodel, seq(0,200,0.1)) #28.1
  penalized$penalty

  #penalized	<- pentrace(finalmodel, c(0,1,2,3,4,5,6,7,8,10,12,14,20, 24, 32,40))
  finalmodel.pen <- update (finalmodel, penalty=penalized$penalty)
  #### bootstrap validation for the penalized model
  set.seed(1)
  val.pen<-validate(finalmodel.pen,method="boot",B=500)
  c_index.pen <- abs(val.pen[1,5])/2 + 0.5
  c_slope.pen<-val.pen[4,5]


  cat("The final model after the shrinkage of coefficients (penalized) is:", fill=TRUE)
  print(finalmodel.pen)
  X$logitRisk<-predict(finalmodel.pen,newx=X)
  X$Risk<-exp(X$logitRisk)/(1+exp(X$logitRisk))

  ##make treatment and study id numeric for the prediction model
  X$TRT01A<-as.numeric(X$TRT01A)-1
  X$TRT01A<-recode(X$TRT01A, "1='1'; 3='2';4='3'; 6='4'")
  ###1=Dimethyl fumerate, 2=Glatiramer acetate, 3=Natalizumab, 4=Placebo
  X$TRT01A<-as.factor(X$TRT01A)
  X$STUDYID<-as.numeric(X$STUDYID)-1
  ## 1=DEFINE, 2=CONFIRM, 3=AFFIRM
  X$STUDYID<-as.factor(X$STUDYID)
  RiskData<-X
  RiskModel<-finalmodel.pen
###remove no needed items
  rm(X)
  rm(penalized)
  rm(finalmodel)
  rm(finalmodel.pen)
  rm(dataset)


