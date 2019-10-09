###################### Script that returns the final selected model   ############################
###################### and the dataset with the predicted risk for each individual ######################                                                                ##############



  ##loading libraries
  #library(glmnet)
  #library(Hmisc)
  #library(rms)
  #library(glmpath)

  ####data
  dataset<-MSrelapse[which(MSrelapse$STUDYID!="ADVANCE"),]
  #drop STUDYID, USIBJID (no needed for the matrix), drop TRT01A (blinded to treatment)
  ##### drop RELAPSE1year and NAs - we need only the analysis for 2 years
  todrop<-c("RELAPSE1year")
  X<-dataset[ , !(names(dataset) %in% todrop)]
  X<-na.omit(X)

  ###Final model - Model 1 with Penalized Maximum Likelihood Estimation
  finalmodel<-lrm(formula = RELAPSE2year ~ AGE + WEIGHTBL + EDSSBL + RLPS3YR +
                    TRELMOS + PRMSGR + REGION + NHPTMBL + GDLESBL + SFPCSBL +
                    DISTWKBL, data = X, x = TRUE, y = TRUE, linear.predictors = T)
  set.seed(1)
  penalized	<- pentrace(finalmodel, seq(0,200,0.1)) #28.1
  finalmodel.pen <- update (finalmodel, penalty=penalized$penalty)

  cat("The final model after the shrinkage of coefficients (penalized) is:", fill=TRUE)
  print(finalmodel.pen)

  ############################# Missing values for MSCRG study ######################
  MSCRG<-dataset[dataset$STUDYID=="MSCRG",]
  todrop<-c("RELAPSE1year","MCDBL","BVZBL","VFT25BL","SFMCSBL","CEREBRBL","T25FWPC","NHPTMPC","PASATPC","BVTBL")
  MSCRG<-MSCRG[ , !(names(dataset) %in% todrop)]
  ### input NA values to TRELMOS - integer
  todrop<-c("TRT01A","DISTWKBL","SFPCSBL","USUBJID","STUDYID","RELAPSE1year","MCDBL","BVZBL","VFT25BL","SFMCSBL","CEREBRBL","T25FWPC","NHPTMPC","PASATPC","BVTBL")
  missingTREL<-MSrelapse[, !(names(dataset) %in% todrop)]
  missingTREL<-na.omit(missingTREL)
  TRELMOSmodel<-glm(TRELMOS~., data=missingTREL)
  MSCRG$TRELMOS<-predict(TRELMOSmodel,new=MSCRG)
  ### input NA values to SFPCSBL - numeric
  todrop<-c("TRT01A","DISTWKBL","TRELMOS","USUBJID","STUDYID","RELAPSE1year","MCDBL","BVZBL","VFT25BL","SFMCSBL","CEREBRBL","T25FWPC","NHPTMPC","PASATPC","BVTBL")
  missingSF<-MSrelapse[, !(names(dataset) %in% todrop)]
  missingSF<-na.omit(missingSF)
  SFPCSBLmodel<-glm(SFPCSBL~., data=missingSF)
  MSCRG$SFPCSBL<-predict(SFPCSBLmodel,new=MSCRG)
  ### input NA values to DISTWKBL - numeric
  todrop<-c("TRT01A","SFPCSBL","TRELMOS","USUBJID","STUDYID","RELAPSE1year","MCDBL","BVZBL","VFT25BL","SFMCSBL","CEREBRBL","T25FWPC","NHPTMPC","PASATPC","BVTBL")
  missingDIS<-MSrelapse[, !(names(dataset) %in% todrop)]
  missingDIS<-na.omit(missingDIS)
  DISTWKBLmodel<-glm(DISTWKBL~.,family="binomial", data=missingDIS)
  MSCRG$DISTWKBL<-predict(DISTWKBLmodel,new=MSCRG,type="response")

  ############################ Final dataset including risk #############

  FinalDataset<-MSrelapse[which(MSrelapse$STUDYID!="ADVANCE"),]
  todrop<-c("RELAPSE1year")
  FinalDataset<-FinalDataset[ , !(names(dataset) %in% todrop)]
  FinalDataset$TRELMOS[FinalDataset$STUDYID=="MSCRG"]<-MSCRG$TRELMOS
  FinalDataset$SFPCSBL[FinalDataset$STUDYID=="MSCRG"]<-MSCRG$SFPCSBL
  FinalDataset$DISTWKBL[FinalDataset$STUDYID=="MSCRG"]<-1
  tokeep<-c("STUDYID","USUBJID","TRT01A","RELAPSE2year","AGE", "WEIGHTBL","EDSSBL","RLPS3YR",
            "TRELMOS","PRMSGR","REGION","NHPTMBL","GDLESBL","SFPCSBL","DISTWKBL")

  FinalDataset<-FinalDataset[ , (names(FinalDataset) %in% tokeep)]
  X1<-na.omit(FinalDataset)
  ###logit risk and risk for each one of the patients based one the model
  X1$logitRisk<-predict(finalmodel.pen,new=X1)
  X1$Risk<-exp(X1$logitRisk)/(1+exp(X1$logitRisk))

  ##make treatment and study id numeric for the prediction model
  X1$TRT01A<-as.numeric(X1$TRT01A)
  X1$TRT01A<-recode(X1$TRT01A, "1='1'; 2='2';4='3'; 5='4'; 7='5'")
  ###1=Avonex, 2= Dimethyl fumarate, 3=Glatiramer acetate , 4=Natalizumab, 5=Placebo
  X1$TRT01A<-as.factor(X1$TRT01A)
  #recode of studyid
  X1$STUDYID<-as.numeric(X1$STUDYID)-1
  X1$STUDYID<-recode(X1$STUDYID, "1='1'; 2='2';3='3'; 5='4'")
  ## 1=DEFINE, 2=CONFIRM, 3=AFFIRM, 4=MSCRG
  X1$STUDYID<-as.factor(X1$STUDYID)
  RiskData<-X1
  RiskModel<-finalmodel.pen
###remove no needed items
  rm(X)
  rm(X1)
  rm(penalized)
  rm(finalmodel)
  rm(finalmodel.pen)
  rm(dataset)
  rm(DISTWKBLmodel)
  rm(SFPCSBLmodel)
  rm(TRELMOSmodel)
  rm(FinalDataset)
  rm(missingSF)
  rm(missingDIS)
  rm(missingTREL)
  rm(MSCRG)
  rm(tokeep)
