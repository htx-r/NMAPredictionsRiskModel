###############################################################################
################      FUNCTION FOR CHECKING       #########################
######## 5 DIFFERENT MODELS AND THEIR DISCRIMINATION ABILITY #######

RiskModelSelection.fun<-function(dataset,model){

  ##loading libraries
  #library(glmnet)
  #library(Hmisc)
  #library(rms)
  #library(glmpath)

  if (model=="Internal") {

    ################# LASSO preperation ##########################################
    ################################################################
    dataset<-dataset[which(dataset$STUDYID!="ADVANCE"),]

    #drop STUDYID, USIBJID (no needed for the matrix), drop TRT01A (blinded to treatment)
    ##### drop RELAPSE1year - we need only the analysis for 2 years
    todrop<-c("STUDYID","USUBJID","TRT01A","RELAPSE1year")
    X<-dataset[ , !(names(dataset) %in% todrop)]
    ### develop the matrix for the model
    X.both<-model.matrix(X$RELAPSE2year~.-1,data=X)
    #### delete NAs values (LASSO requirement)
    X.both<-na.omit(X.both)
    X<-na.omit(X)
    set.seed(1)
    #############################LASSO###########################################
    ########################################################################
    ###LASSO (alpha=1-default) with 10 cross validations
    cv.fit.both<-cv.glmnet(x=X.both,y=X$RELAPSE2year,family="binomial",type.measure = "auc")
    ###coefficients of variables via LASSO
    cv.coef.both<-coef(cv.fit.both,s="lambda.1se")
    ###############RESULTS
    ### Non zero coefficients -> selected variables
    cv.pf.em.both<-rownames(cv.coef.both)[as.numeric(cv.coef.both)!=0]
    ###selected variables
    cv.pf.em.both
    ##plot lambda
    plot(cv.fit.both)
    ##plot coeff
    plot(coef(cv.fit.both))

  #A. Method of shrinkage. Model with LASSO selected variables and coefficients
    ###logistic model based on LASSO selected variables
    finalmodel<-lrm(RELAPSE2year~AGE+WEIGHTBL+EDSSBL+RLPS3YR+TRELMOS+PRMSGR+REGION+NHPTMBL+GDLESBL+SFPCSBL+DISTWKBL,x=TRUE,y=TRUE,linear.predictors = T,data=X)
    finalmodel
    lassomodel<-finalmodel
    #lasso coefficients
    lassocoef<-as.matrix(cv.coef.both)
    ##insert lasso coefficients to the model
    lassomodel$coefficients[1]<-lassocoef[1]
    lassomodel$coefficients[2]<-lassocoef[2]
    lassomodel$coefficients[3]<-lassocoef[7]
    lassomodel$coefficients[4]<-lassocoef[8]
    lassomodel$coefficients[5]<-lassocoef[12]
    lassomodel$coefficients[6]<-lassocoef[13]
    lassomodel$coefficients[7]<-lassocoef[16]
    lassomodel$coefficients[8]<-lassocoef[17]
    lassomodel$coefficients[9]<-lassocoef[18]
    lassomodel$coefficients[10]<-lassocoef[19]
    lassomodel$coefficients[11]<-lassocoef[20]
    lassomodel$coefficients[12]<-lassocoef[22]
    lassomodel$coefficients[13]<-lassocoef[25]
    lassomodel$coefficients[14]<-lassocoef[29]
    lassomodel$coefficients[15]<-lassocoef[47]

    ##validate the model (we do not need SE's for the validation)
    set.seed(1)
    val.lassomodel<-validate(lassomodel,method="boot",B=500) #### For validation we do not need SE's
    ##also for the predictions we do not need SE's
    c_index.lassomodel <- abs(val.lassomodel[1,5])/2 + 0.5
    c_slope.lassomodel<-val.lassomodel[4,5]
 # B. Method of shrinkage. LASSO selected variables and then only in the selected variables, estimate again the LASSO penalty
    finalmodel<-lrm(RELAPSE2year~AGE+WEIGHTBL+EDSSBL+RLPS3YR+TRELMOS+PRMSGR+REGION+NHPTMBL+GDLESBL+SFPCSBL+DISTWKBL,x=TRUE,y=TRUE,linear.predictors = T,data=X)
    set.seed(123)
    cv <- cv.glmnet(x=finalmodel$x, y=finalmodel$y, alpha=1, family=c("binomial"), type.measure = "auc")
    # Display the best lambda value
    cv$lambda.1se# 58.68
    # Fit the model on the training data
    model.L1 <- glmnet(x=finalmodel$x, y=finalmodel$y, alpha = 1, lambda = cv$lambda.1se, family=c("binomial"))
    # Display regression coefficients
    coef(model.L1)
    finalmodel.lasso2<-finalmodel
    finalmodel.lasso2$coefficients<-as.numeric(coef(model.L1))
    set.seed(1)
    val.lasso2<-validate(finalmodel.lasso2,method="boot",B=500)
    c_index.lasso2 <- abs(val.lasso2[1,5])/2 + 0.5
    c_slope.lasso2<-val.lasso2[4,5]

 # C. Method of shrinkage. LASSO selected variables in a logistic model using Uniform shrinkage with heuristic factor
    ###logistic model based on LASSO selected variables
    finalmodel<-lrm(RELAPSE2year~AGE+WEIGHTBL+EDSSBL+RLPS3YR+TRELMOS+PRMSGR+REGION+NHPTMBL+GDLESBL+SFPCSBL+DISTWKBL,x=TRUE,y=TRUE,linear.predictors = T,data=X)
    # Make shrunk models
    finalmodel.uniform <- finalmodel
    shrinkage <- (finalmodel$stats[3]-finalmodel$stats[4])/finalmodel$stats[3]
    finalmodel.uniform$coef <- shrinkage * finalmodel$coef
    # heuristic estimate (LR - df) / LR

    # Estimate new intercept, with shrunk lp as offset variable, i.e. coef fixed at unity
    finalmodel.uniform$coef[1] <- lrm.fit(y=finalmodel$y, offset= finalmodel$x %*% finalmodel.uniform$coefficients[2:15])$coefficients[1]
    finalmodel.uniform$coef /finalmodel$coef
    signif(finalmodel.uniform$coef,2)


    #### bootstrap validation for the uniform shrunk model
    set.seed(1)
    val.uniform<-validate(finalmodel.uniform,method="boot",B=500)
    c_index.uniform <- abs(val.uniform[1,5])/2 + 0.5
    c_slope.uniform<-val.uniform[4,5]


# D. Method of shrinkage, ridge with penalty ?? that maximizes AIC modified
    set.seed(1)
    penalized	<- pentrace(finalmodel, seq(0,200,0.1)) #28.1
    penalized$penalty
    plot(pentrace(finalmodel, seq(0,200,0.1)))
    #penalized	<- pentrace(finalmodel, c(0,1,2,3,4,5,6,7,8,10,12,14,20, 24, 32,40))
    finalmodel.pen <- update (finalmodel, penalty=penalized$penalty)
    #### bootstrap validation for the penalized model
    set.seed(1)
    val.pen<-validate(finalmodel.pen,method="boot",B=500)
    c_index.pen <- abs(val.pen[1,5])/2 + 0.5
    c_slope.pen<-val.pen[4,5]

 #E. Method shrinkage. Ridge regression to selected via LASSO coefficients
    ####Ridge shrinkage
    set.seed(123)
    cv <- cv.glmnet(x=finalmodel$x, y=finalmodel$y, alpha=0, family=c("binomial"), type.measure = "auc")
    # Display the best lambda value
    cv$lambda.1se# 58.68
    # Fit the model on the training data
    model.L2 <- glmnet(x=finalmodel$x, y=finalmodel$y, alpha = 0, lambda = cv$lambda.1se, family=c("binomial"))
    # Display regression coefficients
    coef(model.L2)
    finalmodel.ridge<-finalmodel
    finalmodel.ridge$coefficients<-coef(model.L2)
    set.seed(1)
    val.ridge<-validate(finalmodel.ridge,method="boot",B=500)
    c_index.ridge <- abs(val.ridge[1,5])/2 + 0.5
    c_slope.ridge<-val.ridge[4,5]

    ###return
    cat("The selected variables are:" , cv.pf.em.both, fill = TRUE)
    cat("The bootstap corrected discrimination of the LASSO model is", c_index.lassomodel, "and the bootstap corrected c-slope is", c_slope.lassomodel, fill=TRUE)
    cat("The bootstap corrected discrimination after the shrinkage of coefficients (LASSO) of the model is", c_index.lasso2, "and the bootstap corrected c-slope is", c_slope.lasso2, fill=TRUE)
    cat("The bootstap corrected discrimination after the shrinkage of coefficients (uniform) of the model is", c_index.uniform, "and the bootstap corrected c-slope is", c_slope.uniform,  ". The heuristic penalty factor that was used is", shrinkage, fill=TRUE)
    cat("The bootstap corrected discrimination after ridge penalization of coefficients that maximizes AIC of the model is", c_index.pen, "and the bootstap corrected c-slope is", c_slope.pen, ". The penalty that was used is", penalized$penalty, fill=TRUE)
    cat("The bootstap corrected discrimination after ridge penalization of coefficients that maximizes AUC of the model is", c_index.ridge, "and the bootstap corrected c-slope is", c_slope.ridge, fill=TRUE)

    discrimination<-c(c_index.lassomodel,c_index.lasso2,c_index.uniform,c_index.pen,c_index.ridge)
    calibration<-c(c_slope.lassomodel,c_slope.lasso2,c_slope.uniform,c_slope.pen,c_slope.ridge)
    return(list(discrimination=discrimination, calibration=calibration))

  }


  if (model=="InternalSplinesSign") {

    ####data
   dataset<-dataset[which(dataset$STUDYID!="ADVANCE"),]

    #drop STUDYID, USIBJID (no needed for the matrix), drop TRT01A (blinded to treatment)
    ##### drop RELAPSE1year - we need only the analysis for 2 years
    todrop<-c("STUDYID","USUBJID","TRT01A","RELAPSE1year")
    X<-dataset[ , !(names(dataset) %in% todrop)]
    X<-na.omit(X)
    X.both<-model.matrix(X$RELAPSE2year~.-1,data=X)

    #### check the linearity assumption
    finalmodelcheck<-lrm(formula = RELAPSE2year ~ rcs(AGE,4) + rcs(WEIGHTBL,4) + rcs(EDSSBL,4) + rcs(RLPS3YR,4) +
                     rcs(TRELMOS,4) + PRMSGR + REGION + rcs(NHPTMBL,4) + GDLESBL + rcs(SFPCSBL,4) +
                       DISTWKBL, data = X, x = TRUE, y = TRUE, linear.predictors = T)
    anova(finalmodelcheck)
    #### final model based on linearity check
    finalmodel<-lrm(formula = RELAPSE2year ~ AGE + WEIGHTBL + EDSSBL + rcs(RLPS3YR,4) +
                      TRELMOS + PRMSGR + REGION + NHPTMBL + GDLESBL + SFPCSBL +
                      DISTWKBL, data = X, x = TRUE, y = TRUE, linear.predictors = T)


##A. method of shrinkage. LASSO shrinkage penalty only to the selected variables + non-linear variables of the model
    set.seed(123)
    cv <- cv.glmnet(x=finalmodel$x, y=finalmodel$y, alpha=1, family=c("binomial"), type.measure = "auc")
    # Display the best lambda value
    cv$lambda.1se #0.020
    # Fit the model on the training data
    model.L1 <- glmnet(x=finalmodel$x, y=finalmodel$y, alpha = 1, lambda = cv$lambda.1se, family=c("binomial"))
    # Display regression coefficients
    coef(model.L1)
    finalmodel.lasso2<-finalmodel
    finalmodel.lasso2$coefficients<-as.numeric(coef(model.L1))
    val.lasso2<-validate(finalmodel.lasso2,method="boot",B=500)
    c_index.lasso2 <- abs(val.lasso2[1,5])/2 + 0.5
    c_slope.lasso2<-val.lasso2[4,5]
##B. method of shrinkage. Uniform shrinkage method
    # Make shrunk models

    finalmodel.uniform <- finalmodel
    shrinkage <- (finalmodel$stats[3]-finalmodel$stats[4])/finalmodel$stats[3]
    finalmodel.uniform$coef <- shrinkage * finalmodel$coef
    #heuristic estimate (LR - df) / LR =0.90

    # Estimate new intercept, with shrunk lp as offset variable, i.e. coef fixed at unity
    finalmodel.uniform$coef[1] <- lrm.fit(y=finalmodel$y, offset= finalmodel$x %*% finalmodel.uniform$coefficients[2:17])$coefficients[1]
    finalmodel.uniform$coef /finalmodel$coef
    signif(finalmodel.uniform$coef,2)
    finalmodel.uniform

    #### bootstrap validation for the shrunk model
    set.seed(1)
    val.uniform<-validate(finalmodel.uniform,method="boot",B=500)
    c_index.uniform <- abs(val.uniform[1,5])/2 + 0.5
    c_slope.uniform<-val.uniform[4,5]
###C. Method of shrinkage. Make a ridge penalized model based on a modified AIC criterion
    set.seed(1)
    penalized	<- pentrace(finalmodel, seq(0,200,0.1)) #0.1
    finalmodel.pen <- update (finalmodel, penalty=penalized$penalty)
    #### bootstrap validation for the penalized model
    set.seed(1)
    val.pen<-validate(finalmodel.pen,method="boot",B=500)
    c_index.pen <- abs(val.pen[1,5])/2 + 0.5
    c_slope.pen<-val.pen[4,5]


## D. Method of shrinkage. Ridge shrinkage based on AUC for optimal penalty ??
    set.seed(123)
    cv <- cv.glmnet(x=finalmodel$x, y=finalmodel$y, alpha=0, family=c("binomial"), type.measure = "auc")
    # Display the best lambda value
    cv$lambda.1se # 61.98
    # Fit the model on the training data
    model.L2 <- glmnet(x=finalmodel$x, y=finalmodel$y, alpha = 0, lambda = cv$lambda.1se, family=c("binomial"))

    # Display regression coefficients
    coef(model.L2)
    finalmodel.ridge<-finalmodel
    finalmodel.ridge$coefficients<-as.numeric(coef(model.L2))
    set.seed(1)
    val.ridge<-validate(finalmodel.ridge,method="boot",B=500)
    c_index.ridge <- abs(val.ridge[1,5])/2 + 0.5
    c_slope.ridge<-val.ridge[4,5]


    cat("The bootstap corrected discrimination of the model after LASSO shrinkage to selected variables only is", c_index.lasso2, "and the bootstap corrected c-slope is", c_slope.lasso2, fill=TRUE)
    cat("The bootstap corrected discrimination after the shrinkage of coefficients (uniform) of the model is", c_index.uniform, "and the bootstap corrected c-slope is", c_slope.uniform, ". The heuristic factora s that was used is",shrinkage,fill=TRUE)
    cat("The bootstap corrected discrimination after ridge penalization of coefficients that maximizes AIC of the model is", c_index.pen, "and the bootstap corrected c-slope is", c_slope.pen, ". The penalty that was used is", penalized$penalty, fill=TRUE)
    cat("The bootstap corrected discrimination after ridge penalization of coefficients that maximizes AUC of the model is", c_index.ridge, "and the bootstap corrected c-slope is", c_slope.ridge, fill=TRUE)

    discrimination<-c(c_index.lasso2,c_index.uniform,c_index.pen,c_index.ridge)
    calibration<-c(c_slope.lasso2,c_slope.uniform,c_slope.pen,c_slope.ridge)
    return(list(discrimination=discrimination, calibration=calibration))


  }

  if (model=="InternalInteractions") {

    ################# DATA ##########################################
    ################################################################
    dataset<-dataset[which(dataset$STUDYID!="ADVANCE"),]

    ####keep only the variables from internal LASSO selection & the outcome
    tokeep<-c("AGE","WEIGHTBL","EDSSBL","RLPS3YR","TRELMOS","PRMSGR","REGION","NHPTMBL","GDLESBL","SFPCSBL","DISTWKBL","RELAPSE2year")
    X<-dataset[ , (names(dataset) %in% tokeep)]
    ### develop the matrix for the model
    X.both<-model.matrix(X$RELAPSE2year~.^2,data=X)
    #### delete NAs values (LASSO requirement)
    X.both<-na.omit(X.both[,-1])
    X<-na.omit(X)
    set.seed(1)
    #############################LASSO###########################################
    ########################################################################
    ###LASSO (alpha=1-default) with 10 cross validations
    cv.fit.both<-cv.glmnet(x=X.both,y=X$RELAPSE2year,family="binomial")
    ###coefficients of variables via LASSO
    cv.coef.both<-coef(cv.fit.both,s="lambda.1se")
    ###############RESULTS
    ### Non zero coefficients -> selected variables

    cv.pf.em.both<-rownames(cv.coef.both)[as.numeric(cv.coef.both)!=0]

    ###selected variables
    cv.pf.em.both
    ##plot lambda
    plot(cv.fit.both)
    ##plot coeff
    plot(coef(cv.fit.both))
    ###logistic model based on LASSO selected variables
    finalmodel<-lrm(RELAPSE2year~AGE+WEIGHTBL+EDSSBL+RLPS3YR+TRELMOS+PRMSGR+REGION+NHPTMBL+GDLESBL+SFPCSBL+DISTWKBL+AGE*SFPCSBL+WEIGHTBL*TRELMOS+WEIGHTBL*SFPCSBL+EDSSBL*RLPS3YR+EDSSBL*PRMSGR+EDSSBL*REGION+EDSSBL*GDLESBL+PRMSGR*REGION+REGION*NHPTMBL+REGION*GDLESBL+SFPCSBL*DISTWKBL,linear.predictors = T,x=TRUE,y=TRUE,data=X)

    #################### Make Shrunk coefficients ##############################
##A. method of shrinkage. LASSO shrinkage penalty only to the selected variables + non-linear variables of the model
    set.seed(123)
    cv <- cv.glmnet(x=finalmodel$x, y=finalmodel$y, alpha=1, family=c("binomial"), type.measure = "auc")
    # Display the best lambda value
    cv$lambda.1se #0.044
    # Fit the model on the training data
    model.L1 <- glmnet(x=finalmodel$x, y=finalmodel$y, alpha = 1, lambda = cv$lambda.1se, family=c("binomial"))
    # Display regression coefficients
    coef(model.L1)
    finalmodel.lasso2<-finalmodel
    finalmodel.lasso2$coefficients<-as.numeric(coef(model.L1))
    set.seed(1)
    val.lasso2<-validate(finalmodel.lasso2,method="boot",B=500)
    c_index.lasso2 <- abs(val.lasso2[1,5])/2 + 0.5
    c_slope.lasso2<-val.lasso2[4,5]
##B. Method of shrinkage. Use uniform shrinkage for the selected via LASSO coefficients

    # Make shrunk models
    finalmodel.uniform<- finalmodel
    shrinkage <- (finalmodel$stats[3]-finalmodel$stats[4])/finalmodel$stats[3]
    finalmodel.uniform$coef <- shrinkage * finalmodel$coef
    #heuristic estimate (LR - df) / LR = 0.79

    # Estimate new intercept, with shrunk lp as offset variable, i.e. coef fixed at unity
    finalmodel.uniform$coef[1] <- lrm.fit(y=finalmodel$y, offset= finalmodel$x %*% finalmodel.uniform$coefficients[2:38])$coefficients[1]
    finalmodel.uniform$coef /finalmodel$coef
    signif(finalmodel.uniform$coef,2)
    finalmodel.uniform

    #### bootstrap validation for the shrunk model
    set.seed(1)
    val.uniform<-validate(finalmodel.uniform,method="boot",B=500)
    c_index.uniform <- abs(val.uniform[1,5])/2 + 0.5
    c_slope.uniform<-val.uniform[4,5]
####C. Method of shrinkage. Make a ridge penalized model based on a modified AIC criterion
    # Make a penalized model
    set.seed(1)
    penalized	<- pentrace(finalmodel, seq(0,300,0.1))
    finalmodel.pen <- update (finalmodel, penalty=penalized$penalty) #238.3
    #### bootstrap validation for the penalized model
    set.seed(1)
    val.pen<-validate(finalmodel.pen,method="boot",B=500)
    c_index.pen <- abs(val.pen[1,5])/2 + 0.5
    c_slope.pen<-val.pen[4,5]


## D. Method of shrinkage. Ridge shrinkage based on AUC for optimal penalty ??

    ####Ridge shrinkage
    set.seed(123)
    cv <- cv.glmnet(x=finalmodel$x, y=finalmodel$y, alpha=0, family=c("binomial"), type.measure = "auc")
    # Display the best lambda value
    cv$lambda.1se # 65
    # Fit the model on the training data
    model.L2 <- glmnet(x=finalmodel$x, y=finalmodel$y, alpha = 0, lambda = cv$lambda.1se, family=c("binomial"))

    # Display regression coefficients
    coef(model.L2)
    finalmodel.ridge<-finalmodel
    finalmodel.ridge$coefficients<-as.numeric(coef(model.L2))
    set.seed(1)
    val.ridge<-validate(finalmodel.ridge,method="boot",B=500)
    c_index.ridge <- abs(val.ridge[1,5])/2 + 0.5
    c_slope.ridge<-val.ridge[4,5]


    cat("The bootstap corrected discrimination of the model after LASSO shrinkage to selected variables only is", c_index.lasso2, "and the bootstap corrected c-slope is", c_slope.lasso2, fill=TRUE)
    cat("The bootstap corrected discrimination after the shrinkage of coefficients (uniform) of the model is", c_index.uniform, "and the bootstap corrected c-slope is", c_slope.uniform, ". The heuristic factor s that was used is", shrinkage,  fill=TRUE)
    cat("The bootstap corrected discrimination after ridge penalization of coefficients that maximizes AIC of the model is", c_index.pen, "and the bootstap corrected c-slope is", c_slope.pen, ". The penalty that was used is", penalized$penalty, fill=TRUE)
    cat("The bootstap corrected discrimination after ridge penalization of coefficients that maximizes AUC of the model is", c_index.ridge, "and the bootstap corrected c-slope is", c_slope.ridge,fill=TRUE)

    discrimination<-c(c_index.lasso2,c_index.uniform,c_index.pen,c_index.ridge)
    calibration<-c(c_slope.lasso2,c_slope.uniform,c_slope.pen,c_slope.ridge)
    return(list(discrimination=discrimination, calibration=calibration))

  }

  if (model=="Fabio") {
    ###data
   dataset<-dataset[which(dataset$STUDYID!="ADVANCE"),]

    #drop STUDYID, USIBJID (no needed for the matrix), drop TRT01A (blinded to treatment)
    ##### drop RELAPSE1year - we need only the analysis for 2 years
    todrop<-c("STUDYID","USUBJID","TRT01A","RELAPSE1year")
    X<-dataset[ , !(names(dataset) %in% todrop)]
    X<-na.omit(X)
    X.both<-model.matrix(X$RELAPSE2year~.-1,data=X)
    finalmodel<-lrm(RELAPSE2year~AGE+SEX+RACE+EDSSBL+ONSYRS+RLPS3YR+TRELMOS+PRMSGR+T25FWABL+NHPTMBL+PASATABL+VFT25BL+SFPCSBL+SFMCSBL,x=TRUE,y=TRUE,linear.predictors = TRUE,data=X)
    finalmodel

    #################### Make Shrunk coefficients ##############################
##A. method of shrinkage. LASSO shrinkage penalty only to the selected variables + non-linear variables of the model
    set.seed(123)
    cv <- cv.glmnet(x=finalmodel$x, y=finalmodel$y, alpha=1, family=c("binomial"), type.measure = "auc")
    # Display the best lambda value
    cv$lambda.1se #0.020
    # Fit the model on the training data
    model.L1 <- glmnet(x=finalmodel$x, y=finalmodel$y, alpha = 1, lambda = cv$lambda.1se, family=c("binomial"))
    # Display regression coefficients
    coef(model.L1)
    finalmodel.lasso2<-finalmodel
    finalmodel.lasso2$coefficients<-as.numeric(coef(model.L1))
    set.seed(1)
    val.lasso2<-validate(finalmodel.lasso2,method="boot",B=500)
    c_index.lasso2 <- abs(val.lasso2[1,5])/2 + 0.5
    c_slope.lasso2<-val.lasso2[4,5]
##B. Method of shrinkage. Use uniform shrinkage for the selected via LASSO coefficients

    # Make shrunk models
    finalmodel.uniform <- finalmodel
    shrinkage <- (finalmodel$stats[3]-finalmodel$stats[4])/finalmodel$stats[3]
    finalmodel.uniform$coef <- shrinkage * finalmodel$coef
    #  heuristic estimate (LR - df) / LR = 0.92

    # Estimate new intercept, with shrunk lp as offset variable, i.e. coef fixed at unity
    finalmodel.uniform$coef[1] <- lrm.fit(y=finalmodel$y, offset= finalmodel$x %*% finalmodel.uniform$coefficients[2:15])$coefficients[1]
    finalmodel.uniform$coef /finalmodel$coef
    signif(finalmodel.uniform$coef,2)
    finalmodel.uniform

    #### bootstrap validation for the shrunk model
    set.seed(1)
    val.uniform<-validate(finalmodel.uniform,method="boot",B=500)
    c_index.uniform <- abs(val.uniform[1,5])/2 + 0.5
    c_slope.uniform<-val.uniform[4,5]
###C. Method of shrinkage. Make a ridge penalized model based on a modified AIC criterion
    # Make a penalized model
    set.seed(1)
    penalized	<-  pentrace(finalmodel, seq(0,200,0.1)) #44.4
    finalmodel.pen <- update (finalmodel, penalty=penalized$penalty)
    #### bootstrap validation for the penalized model
    set.seed(1)
    val.pen<-validate(finalmodel.pen,method="boot",B=500)
    c_index.pen <- abs(val.pen[1,5])/2 + 0.5
    c_slope.pen<-val.pen[4,5]
## D. Method of shrinkage. Ridge shrinkage based on AUC for optimal penalty ??

    ####Ridge shrinkage
    set.seed(123)
    cv <- cv.glmnet(x=finalmodel$x, y=finalmodel$y, alpha=0, family=c("binomial"), type.measure = "auc")
    # Display the best lambda value
    cv$lambda.1se # 6.77
    # Fit the model on the training data
    model.L2 <- glmnet(x=finalmodel$x, y=finalmodel$y, alpha = 0, lambda = cv$lambda.1se, family=c("binomial"))

    # Display regression coefficients
    coef(model.L2)
    finalmodel.ridge<-finalmodel
    finalmodel.ridge$coefficients<-as.numeric(coef(model.L2))
    set.seed(1)
    val.ridge<-validate(finalmodel.ridge,method="boot",B=500)
    c_index.ridge <- abs(val.ridge[1,5])/2 + 0.5
    c_slope.ridge<-val.ridge[4,5]


    cat("The bootstap corrected discrimination of the model after LASSO shrinkage to selected variables only is", c_index.lasso2, "and the bootstap corrected c-slope is", c_slope.lasso2, fill=TRUE)
    cat("The bootstap corrected discrimination after the shrinkage of coefficients (uniform) of the model is", c_index.uniform, "and the bootstap corrected c-slope is", c_slope.uniform,". The heuristic factor s that was used is", penalized$penalty, fill=TRUE)
    cat("The bootstap corrected discrimination after ridge penalization of coefficients that maximizes AIC of the model is", c_index.pen, "and the bootstap corrected c-slope is", c_slope.pen, ". The penalty that was used is", penalized$penalty, fill=TRUE)
    cat("The bootstap corrected discrimination after ridge penalization of coefficients that maximizes AUC of the model is", c_index.ridge, "and the bootstap corrected c-slope is", c_slope.ridge,  fill=TRUE)

    discrimination<-c(c_index.lasso2,c_index.uniform,c_index.pen,c_index.ridge)
    calibration<-c(c_slope.lasso2,c_slope.uniform,c_slope.pen,c_slope.ridge)
    return(list(discrimination=discrimination, calibration=calibration))


  }

}




