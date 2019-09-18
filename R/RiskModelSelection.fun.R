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
    X.both<-model.matrix(X$RELAPSE2year~.,data=X)
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
    finalmodel<-lrm(RELAPSE2year~AGE+WEIGHTBL+EDSSBL+RLPS3YR+TRELMOS+PRMSGR+REGION+SENSORBL+NHPTMBL+GDLESBL+SFPCSBL+DISTWKBL,x=TRUE,y=TRUE,linear.predictors = T,data=X)
    finalmodel
    val<-validate(finalmodel,method="boot",B=500)
    c_index <- abs(val[1,5])/2 + 0.5
    c_slope<-val[4,5]
    ##check the assumption of linearity
    ##best choice of knots and test for linearity
    #a<-lrm(RELAPSE2year~rcs(AGE,5)+rcs(WEIGHTBL,5)+rcs(EDSSBL,5)+rcs(RLPS3YR,5)+rcs(TRELMOS,5)+PRMSGR+REGION+rcs(NHPTDHBL,5)+rcs(GDLESBL,5)+rcs(SFPCSBL,5)+SENSORBL+DISTWKBL,data=X)
    #anova(a)##relapse3years
    #b<-lrm(RELAPSE2year~rcs(AGE,4)+rcs(WEIGHTBL,4)+rcs(EDSSBL,4)+rcs(RLPS3YR,4)+rcs(TRELMOS,4)+PRMSGR+REGION+rcs(NHPTDHBL,4)+rcs(GDLESBL,4)+rcs(SFPCSBL,4)+SENSORBL+DISTWKBL,data=X)
    #anova(b)##relapse3years non linear
    #AIC(a)
    #AIC(b) ###best choice of knots is knots=4

    #################### Make Shrunk coefficients ##############################

    # Make shrunk models
    finalmodel.shrunk <- finalmodel
    finalmodel.shrunk$coef <- val[4,5] * finalmodel.shrunk$coef  # use result from bootstrapping
    shrinkage <- (finalmodel$stats[3]-finalmodel$stats[4])/finalmodel$stats[3]
    finalmodel.shrunk$coef <- shrinkage * finalmodel$coef
    # val.full[4,5] is shrinkage factor; heuristic estimate (LR - df) / LR =0.89

    # Estimate new intercept, with shrunk lp as offset variable, i.e. coef fixed at unity
    finalmodel.shrunk$coef[1] <- lrm.fit(y=finalmodel$y, offset= finalmodel$x %*% finalmodel.shrunk$coefficients[2:18])$coefficients[1]
    finalmodel.shrunk$coef /finalmodel$coef
    signif(finalmodel.shrunk$coef,2)
    finalmodel.shrunk

    #### bootstrap validation for the shrunk model
    val.shrunk<-validate(finalmodel.shrunk,method="boot",B=500)
    c_index.shrunk <- abs(val.shrunk[1,5])/2 + 0.5
    c_slope.shrunk<-val.shrunk[4,5]
    # Make a penalized model

    penalized	<- pentrace(finalmodel, seq(0,100,0.1)) #37.7

    #penalized	<- pentrace(finalmodel, c(0,1,2,3,4,5,6,7,8,10,12,14,20, 24, 32,40))
    finalmodel.pen <- update (finalmodel, penalty=penalized$penalty)
    #### bootstrap validation for the penalized model
    val.pen<-validate(finalmodel.pen,method="boot",B=500)
    c_index.pen <- abs(val.pen[1,5])/2 + 0.5
    c_slope.pen<-val.pen[4,5]

    ####other penalized way - not sure about that!!!
    set.seed(1)
    cv.penalty<-cv.glmnet(x=finalmodel$x,y=finalmodel$y,family="binomial", alpha=0,type.measure = "auc")
    best.lambdase=cv.penalty$lambda.1se #58.68
    pen2<-pentrace(finalmodel,c(2*best.lambdase))
    finalmodel.pen2 <- update (finalmodel, penalty=pen2$penalty) #117.36

    #### bootstrap validation for the penalized model
    val.pen2<-validate(finalmodel.pen2,method="boot",B=500)
    c_index.pen2 <- abs(val.pen2[1,5])/2 + 0.5
    c_slope.pen2<-val.pen2[4,5]


    ####Ridge shrinkage
    set.seed(123)
    cv <- cv.glmnet(x=finalmodel$x, y=finalmodel$y, alpha=0, family=c("binomial"), type.measure = "auc")
    # Display the best lambda value
    cv$lambda.min # 0.2013
    # Fit the model on the training data
    model.L2 <- glmnet(x=finalmodel$x, y=finalmodel$y, alpha = 0, lambda = cv$lambda.min, family=c("binomial"),)

    # Display regression coefficients
    coef(model.L2)
    finalmodel.ridge<-finalmodel
    finalmodel.ridge$coefficients<-as.numeric(coef(model.L2))
    val.ridge<-validate(finalmodel.ridge,method="boot",B=500)
    c_index.ridge <- abs(val.ridge[1,5])/2 + 0.5
    c_slope.ridge<-val.ridge[4,5]

    #####LASSO shrinkage

    set.seed(123)
    cv <- cv.glmnet(x=finalmodel$x, y=finalmodel$y, family=c("binomial"), type.measure = "auc")
    # Display the best lambda value
    cv$lambda.min # 0.0076
    # Fit the model on the training data
    model.L1 <- glmnet(x=finalmodel$x, y=finalmodel$y,lambda = cv$lambda.min, family=c("binomial"),)

    # Display regression coefficients
    coef(model.L1)
    finalmodel.lasso<-finalmodel
    finalmodel.lasso$coefficients<-as.numeric(coef(model.L1))
    val.lasso<-validate(finalmodel.lasso,method="boot",B=500)
    c_index.lasso <- abs(val.lasso[1,5])/2 + 0.5
    c_slope.lasso<-val.lasso[4,5]


    #######EPV
    df<-finalmodel[["stats"]][["d.f."]]
    events<- nrow(X[which(X$RELAPSE2year==1),])
    ####return back
    cat("The selected variables are:" , cv.pf.em.both, fill = TRUE)
    EPV<-events/df
    cat("The EPV of the model is", EPV, fill=TRUE)
    if (EPV>10){cat("The EPV of the model is >10, as recommended.", fill=TRUE)}
    if (EPV<10){cat("The EPV of the model is <10, possibility of overfitting problem.", fill = TRUE)}
    #cat("The final model is:",fill=TRUE)
    cat("The bootstap corrected discrimination of the model is", c_index, "and the bootstap corrected c-slope is", c_slope, fill=TRUE)
    cat("The bootstap corrected discrimination after the shrinkage of coefficients (uniform) of the model is", c_index.shrunk, "and the bootstap corrected c-slope is", c_slope.shrunk, fill=TRUE)
    cat("The bootstap corrected discrimination after penalization of coefficients (Penalized ML) of the model is", c_index.pen, "and the bootstap corrected c-slope is", c_slope.pen, ". The penalty that was used is", penalized$penalty, fill=TRUE)

    discrimination<-c(c_index,c_index.shrunk,c_index.pen,c_index.pen2,c_index.ridge,c_index.lasso)
    calibration<-c(c_slope,c_slope.shrunk,c_slope.pen,c_slope.pen2,c_slope.ridge,c_slope.lasso)
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
    X.both<-model.matrix(X$RELAPSE2year~.,data=X)
    #### delete NAs values (LASSO requirement)
    X.both<-na.omit(X.both[,-1])

    finalmodel<-lrm(formula = RELAPSE2year ~ AGE + WEIGHTBL + EDSSBL + rcs(RLPS3YR,4) +
                      TRELMOS + PRMSGR + REGION + NHPTMBL + GDLESBL + SFPCSBL +
                      SENSORBL + DISTWKBL, data = X, x = TRUE, y = TRUE, linear.predictors = T)
    set.seed(1)
    val<-validate(finalmodel,method="boot",B=500)
    c_index <- abs(val[1,5])/2 + 0.5
    c_slope<-val[4,5]
    ##check the assumption of linearity
    ##best choice of knots and test for linearity
    #a<-lrm(RELAPSE2year~rcs(AGE,5)+rcs(WEIGHTBL,5)+rcs(EDSSBL,5)+rcs(RLPS3YR,5)+rcs(TRELMOS,5)+PRMSGR+REGION+rcs(NHPTDHBL,5)+rcs(GDLESBL,5)+rcs(SFPCSBL,5)+SENSORBL+DISTWKBL,data=X)
    #anova(a)##relapse3years
    #b<-lrm(RELAPSE2year~rcs(AGE,4)+rcs(WEIGHTBL,4)+rcs(EDSSBL,4)+rcs(RLPS3YR,4)+rcs(TRELMOS,4)+PRMSGR+REGION+rcs(NHPTDHBL,4)+rcs(GDLESBL,4)+rcs(SFPCSBL,4)+SENSORBL+DISTWKBL,data=X)
    anova(finalmodel)##relapse3years non linear
    #AIC(a)
    #AIC(b) ###best choice of knots is knots=4
    #################### Make Shrunk coefficients ##############################

    # Make shrunk models
    finalmodel.shrunk <- finalmodel
    finalmodel.shrunk$coef <- val[4,5] * finalmodel.shrunk$coef  # use result from bootstrapping
    shrinkage <- (finalmodel$stats[3]-finalmodel$stats[4])/finalmodel$stats[3]
    finalmodel.shrunk$coef <- shrinkage * finalmodel$coef
    # val.full[4,5] is shrinkage factor; heuristic estimate (LR - df) / LR =0.89

    # Estimate new intercept, with shrunk lp as offset variable, i.e. coef fixed at unity
    finalmodel.shrunk$coef[1] <- lrm.fit(y=finalmodel$y, offset= finalmodel$x %*% finalmodel.shrunk$coefficients[2:20])$coefficients[1]
    finalmodel.shrunk$coef /finalmodel$coef
    signif(finalmodel.shrunk$coef,2)
    finalmodel.shrunk

    #### bootstrap validation for the shrunk model

    val.shrunk<-validate(finalmodel.shrunk,method="boot",B=500)
    c_index.shrunk <- abs(val.shrunk[1,5])/2 + 0.5
    c_slope.shrunk<-val.shrunk[4,5]
    # Make a penalized model
    penalized	<- pentrace(finalmodel, seq(0,100,0.1)) #0.1
    finalmodel.pen <- update (finalmodel, penalty=penalized$penalty)
    #### bootstrap validation for the penalized model
    set.seed(1)
    val.pen<-validate(finalmodel.pen,method="boot",B=500)
    c_index.pen <- abs(val.pen[1,5])/2 + 0.5
    c_slope.pen<-val.pen[4,5]

    ####other penalized way - not sure about that!!!
    set.seed(1)
    cv.penalty<-cv.glmnet(x=finalmodel$x,y=finalmodel$y,family="binomial", alpha=0,type.measure = "auc")
    best.lambdase=cv.penalty$lambda.1se #58.68
    pen2<-pentrace(finalmodel,c(2*best.lambdase))
    finalmodel.pen2 <- update (finalmodel, penalty=pen2$penalty) #0

    #### bootstrap validation for the penalized model
    val.pen2<-validate(finalmodel.pen2,method="boot",B=500)
    c_index.pen2 <- abs(val.pen2[1,5])/2 + 0.5
    c_slope.pen2<-val.pen2[4,5]


    ####Ridge shrinkage
    set.seed(123)
    cv <- cv.glmnet(x=finalmodel$x, y=finalmodel$y, alpha=0, family=c("binomial"), type.measure = "auc")
    # Display the best lambda value
    cv$lambda.min # 0.2013
    # Fit the model on the training data
    model.L2 <- glmnet(x=finalmodel$x, y=finalmodel$y, alpha = 0, lambda = cv$lambda.min, family=c("binomial"),)

    # Display regression coefficients
    coef(model.L2)
    finalmodel.ridge<-finalmodel
    finalmodel.ridge$coefficients<-as.numeric(coef(model.L2))
    val.ridge<-validate(finalmodel.ridge,method="boot",B=500)
    c_index.ridge <- abs(val.ridge[1,5])/2 + 0.5
    c_slope.ridge<-val.ridge[4,5]

    #####LASSO shrinkage

    set.seed(123)
    cv <- cv.glmnet(x=finalmodel$x, y=finalmodel$y, family=c("binomial"), type.measure = "auc")
    # Display the best lambda value
    cv$lambda.min # 0.0076
    # Fit the model on the training data
    model.L1 <- glmnet(x=finalmodel$x, y=finalmodel$y,lambda = cv$lambda.min, family=c("binomial"),)

    # Display regression coefficients
    coef(model.L1)
    finalmodel.lasso<-finalmodel
    finalmodel.lasso$coefficients<-as.numeric(coef(model.L1))
    val.lasso<-validate(finalmodel.lasso,method="boot",B=500)
    c_index.lasso <- abs(val.lasso[1,5])/2 + 0.5
    c_slope.lasso<-val.lasso[4,5]


    #######EPV

    df<-finalmodel[["stats"]][["d.f."]]
    events<- nrow(X[which(X$RELAPSE2year==1),])
    ####return back
    EPV<-events/df
    cat("The EPV of the model is", EPV, fill=TRUE)
    if (EPV>10){cat("The EPV of the model is >10, as recommended.", fill=TRUE)}
    if (EPV<10){cat("The EPV of the model is <10, possibility of overfitting problem.", fill = TRUE)}
    #cat("The final model is:",fill=TRUE)
    cat("The bootstap corrected discrimination of the model is", c_index, "and the bootstap corrected c-slope is", c_slope, fill=TRUE)
    cat("The bootstap corrected discrimination after the shrinkage of coefficients (uniform) of the model is", c_index.shrunk, fill=TRUE)
    cat("The bootstap corrected discrimination after penalization of coefficients (Penalized ML) of the model is", c_index.pen, "and the bootstap corrected c-slope is", c_slope.pen, ". The penalty that was used is", penalized$penalty, fill=TRUE)
    discrimination<-c(c_index,c_index.shrunk,c_index.pen,c_index.pen2,c_index.ridge,c_index.lasso)
    calibration<-c(c_slope,c_slope.shrunk,c_slope.pen,c_slope.pen2,c_slope.ridge,c_slope.lasso)
    return(list(discrimination=discrimination, calibration=calibration))



  }

  if (model=="InternalInteractions") {

    ################# DATA ##########################################
    ################################################################
    dataset<-dataset[which(dataset$STUDYID!="ADVANCE"),]

    ####keep only the variables from internal LASSO selection & the outcome
    tokeep<-c("AGE","WEIGHTBL","EDSSBL","RLPS3YR","TRELMOS","PRMSGR","REGION","NHPTMBL","GDLESBL","SFPCSBL","SENSORBL","DISTWKBL","RELAPSE2year")
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
    finalmodel<-lrm(RELAPSE2year~AGE+WEIGHTBL+EDSSBL+RLPS3YR+TRELMOS+PRMSGR+REGION+NHPTMBL+GDLESBL+SFPCSBL+SENSORBL+DISTWKBL+AGE*SFPCSBL+WEIGHTBL*TRELMOS+WEIGHTBL*SFPCSBL+EDSSBL*RLPS3YR+EDSSBL*PRMSGR+EDSSBL*REGION+EDSSBL*GDLESBL+PRMSGR*REGION+REGION*NHPTMBL+REGION*GDLESBL+GDLESBL*SENSORBL+SFPCSBL*DISTWKBL,linear.predictors = T,x=TRUE,y=TRUE,data=X)
    finalmodel
    val<-validate(finalmodel,method="boot",B=500)
    c_index <- abs(val[1,5])/2 + 0.5
    c_slope<-val[4,5]
    ##check the assumption of linearity
    ##best choice of knots and test for linearity
    #a<-lrm(RELAPSE2year~rcs(AGE,5)+rcs(WEIGHTBL,5)+rcs(EDSSBL,5)+rcs(RLPS3YR,5)+rcs(TRELMOS,5)+PRMSGR+REGION+rcs(NHPTDHBL,5)+rcs(GDLESBL,5)+rcs(SFPCSBL,5)+SENSORBL+DISTWKBL,data=X)
    #anova(a)##relapse3years
    #b<-lrm(RELAPSE2year~rcs(AGE,4)+rcs(WEIGHTBL,4)+rcs(EDSSBL,4)+rcs(RLPS3YR,4)+rcs(TRELMOS,4)+PRMSGR+REGION+rcs(NHPTDHBL,4)+rcs(GDLESBL,4)+rcs(SFPCSBL,4)+SENSORBL+DISTWKBL,data=X)
    #anova(b)##relapse3years non linear
    #AIC(a)
    #AIC(b) ###best choice of knots is knots=4


    #################### Make Shrunk coefficients ##############################

    # Make shrunk models
    finalmodel.shrunk <- finalmodel
    finalmodel.shrunk$coef <- val[4,5] * finalmodel.shrunk$coef  # use result from bootstrapping
    shrinkage <- (finalmodel$stats[3]-finalmodel$stats[4])/finalmodel$stats[3]
    finalmodel.shrunk$coef <- shrinkage * finalmodel$coef
    # val.full[4,5] is shrinkage factor; heuristic estimate (LR - df) / LR = 0.76

    # Estimate new intercept, with shrunk lp as offset variable, i.e. coef fixed at unity
    finalmodel.shrunk$coef[1] <- lrm.fit(y=finalmodel$y, offset= finalmodel$x %*% finalmodel.shrunk$coefficients[2:44])$coefficients[1]
    finalmodel.shrunk$coef /finalmodel$coef
    signif(finalmodel.shrunk$coef,2)
    finalmodel.shrunk

    #### bootstrap validation for the shrunk model
    val.shrunk<-validate(finalmodel.shrunk,method="boot",B=500)
    c_index.shrunk <- abs(val.shrunk[1,5])/2 + 0.5
    c_slope.shrunk<-val.shrunk[4,5]

    # Make a penalized model
    penalized	<- pentrace(finalmodel, seq(0,100,0.1))
    finalmodel.pen <- update (finalmodel, penalty=penalized$penalty) #100
    #### bootstrap validation for the penalized model
    val.pen<-validate(finalmodel.pen,method="boot",B=500)
    c_index.pen <- abs(val.pen[1,5])/2 + 0.5
    c_slope.pen<-val.pen[4,5]

    ####other penalized way - not sure about that!!!
    set.seed(1)
    cv.penalty<-cv.glmnet(x=finalmodel$x,y=finalmodel$y,family="binomial", alpha=0,type.measure = "auc")
    best.lambdase=cv.penalty$lambda.1se #71.38
    pen2<-pentrace(finalmodel,c(2*best.lambdase))
    finalmodel.pen2 <- update (finalmodel, penalty=pen2$penalty) #142.7733

    #### bootstrap validation for the penalized model
    val.pen2<-validate(finalmodel.pen2,method="boot",B=500)
    c_index.pen2 <- abs(val.pen2[1,5])/2 + 0.5
    c_slope.pen2<-val.pen2[4,5]


    ####Ridge shrinkage
    set.seed(123)
    cv <- cv.glmnet(x=finalmodel$x, y=finalmodel$y, alpha=0, family=c("binomial"), type.measure = "auc")
    # Display the best lambda value
    cv$lambda.min # 0.747858
    # Fit the model on the training data
    model.L2 <- glmnet(x=finalmodel$x, y=finalmodel$y, alpha = 0, lambda = cv$lambda.min, family=c("binomial"),)

    # Display regression coefficients
    coef(model.L2)
    finalmodel.ridge<-finalmodel
    finalmodel.ridge$coefficients<-as.numeric(coef(model.L2))
    val.ridge<-validate(finalmodel.ridge,method="boot",B=500)
    c_index.ridge <- abs(val.ridge[1,5])/2 + 0.5
    c_slope.ridge<-val.ridge[4,5]

    #####LASSO shrinkage

    set.seed(123)
    cv <- cv.glmnet(x=finalmodel$x, y=finalmodel$y, family=c("binomial"), type.measure = "auc")
    # Display the best lambda value
    cv$lambda.min # 0.012
    # Fit the model on the training data
    model.L1 <- glmnet(x=finalmodel$x, y=finalmodel$y,lambda = cv$lambda.min, family=c("binomial"),)

    # Display regression coefficients
    coef(model.L1)
    finalmodel.lasso<-finalmodel
    finalmodel.lasso$coefficients<-as.numeric(coef(model.L1))
    val.lasso<-validate(finalmodel.lasso,method="boot",B=500)
    c_index.lasso <- abs(val.lasso[1,5])/2 + 0.5
    c_slope.lasso<-val.lasso[4,5]

    #######EPV

    df<-finalmodel[["stats"]][["d.f."]]
    events<- nrow(X[which(X$RELAPSE2year==1),])
    ####return back
    cat("The selected variables are:" , cv.pf.em.both, fill = TRUE)
    EPV<-events/df
    cat("The EPV of the model is", EPV, fill=TRUE)
    if (EPV>10){cat("The EPV of the model is >10, as recommended.", fill=TRUE)}
    if (EPV<10){cat("The EPV of the model is <10, possibility of overfitting problem.", fill = TRUE)}
    #cat("The final model is:",fill=TRUE)
    cat("The bootstap corrected discrimination of the model is", c_index, "and the bootstap corrected c-slope is", c_slope, fill=TRUE)
    cat("The bootstap corrected discrimination after the shrinkage of coefficients (uniform) of the model is", c_index.shrunk, "and the bootstap corrected c-slope is", c_slope.shrunk, fill=TRUE)
    cat("The bootstap corrected discrimination after penalization of coefficients (Penalized ML) of the model is", c_index.pen, "and the bootstap corrected c-slope is", c_slope.pen,". The penalty that was used is", penalized$penalty, fill=TRUE)


    discrimination<-c(c_index,c_index.shrunk,c_index.pen,c_index.pen2,c_index.ridge,c_index.lasso)
    calibration<-c(c_slope,c_slope.shrunk,c_slope.pen,c_slope.pen2,c_slope.ridge,c_slope.lasso)
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
    X.both<-model.matrix(X$RELAPSE2year~.,data=X)
    #### delete NAs values (LASSO requirement)
    X.both<-na.omit(X.both[,-1])
    finalmodel<-lrm(RELAPSE2year~AGE+SEX+RACE+EDSSBL+ONSYRS+RLPS3YR+TRELMOS+PRMSGR+T25FWABL+NHPTMBL+PASATABL+VFT25BL+SFPCSBL+SFMCSBL,x=TRUE,y=TRUE,linear.predictors = TRUE,data=X)
    finalmodel
    set.seed(1)
    val<-validate(finalmodel,method="boot",B=500)
    c_index <- abs(val[1,5])/2 + 0.5
    c_slope<-val[4,5]
    ##check the assumption of linearity
    ##best choice of knots and test for linearity
    #a<-lrm(RELAPSE2year~rcs(AGE,5)+rcs(WEIGHTBL,5)+rcs(EDSSBL,5)+rcs(RLPS3YR,5)+rcs(TRELMOS,5)+PRMSGR+REGION+rcs(NHPTDHBL,5)+rcs(GDLESBL,5)+rcs(SFPCSBL,5)+SENSORBL+DISTWKBL,data=X)
    #anova(a)##relapse3years
    #b<-lrm(RELAPSE2year~rcs(AGE,4)+rcs(WEIGHTBL,4)+rcs(EDSSBL,4)+rcs(RLPS3YR,4)+rcs(TRELMOS,4)+PRMSGR+REGION+rcs(NHPTDHBL,4)+rcs(GDLESBL,4)+rcs(SFPCSBL,4)+SENSORBL+DISTWKBL,data=X)
    #anova(b)##relapse3years non linear
    #AIC(a)
    #AIC(b) ###best choice of knots is knots=4

    #################### Make Shrunk coefficients ##############################

    # Make shrunk models
    finalmodel.shrunk <- finalmodel
    finalmodel.shrunk$coef <- val[4,5] * finalmodel.shrunk$coef  # use result from bootstrapping
    shrinkage <- (finalmodel$stats[3]-finalmodel$stats[4])/finalmodel$stats[3]
    finalmodel.shrunk$coef <- shrinkage * finalmodel$coef
    # val.full[4,5] is shrinkage factor; heuristic estimate (LR - df) / LR = (62.6-8)/62.6=0.88

    # Estimate new intercept, with shrunk lp as offset variable, i.e. coef fixed at unity
    finalmodel.shrunk$coef[1] <- lrm.fit(y=finalmodel$y, offset= finalmodel$x %*% finalmodel.shrunk$coefficients[2:15])$coefficients[1]
    finalmodel.shrunk$coef /finalmodel$coef
    signif(finalmodel.shrunk$coef,2)
    finalmodel.shrunk

    #### bootstrap validation for the shrunk model
    val.shrunk<-validate(finalmodel.shrunk,method="boot",B=500)
    c_index.shrunk <- abs(val.shrunk[1,5])/2 + 0.5
    c_slope.shrunk<-val.shrunk[4,5]
    # Make a penalized model
    penalized	<-  pentrace(finalmodel, seq(0,100,0.1)) #70.3
    finalmodel.pen <- update (finalmodel, penalty=penalized$penalty)
    #### bootstrap validation for the penalized model
    val.pen<-validate(finalmodel.pen,method="boot",B=500)
    c_index.pen <- abs(val.pen[1,5])/2 + 0.5
    c_slope.pen<-val.pen[4,5]
    ####other penalized way - not sure about that!!!
    set.seed(1)
    cv.penalty<-cv.glmnet(x=finalmodel$x,y=finalmodel$y,family="binomial", alpha=0,type.measure = "auc")
    best.lambdase=cv.penalty$lambda.1se #4.76
    pen2<-pentrace(finalmodel,c(2*best.lambdase))
    finalmodel.pen2 <- update (finalmodel, penalty=pen2$penalty) #9.52

    #### bootstrap validation for the penalized model
    val.pen2<-validate(finalmodel.pen2,method="boot",B=500)
    c_index.pen2 <- abs(val.pen2[1,5])/2 + 0.5
    c_slope.pen2<-val.pen2[4,5]


    ####Ridge shrinkage
    set.seed(123)
    cv <- cv.glmnet(x=finalmodel$x, y=finalmodel$y, alpha=0, family=c("binomial"), type.measure = "auc")
    # Display the best lambda value
    cv$lambda.min # 0.1264
    # Fit the model on the training data
    model.L2 <- glmnet(x=finalmodel$x, y=finalmodel$y, alpha = 0, lambda = cv$lambda.min, family=c("binomial"),)

    # Display regression coefficients
    coef(model.L2)
    finalmodel.ridge<-finalmodel
    finalmodel.ridge$coefficients<-as.numeric(coef(model.L2))
    val.ridge<-validate(finalmodel.ridge,method="boot",B=500)
    c_index.ridge <- abs(val.ridge[1,5])/2 + 0.5
    c_slope.ridge<-val.ridge[4,5]

    #####LASSO shrinkage

    set.seed(123)
    cv <- cv.glmnet(x=finalmodel$x, y=finalmodel$y, family=c("binomial"), type.measure = "auc")
    # Display the best lambda value
    cv$lambda.min # 0.005
    # Fit the model on the training data
    model.L1 <- glmnet(x=finalmodel$x, y=finalmodel$y,lambda = cv$lambda.min, family=c("binomial"),)

    # Display regression coefficients
    coef(model.L1)
    finalmodel.lasso<-finalmodel
    finalmodel.lasso$coefficients<-as.numeric(coef(model.L1))
    val.lasso<-validate(finalmodel.lasso,method="boot",B=500)
    c_index.lasso <- abs(val.lasso[1,5])/2 + 0.5
    c_slope.lasso<-val.lasso[4,5]

    ###EPV calculation
    df<-finalmodel[["stats"]][["d.f."]]
    events<- nrow(X[which(X$RELAPSE2year==1),])
    ####return back
    EPV<-events/df
    cat("The EPV of the model is", EPV, fill=TRUE)
    if (EPV>10){cat("The EPV of the model is >10, as recommended.", fill=TRUE)}
    if (EPV<10){cat("The EPV of the model is <10, possibility of overfitting problem.", fill = TRUE)}
    #cat("The final model is:",fill=TRUE)
    cat("The bootstap corrected discrimination of the model is", c_index, "and the bootstap corrected c-slope is", c_slope,fill=TRUE)
    cat("The bootstap corrected discrimination after the shrinkage of coefficients (uniform) of the model is", c_index.shrunk, "and the bootstap corrected c-slope is", c_slope.shrunk, fill=TRUE)
    cat("The bootstap corrected discrimination after penalization of coefficients (Penalized ML) of the model is", c_index.pen, "and the bootstap corrected c-slope is", c_slope.pen,".The penalty that was used is", penalized$penalty, fill=TRUE)
    discrimination<-c(c_index,c_index.shrunk,c_index.pen,c_index.pen2,c_index.ridge,c_index.lasso)
    calibration<-c(c_slope,c_slope.shrunk,c_slope.pen,c_slope.pen2,c_slope.ridge,c_slope.lasso)

    return(list(discrimination=discrimination, calibration=calibration))


  }

}
