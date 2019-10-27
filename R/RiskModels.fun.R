###############################################################################
################      FUNCTION FOR CHECKING       #########################
######## 5 DIFFERENT MODELS AND THEIR DISCRIMINATION ABILITY #######

RiskModels.fun<-function(dataset,model){

  ##loading libraries
  #library(glmnet)
  #library(Hmisc)
  #library(rms)
  #library(glmpath)

  if (model=="LASSOModel") {

    ################# LASSO preperation ##########################################
    ################################################################


    #drop STUDYID, USIBJID (no needed for the matrix), drop TRT01A (blinded to treatment)
    todrop<-c("STUDYID","USUBJID","TRT01A")
    X<-dataset[ , !(names(dataset) %in% todrop)]
    ### develop the matrix for the model
    X.both<-model.matrix(X$RELAPSE2year~.-1,data=X)
    #### delete NAs values (LASSO requirement)
    X.both<-na.omit(X.both)
    X<-na.omit(X)
    set.seed(123)
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

  # Method of shrinkage. Model with LASSO selected variables and coefficients
    ###logistic model based on LASSO selected variables
    finalmodel<-lrm(RELAPSE2year~AGE+WEIGHTBL+EDSSBL+RLPS1YR+PRMSGR+REGION+GDVOLBL+SFPCSBL+DISTWKBL,x=TRUE,y=TRUE,linear.predictors = T,data=X)
    finalmodel
    lassomodel<-finalmodel
    #lasso coefficients
    lassocoef<-as.matrix(cv.coef.both)
    ##insert lasso coefficients to the model
    lassomodel$coefficients[1]<-lassocoef[1]
    lassomodel$coefficients[2]<-lassocoef[2]
    lassomodel$coefficients[3]<-lassocoef[7]
    lassomodel$coefficients[4]<-lassocoef[8]
    lassomodel$coefficients[5]<-lassocoef[11]
    lassomodel$coefficients[6]<-lassocoef[15]
    lassomodel$coefficients[7]<-lassocoef[16]
    lassomodel$coefficients[8]<-lassocoef[17]
    lassomodel$coefficients[9]<-lassocoef[18]
    lassomodel$coefficients[10]<-lassocoef[19]
    lassomodel$coefficients[11]<-lassocoef[23]
    lassomodel$coefficients[12]<-lassocoef[25]
    lassomodel$coefficients[13]<-lassocoef[44]


    ##validate the model (we do not need SE's for the validation)
    set.seed(1)
    val.lassomodel<-validate(lassomodel,method="boot",B=500) #### For validation we do not need SE's
    ##also for the predictions we do not need SE's
    c_index.lassomodel <- abs(val.lassomodel[1,5])/2 + 0.5
    c_slope.lassomodel<-val.lassomodel[4,5]


    ###return
     cat("The selected variables and their coefficients by LASSO model are:" )
     print(lassomodel)
     cat("The bootstap corrected discrimination of the LASSO model is", c_index.lassomodel, "and the bootstap corrected c-slope is", c_slope.lassomodel, fill=TRUE)

    discrimination<-c(c_index.lassomodel)
    calibration<-c(c_slope.lassomodel)
    return(list(discrimination=discrimination, calibration=calibration, lassomodel=lassomodel))

  }


  if (model=="FabioModel") {


    #drop STUDYID, USIBJID (no needed for the matrix), drop TRT01A (blinded to treatment)

    todrop<-c("STUDYID","USUBJID","TRT01A")
    X<-dataset[ , !(names(dataset) %in% todrop)]
    X<-na.omit(X)

    finalmodel<-lrm(RELAPSE2year~AGE+SEX+EDSSBL+ONSYRS+RACE+RLPS1YR+TRELMOS+PRMSGR+T25FWABL+NHPTMBL+PASATABL+VFT25BL+SFPCSBL+SFMCSBL,x=TRUE,y=TRUE,linear.predictors = TRUE,data=X)
    finalmodel


### Method of shrinkage. Make a PMLE penalized model based on minimizing a modified AIC criterion
    # Make a penalized model
    set.seed(1)
    penalized	<-  pentrace(finalmodel, seq(0,200,0.1))
    finalmodel.pen <- update (finalmodel, penalty=penalized$penalty)
    #### bootstrap validation for the penalized model
    set.seed(1)
    val.pen<-validate(finalmodel.pen,method="boot",B=500)
    c_index.pen <- abs(val.pen[1,5])/2 + 0.5
    c_slope.pen<-val.pen[4,5]


    cat("Pellegrini's et al. model is:")
    print(finalmodel.pen)
    cat("The bootstap corrected discrimination of the model after PMLE of coefficients is", c_index.pen, "and the bootstap corrected c-slope is", c_slope.pen, fill=TRUE)

    discrimination<-c(c_index.pen)
    calibration<-c(c_slope.pen)
    return(list(discrimination=discrimination, calibration=calibration, fabiomodel=finalmodel.pen))


  }

}




