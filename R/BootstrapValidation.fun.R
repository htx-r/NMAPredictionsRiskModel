########################## Function for internal bootstrap validation - corrected for optimism ##################
############################################################################################################
BootstrapValidation.fun <- function(data, samples, alphaElasticNet,modelElasticNet, modelSpecific ){
  resultselasticnet <- matrix(nrow = samples,ncol = 6)
  results2 <- matrix(nrow = samples,ncol = 6)
  set.seed(231398)
  for (i in 1:samples) {
    samp_index <- sample(1:nrow(data), nrow(data), rep=TRUE) # create a sampling index vector

    bs_samp <- data[samp_index,] # index the orignal dataset using the sampling vector to give the bs sample
    bs_samp_modelMatrix<-model.matrix(bs_samp$RELAPSE2year~.-1,data=bs_samp)
    bs_samp_modelMatrix<-bs_samp_modelMatrix[,-2]
################################ Model LASSO ##############################
    cv.fit.both<-cv.glmnet(x=bs_samp_modelMatrix, y=bs_samp$RELAPSE2year, alpha=alphaElasticNet, family="binomial",type.measure = "auc")
    cv.coef.both<-coef(cv.fit.both,s="lambda.1se")
    elasticnetcoef<-as.matrix(cv.coef.both)
    finalmodel<-lrm(RELAPSE2year~.,x=TRUE,y=TRUE,linear.predictors = T,data=bs_samp)
    elasticnetmodel<-finalmodel
        for (j in 1:(finalmodel[["stats"]][["d.f."]]+1)){
          elasticnetmodel$coefficients[j]<-elasticnetcoef[j]
        }

    modelelasticnet<-elasticnetmodel
    lp_bselasticnet <- predict(modelelasticnet,newdata=bs_samp) # predict lp from the bootstrap model in the bs sample
    pr_bselasticnet <- exp(lp_bselasticnet)/(1+exp(lp_bselasticnet))# predict probabilities from the bootstrap model in the bs sample

    lp_testelasticnet<- predict(modelelasticnet, newdata = data) # predict lp from the bootstrap model in the original sample
    pr_testelasticnet<-exp(lp_testelasticnet)/(1+exp(lp_testelasticnet))# predict probabilities from the bootstrap model in the original sample

    # calculate the apparent performance of the bootstrap model in the bs sample
    app_cstat_modelelasticnet <- roc(RELAPSE2year~pr_bselasticnet,data=bs_samp)
    resultselasticnet[i,1] <- as.numeric(app_cstat_modelelasticnet$auc)
    app_citl_modelelasticnet <- glm(RELAPSE2year ~ offset(lp_bselasticnet),family=binomial, data=bs_samp)
    resultselasticnet[i,2] <- summary(app_citl_modelelasticnet)$coefficients[1,1]
    app_cslope_modelelasticnet <- glm(RELAPSE2year ~ lp_bselasticnet,family=binomial(link='logit'), data=bs_samp)
    resultselasticnet[i,3] <- summary(app_cslope_modelelasticnet)$coefficients[2,1]

    # calculate the test performance of the bootstrap model in the original sample
    test_cstat_modelelasticnet <- roc(RELAPSE2year~pr_testelasticnet,data=data)
    resultselasticnet[i,4] <- test_cstat_modelelasticnet$auc
    test_citl_modelelasticnet <- glm(RELAPSE2year ~ offset(lp_testelasticnet),family=binomial, data=data)
    resultselasticnet[i,5] <- summary(test_citl_modelelasticnet)$coefficients[1,1]
    test_cslope_modelelasticnet<- glm(RELAPSE2year ~ lp_testelasticnet,family=binomial, data=data)
    resultselasticnet[i,6] <- summary(test_cslope_modelelasticnet)$coefficients[2,1]

##################### FABIO'S Model ####################################################

    model2<-lrm(RELAPSE2year~AGE+SEX+EDSSBL+ONSYRS+RACE+RLPS1YR+TRELMOS+PRMSGR+T25FWABL+NHPTMBL+PASATABL+VFT25BL+SFPCSBL+SFMCSBL,x=TRUE,y=TRUE,linear.predictors = TRUE,data=bs_samp)
    penalized	<-  pentrace(model2, seq(0,200,0.1))
    model2 <- update (model2, penalty=penalized$penalty)

    lp_bs2<- predict(model2) # predict lp from the bootstrap model in the bs sample
    pr_bs2 <- exp(lp_bs2)/(1+exp(lp_bs2))# predict probabilities from the bootstrap model in the bs sample

    lp_test2 <- predict(model2, newdata = data) # predict lp from the bootstrap model in the original sample
    pr_test2<-exp(lp_test2)/(1+exp(lp_test2))# predict probabilities from the bootstrap model in the original sample

    # calculate the apparent performance of the bootstrap model in the bs sample
    app_cstat_model2<- roc(RELAPSE2year~pr_bs2,data=bs_samp)
    results2[i,1] <- as.numeric(app_cstat_model2$auc)
    app_citl_model2<- glm(RELAPSE2year ~ offset(lp_bs2),family=binomial, data=bs_samp)
    results2[i,2] <- summary(app_citl_model2)$coefficients[1,1]
    app_cslope_model2 <- glm(RELAPSE2year ~ lp_bs2,family=binomial(link='logit'), data=bs_samp)
    results2[i,3] <- summary(app_cslope_model2)$coefficients[2,1]


    # calculate the test performance of the bootstrap model in the original sample
    test_cstat_model2 <- roc(RELAPSE2year~pr_test2,data=data)
    results2[i,4] <- test_cstat_model2$auc
    test_citl_model2 <- glm(RELAPSE2year ~ offset(lp_test2),family=binomial, data=data)
    results2[i,5] <- summary(test_citl_model2)$coefficients[1,1]
    test_cslope_model2 <- glm(RELAPSE2year ~ lp_test2,family=binomial, data=data)
    results2[i,6] <- summary(test_cslope_model2)$coefficients[2,1]



########################################################################################


  }
  results2elasticnet <- as.data.frame(resultselasticnet)
  results2model2 <- as.data.frame(results2)
  colnames(results2elasticnet) <- c("app_c_stat","app_citl","app_c_slope","test_c_stat","test_citl","test_c_slope")
  colnames(results2model2) <- c("app_c_stat","app_citl","app_c_slope","test_c_stat","test_citl","test_c_slope")

  #####for LASSO
  # optimism adjusted statistics
  lpL<- predict(modelElasticNet, newdata=data) # predict lp from the bootstrap model in the bs sample
  prL <- exp(lpL)/(1+exp(lpL))# predict probabilities from the bootstrap model in the bs sample


  appL<-roc(RELAPSE2year~prL, data=data1)
  apparent_dis_elnet<-as.numeric(appL$auc)
  C_index_correctedelnet<-apparent_dis_elnet- (mean(results2elasticnet$app_c_stat)-mean(results2elasticnet$test_c_stat)) # c-stat
  app_cslopeL <- glm(RELAPSE2year ~ lpL,family=binomial(link='logit'), data=data1)
  apparent_cal_elnet <- summary(app_cslopeL)$coefficients[2,1]
  Calibration_slope_correctedelnet<-apparent_cal_elnet - (mean(results2elasticnet$app_c_slope)-mean(results2elasticnet$test_c_slope)) # c-slope


  ### For Fabio's
  # optimism adjusted statistics
  lpF<- predict(modelSpecific, newdata=data1 ) # predict lp from the bootstrap model in the bs sample
  prF <- exp(lpF)/(1+exp(lpF))# predict probabilities from the bootstrap model in the bs sample


  appF<-roc(RELAPSE2year~prF, data=data1)
  apparent_dis_modelSpecific<-as.numeric(appF$auc)
  C_index_correctedmodelSpecific<-apparent_dis_modelSpecific - (mean(results2model2$app_c_stat)-mean(results2model2$test_c_stat)) # c-stat
  app_cslopeF <- glm(RELAPSE2year ~ lpF,family=binomial(link='logit'), data=data1)
  apparent_cal_modelSpecific<- summary(app_cslopeF)$coefficients[2,1]
  Calibration_slope_correctedmodelSpecific<-apparent_cal_modelSpecific- (mean(results2model2$app_c_slope)-mean(results2model2$test_c_slope)) # c-slope

  ###creation of a 2x1 table with both model's discriminations
  ElNetModel<-NA
  ElNetModel$discrimination<-as.data.frame(C_index_correctedelnet, row.names=c("LASSO model"))
  colnames(ElNetModel$discrimination)<-c("discrimination")
  ModelSpecific<-NA
  ModelSpecific$discrimination<-as.data.frame(C_index_correctedmodelSpecific, row.names=c("Pellggrini model with PMLE shrinkage method"))
  colnames(ModelSpecific$discrimination)<-c("discrimination")

  PerformanceTabledis<-rbind(ElNetModel$discrimination,ModelSpecific$discrimination)

  ###creation of a 2x1 table with both model's calibrations
  ElNetModel$calibration<-as.data.frame(Calibration_slope_correctedelnet, row.names=c("LASSO model"))
  colnames(ElNetModel$calibration)<-c("calibration")
  ModelSpecific$calibration<-as.data.frame(Calibration_slope_correctedmodelSpecific, row.names=c("Fabio model with PMLE shrinkage method"))
  colnames(ModelSpecific$calibration)<-c("calibration")

  PerformanceTablecal<-rbind(ElNetModel$calibration,ModelSpecific$calibration)

  ### 2x2 table with discrimination and calibration for all the models
  PerformanceTable<-cbind(PerformanceTabledis,PerformanceTablecal)


  return(list(results2elasticnet,results2model2,  apparent_dis_elnet, apparent_dis_modelSpecific,apparent_cal_elnet, apparent_cal_modelSpecific,  PerformanceTable ))
}

