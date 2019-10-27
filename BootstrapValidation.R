ApparentLASSO_discrimination<-LASSOModel[["lassomodel"]][["stats"]][6]
ApparentFabio_discrimination<-FabioModel[["fabiomodel"]][["stats"]][6]

data1<-na.omit(MSrelapse)
todrop<-c("STUDYID","USUBJID","TRT01A")
data1<-data1[ , !(names(data1) %in% todrop)]

manual_boot <- function(data,samples){

  resultslasso <- matrix(nrow = samples,ncol = 6)
  resultsFabio <- matrix(nrow = samples,ncol = 6)
  set.seed(231398)
  for (i in 1:samples) {
    samp_index <- sample(1:nrow(data), nrow(data), rep=TRUE) # create a sampling index vector

    bs_samp <- data[samp_index,] # index the orignal dataset using the sampling vector to give the bs sample
    bs_samp_modelMatrix<-model.matrix(bs_samp$RELAPSE2year~.-1,data=bs_samp)
################################ Model LASSO ##############################
    cv.fit.both<-cv.glmnet(x=bs_samp_modelMatrix,y=bs_samp$RELAPSE2year,family="binomial",type.measure = "auc")
    cv.coef.both<-coef(cv.fit.both,s="lambda.1se")
    lassocoef<-as.matrix(cv.coef.both)[-3]
    finalmodel<-lrm(RELAPSE2year~.,x=TRUE,y=TRUE,linear.predictors = T,data=bs_samp)
    lassomodel<-finalmodel
        for (j in 1:(finalmodel[["stats"]][["d.f."]]+1)){
          lassomodel$coefficients[j]<-lassocoef[j]
        }

    modellasso<-lassomodel
    lp_bslasso <- predict(modellasso) # predict lp from the bootstrap model in the bs sample
    pr_bslasso <- exp(lp_bslasso)/(1+exp(lp_bslasso))# predict probabilities from the bootstrap model in the bs sample

    lp_testlasso <- predict(modellasso, newdata = data) # predict lp from the bootstrap model in the original sample
    pr_testlasso <-exp(lp_testlasso)/(1+exp(lp_testlasso))# predict probabilities from the bootstrap model in the original sample

    # calculate the apparent performance of the bootstrap model in the bs sample
    app_cstat_modellasso <- roc(RELAPSE2year~pr_bslasso,data=bs_samp)
    resultslasso[i,1] <- as.numeric(app_cstat_modellasso$auc)
    app_citl_modellasso <- glm(RELAPSE2year ~ offset(lp_bslasso),family=binomial, data=bs_samp)
    resultslasso[i,2] <- summary(app_citl_modellasso)$coefficients[1,1]
    app_cslope_modellasso <- glm(RELAPSE2year ~ lp_bslasso,family=binomial(link='logit'), data=bs_samp)
    resultslasso[i,3] <- summary(app_cslope_modellasso)$coefficients[2,1]

    # calculate the test performance of the bootstrap model in the original sample
    test_cstat_modellasso <- roc(RELAPSE2year~pr_testlasso,data=data)
    resultslasso[i,4] <- test_cstat_modellasso$auc
    test_citl_modellasso <- glm(RELAPSE2year ~ offset(lp_testlasso),family=binomial, data=data)
    resultslasso[i,5] <- summary(test_citl_modellasso)$coefficients[1,1]
    test_cslope_modellasso <- glm(RELAPSE2year ~ lp_testlasso,family=binomial, data=data)
    resultslasso[i,6] <- summary(test_cslope_modellasso)$coefficients[2,1]

##################### FABIO'S Model ####################################################

    modelFabio1<-lrm(RELAPSE2year~AGE+SEX+EDSSBL+ONSYRS+RACE+RLPS1YR+TRELMOS+PRMSGR+T25FWABL+NHPTMBL+PASATABL+VFT25BL+SFPCSBL+SFMCSBL,x=TRUE,y=TRUE,linear.predictors = TRUE,data=bs_samp)
    penalized	<-  pentrace(modelFabio1, seq(0,200,0.1))
    modelFabio <- update (modelFabio1, penalty=penalized$penalty)

    lp_bsFabio<- predict(modelFabio) # predict lp from the bootstrap model in the bs sample
    pr_bsFabio <- exp(lp_bsFabio)/(1+exp(lp_bsFabio))# predict probabilities from the bootstrap model in the bs sample

    lp_testFabio <- predict(modelFabio, newdata = data) # predict lp from the bootstrap model in the original sample
    pr_testFabio<-exp(lp_testFabio)/(1+exp(lp_testFabio))# predict probabilities from the bootstrap model in the original sample

    # calculate the apparent performance of the bootstrap model in the bs sample
    app_cstat_modelFabio <- roc(RELAPSE2year~pr_bsFabio,data=bs_samp)
    resultsFabio[i,1] <- as.numeric(app_cstat_modelFabio$auc)
    app_citl_modelFabio <- glm(RELAPSE2year ~ offset(lp_bsFabio),family=binomial, data=bs_samp)
    resultsFabio[i,2] <- summary(app_citl_modelFabio)$coefficients[1,1]
    app_cslope_modelFabio <- glm(RELAPSE2year ~ lp_bsFabio,family=binomial(link='logit'), data=bs_samp)
    resultsFabio[i,3] <- summary(app_cslope_modelFabio)$coefficients[2,1]

    # calculate the test performance of the bootstrap model in the original sample
    test_cstat_modelFabio <- roc(RELAPSE2year~pr_testFabio,data=data)
    resultsFabio[i,4] <- test_cstat_modelFabio$auc
    test_citl_modelFabio <- glm(RELAPSE2year ~ offset(lp_testFabio),family=binomial, data=data)
    resultsFabio[i,5] <- summary(test_citl_modelFabio)$coefficients[1,1]
    test_cslope_modelFabio <- glm(RELAPSE2year ~ lp_testFabio,family=binomial, data=data)
    resultsFabio[i,6] <- summary(test_cslope_modelFabio)$coefficients[2,1]



########################################################################################


  }
  results2lasso <- as.data.frame(resultslasso)
  results2Fabio <- as.data.frame(resultsFabio)
  colnames(results2lasso) <- c("app_c_stat","app_citl","app_c_slope","test_c_stat","test_citl","test_c_slope")
  colnames(results2Fabio) <- c("app_c_stat","app_citl","app_c_slope","test_c_stat","test_citl","test_c_slope")

  return(list(results2lasso,results2Fabio))
}

boot_results <- manual_boot(data=data1,samples=500)
boot_results_LASSO<-boot_results[[1]]
boot_results_Fabio<-boot_results[[2]]

#####for LASSO
# optimism adjusted statistics
C_index_correctedLASSO<-ApparentLASSO_discrimination - (mean(boot_results[[1]]$app_c_stat)-mean(boot_results[[1]]$test_c_stat)) # c-stat
Calibration_slope_correctedLASSO<-1- (mean(boot_results[[1]]$app_c_slope)-mean(boot_results[[1]]$test_c_slope)) # c-slope

#optimism on its own
C_index_optimism_LASSO<-(mean(boot_results[[1]]$app_c_stat)-mean(boot_results[[1]]$test_c_stat))
Calibration_slope_optimismLASSO<-(mean(boot_results[[1]]$app_c_slope)-mean(boot_results[[1]]$test_c_slope)) # c-slope
### For Fabio's
# optimism adjusted statistics
C_index_correctedFabio<-ApparentLASSO_discrimination - (mean(boot_results[[2]]$app_c_stat)-mean(boot_results[[2]]$test_c_stat)) # c-stat
Calibration_slope_correctedFabio<-1- (mean(boot_results[[2]]$app_c_slope)-mean(boot_results[[2]]$test_c_slope)) # c-slope

#optimism on its own
C_index_optimism_Fabio<-(mean(boot_results[[2]]$app_c_stat)-mean(boot_results[[2]]$test_c_stat))
Calibration_slope_optimismFabio<-(mean(boot_results[[2]]$app_c_slope)-mean(boot_results[[2]]$test_c_slope)) # c-slope



###creation of a 2x1 table with both model's discriminations
LASSOModel$discrimination<-as.data.frame(C_index_correctedLASSO, row.names=c("LASSO model"))
colnames(LASSOModel$discrimination)<-c("discrimination")
FabioModel$discrimination<-as.data.frame(C_index_correctedFabio, row.names=c("Fabio model with PMLE shrinkage method"))
colnames(FabioModel$discrimination)<-c("discrimination")

PerformanceTabledis<-rbind(LASSOModel$discrimination,FabioModel$discrimination)

###creation of a 2x1 table with both model's calibrations
LASSOModel$calibration<-as.data.frame(Calibration_slope_correctedLASSO, row.names=c("LASSO model"))
colnames(LASSOModel$calibration)<-c("calibration")
FabioModel$calibration<-as.data.frame(Calibration_slope_correctedFabio, row.names=c("Fabio model with PMLE shrinkage method"))
colnames(FabioModel$calibration)<-c("calibration")

PerformanceTablecal<-rbind(LASSOModel$calibration,FabioModel$calibration)

### 2x2 table with discrimination and calibration for all the models
PerformanceTable<-cbind(PerformanceTabledis,PerformanceTablecal)
