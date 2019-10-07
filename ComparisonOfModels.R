########################## Script for comparing the prognostic models #####################

###creation of a 12x1 table with all model's discriminations
InternalModel$discrimination<-as.data.frame(InternalModel$discrimination, row.names=(c("Model 1 with LASSO coefficients and shirnkage", "Model 1 LASSO 2", "Model 1 uniform shrinkage", "Model 1 ridge AIC", "Model 1 ridge AUC")))
colnames(InternalModel$discrimination)<-c("discrimination")
InternalSplinesSignModel$discrimination<-as.data.frame(InternalSplinesSignModel$discrimination, row.names=(c("Model 2 LASSO 2", "Model 2 uniform shrinkage", "Model 2 ridge AIC","Model 2 ridge AUC")))
colnames(InternalSplinesSignModel$discrimination)<-c("discrimination")
InternalInteractionsModel$discrimination<-as.data.frame(InternalInteractionsModel$discrimination,row.names=(c("Model 3 LASSO 2", "Model 3 uniform shrinkage", "Model 3 ridge AIC","Model 3 ridge AUC")))
colnames(InternalInteractionsModel$discrimination)<-c("discrimination")
FabioModel$discrimination<-as.data.frame(FabioModel$discrimination, row.names=(c("Model 4 LASSO 2", "Model 4 uniform shrinkage", "Model 4 ridge AIC","Model 4 ridge AUC")))
colnames(FabioModel$discrimination)<-c("discrimination")

comparisonofmodelstabledis<-rbind(InternalModel$discrimination,InternalSplinesSignModel$discrimination,InternalInteractionsModel$discrimination,FabioModel$discrimination)

###creation of a 12x1 table with all model's calibrations
InternalModel$calibration<-as.data.frame(InternalModel$calibration, row.names=(c("Model 1 with LASSO coefficients and shirnkage", "Model 1 LASSO 2", "Model 1 uniform shrinkage", "Model 1 ridge AIC", "Model 1 ridge AUC")))
colnames(InternalModel$calibration)<-c("calibration")
InternalSplinesSignModel$calibration<-as.data.frame(InternalSplinesSignModel$calibration, row.names=(c("Model 2 LASSO 2", "Model 2 uniform shrinkage", "Model 2 ridge AIC","Model 2 ridge AUC")))
colnames(InternalSplinesSignModel$calibration)<-c("calibration")
InternalInteractionsModel$calibration<-as.data.frame(InternalInteractionsModel$calibration,row.names=(c("Model 3 LASSO 2", "Model 3 uniform shrinkage", "Model 3 ridge AIC","Model 3 ridge AUC")))
colnames(InternalInteractionsModel$calibration)<-c("calibration")
FabioModel$calibration<-as.data.frame(FabioModel$calibration, row.names=(c("Model 4 LASSO 2", "Model 4 uniform shrinkage", "Model 4 ridge AIC","Model 4 ridge AUC")))
colnames(FabioModel$calibration)<-c("calibration")

comparisonofmodelstablecal<-rbind(InternalModel$calibration,InternalSplinesSignModel$calibration,InternalInteractionsModel$calibration,FabioModel$calibration)

### 12x2 table with discrimination and calibration for all the models and all shrinkage approaches
ComparisonOfModelsTable<-cbind(comparisonofmodelstabledis,comparisonofmodelstablecal)

rm(comparisonofmodelstablecal)
rm(comparisonofmodelstabledis)
