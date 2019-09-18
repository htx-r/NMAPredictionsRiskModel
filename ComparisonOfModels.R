########################## Script for comparing the prognostic models #####################

###creation of a 12x1 table with all model's discriminations
InternalModel$discrimination<-as.data.frame(InternalModel$discrimination, row.names=(c("Model 1 without shrinkage", "Model 1 uniform", "Model 1 penalized", "Model 1 penalized 2", "Model 1 ridge", "Model 1 lasso")))
colnames(InternalModel$discrimination)<-c("discrimination")
InternalSplinesSignModel$discrimination<-as.data.frame(InternalSplinesSignModel$discrimination, row.names=(c("Model 2 without shrinkage", "Model 2 uniform", "Model 2 penalized","Model 2 penalized 2", "Model 2 ridge", "Model 2 lasso")))
colnames(InternalSplinesSignModel$discrimination)<-c("discrimination")
InternalInteractionsModel$discrimination<-as.data.frame(InternalInteractionsModel$discrimination,row.names=(c("Model 3 without shrinkage", "Model 3 uniform", "Model 3 penalized","Model 3 penalized 2", "Model 3 ridge", "Model 3 lasso")))
colnames(InternalInteractionsModel$discrimination)<-c("discrimination")
FabioModel$discrimination<-as.data.frame(FabioModel$discrimination, row.names=(c("Model 4 without shrinkage", "Model 4 uniform", "Model 4 penalized", "Model 4 penalized 2", "Model 4 ridge", "Model 4 lasso")))
colnames(FabioModel$discrimination)<-c("discrimination")

comparisonofmodelstabledis<-rbind(InternalModel$discrimination,InternalSplinesSignModel$discrimination,InternalInteractionsModel$discrimination,FabioModel$discrimination)

###creation of a 12x1 table with all model's calibrations
InternalModel$calibration<-as.data.frame(InternalModel$calibration, row.names=(c("Model 1 without shrinkage", "Model 1 uniform", "Model 1 penalized", "Model 1 penalized 2", "Model 1 ridge", "Model 1 lasso")))
colnames(InternalModel$calibration)<-c("calibration")
InternalSplinesSignModel$calibration<-as.data.frame(InternalSplinesSignModel$calibration, row.names=(c("Model 2 without shrinkage", "Model 2 uniform", "Model 2 penalized", "Model 2 penalized 2", "Model 2 ridge", "Model 2 lasso")))
colnames(InternalSplinesSignModel$calibration)<-c("calibration")
InternalInteractionsModel$calibration<-as.data.frame(InternalInteractionsModel$calibration,row.names=(c("Model 3 without shrinkage", "Model 3 uniform", "Model 3 penalized", "Model 3 penalized 2", "Model 3 ridge", "Model 3 lasso")))
colnames(InternalInteractionsModel$calibration)<-c("calibration")
FabioModel$calibration<-as.data.frame(FabioModel$calibration, row.names=(c("Model 4 without shrinkage", "Model 4 uniform", "Model 4 penalized", "Model 4 penalized 2", "Model 4 ridge", "Model 4 lasso")))
colnames(FabioModel$calibration)<-c("calibration")

comparisonofmodelstablecal<-rbind(InternalModel$calibration,InternalSplinesSignModel$calibration,InternalInteractionsModel$calibration,FabioModel$calibration)

### 12x2 table with discrimination and calibration for all the models and all shrinkage approaches
ComparisonOfModelsTable<-cbind(comparisonofmodelstabledis,comparisonofmodelstablecal)

rm(comparisonofmodelstablecal)
rm(comparisonofmodelstabledis)
