####################################################################
################## Function to add columns arm and meanRisk ############
################## that are needed for the IPD NMR model ################

DataForIPDNMR.fun<-function(){
RiskData$meanlogitRisk<-NA
RiskData$meanlogitRisk[RiskData$STUDYID==1]<-mean(RiskData$logitRisk[RiskData$STUDYID==1])
RiskData$meanlogitRisk[RiskData$STUDYID==2]<-mean(RiskData$logitRisk[RiskData$STUDYID==2])
RiskData$meanlogitRisk[RiskData$STUDYID==3]<-mean(RiskData$logitRisk[RiskData$STUDYID==3])

RiskData$meanRisk<-NA
RiskData$meanRisk[RiskData$STUDYID==1]<-mean(RiskData$Risk[RiskData$STUDYID==1])
RiskData$meanRisk[RiskData$STUDYID==2]<-mean(RiskData$Risk[RiskData$STUDYID==2])
RiskData$meanRisk[RiskData$STUDYID==3]<-mean(RiskData$Risk[RiskData$STUDYID==3])

RiskData$arm<-NA
RiskData$arm[RiskData$STUDYID==1 & RiskData$TRT01A==1]<-1
RiskData$arm[RiskData$STUDYID==1 & RiskData$TRT01A==4]<-2
RiskData$arm[RiskData$STUDYID==2 & RiskData$TRT01A==1]<-1
RiskData$arm[RiskData$STUDYID==2 & RiskData$TRT01A==2]<-2
RiskData$arm[RiskData$STUDYID==2 & RiskData$TRT01A==4]<-3
RiskData$arm[RiskData$STUDYID==3 & RiskData$TRT01A==3]<-1
RiskData$arm[RiskData$STUDYID==3 & RiskData$TRT01A==4]<-2

return(RiskData)
}