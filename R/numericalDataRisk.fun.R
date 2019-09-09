###########################################################################################################################
##############  function that selects variables of interest and ####################
##############         recodes some of factors        #############################
#####################################################################################################################################################

numericalDataRisk.fun=function(dataset){
  library(car)
  ##################################HANDLING VARIABLES##############################################
  ##keep only needed variables
  keep<-c("STUDYID","USUBJID","AGE","SEX","RACE","TRT01A","HEIGHTBL","WEIGHTBL","EDSSBL","ONSYRS","DIAGYRS",
          "DOMIHAND","RLPS3YR","TRELMOS", "MCDBL","PRMSGR","REGION","T25FWABL","NHPTMBL","NHPTDHBL","NHPTNHBL",
          "PASATABL","MSFCBL","GDLESBL","T2VOLBL","T1VOLBL","BVZBL","VFT100BL","VFT25BL",
          "VFT125BL","SFPCSBL","SFMCSBL","VISUALBL","BRAINBL","PYRAMIBL","SENSORBL","BOWLBLBL",
          "CEREBRBL","DISTWKBL","T25FWP1","NHPTMP1","NHPTDHP1","NHPTNHP1","PASATP1","T25FWPC","NHPTMPC","NHPTDHPC",
          "NHPTNHPC","PASATPC","RELAPSE1year","RELAPSE2year")
  MSrelapse<-dataset[,keep]
  
  ##################################HANDLING VARIABLES############################
  ###########################RECODE VARIABLES and make them factors############################
  MSrelapse$SEX<-recode(MSrelapse$SEX, "'M'=1; 'F'=0")
  MSrelapse$RACE<-recode(MSrelapse$RACE, "'WHITE'=1; 'NON-WHITE'=0")
  MSrelapse$DOMIHAND<-recode(MSrelapse$DOMIHAND, "'Left'=1;'LEFT'=1;  'Right'=0; 'RIGHT'=0;" )
  MSrelapse$DOMIHAND[which(MSrelapse$DOMIHAND=="")]<-NA
  ###MCDBL categories instead of 1, 2, 3, 4 there are  1, 2, >=3
  MSrelapse$MCDBL[which(MSrelapse$MCDBL==4)]<-3
  MSrelapse$MCDBL<-as.factor(MSrelapse$MCDBL)
  MSrelapse$PRMSGR<-as.factor(MSrelapse$PRMSGR)
  MSrelapse$REGION<-recode(MSrelapse$REGION, "'Eastern Europe'=1;'India'=2;  'North America'=3; 'ROW'=4;'Western Europe'=5 " )
  ###VISUALBL categories instead of 0, 1, 2, 3, 4, 5, 6 there are 0, 1, 2, >=3
  MSrelapse$VISUALBL[which(MSrelapse$VISUALBL==4 | MSrelapse$VISUALBL==5 | MSrelapse$VISUALBL==6)]<-3
  MSrelapse$VISUALBL<-as.factor(MSrelapse$VISUALBL)
  ###BRAINBL categories instead of 0, 1, 2, 3, 4 there are 0, 1, >=2
  MSrelapse$BRAINBL[which(MSrelapse$BRAINBL==3 | MSrelapse$BRAINBL==4)]<-2
  MSrelapse$BRAINBL<-as.factor(MSrelapse$BRAINBL) 
  ###PYRAMIBL categories instead of 0, 1, 2, 3, 4, 5, 6 there are 0, 1, 2, >=3
  MSrelapse$PYRAMIBL[which(MSrelapse$PYRAMIBL==4 | MSrelapse$PYRAMIBL==5)]<-3
  MSrelapse$PYRAMIBL<-as.factor(MSrelapse$PYRAMIBL)

  ###SENSOR BL categories instead of 0, 1, 2, 3, 4, 5, 6 there are 0, 1, 2, >=3
  MSrelapse$SENSORBL[which(MSrelapse$SENSORBL==4 | MSrelapse$SENSORBL==5 | MSrelapse$SENSORBL==6)]<-3
  MSrelapse$SENSORBL<-as.factor(MSrelapse$SENSORBL)
  ###BOWLBLBLL categories instead of 0, 1, 2, 3, 4, 5, 6 there are 0, 1, 2, >=3
   MSrelapse$BOWLBLBL[which(MSrelapse$BOWLBLBL==4 | MSrelapse$BOWLBLBL==5 | MSrelapse$BOWLBLBL==6)]<-3
  MSrelapse$BOWLBLBL<-as.factor(MSrelapse$BOWLBLBL)
  ###CEREBRBL categories instead of 0 ,1, 2, 3, 4 there are 0, 1, >=2
   MSrelapse$CEREBRBL[which(MSrelapse$CEREBRBL==3 | MSrelapse$CEREBRBL==4)]<-2
  MSrelapse$CEREBRBL<-as.factor(MSrelapse$CEREBRBL)
  ####GDLESBL
  MSrelapse$GDLESBL[which(MSrelapse$GDLESBL>=1)]<-1
  MSrelapse$GDLESBL<-as.factor(MSrelapse$GDLESBL)
  ###DISTWKBL
  MSrelapse$DISTWKBL<-as.numeric(MSrelapse$DISTWKBL)
  MSrelapse$DISTWKBL[which(MSrelapse$DISTWKBL==1)]<-NA
  MSrelapse$DISTWKBL[which(MSrelapse$DISTWKBL==2 | MSrelapse$DISTWKBL==3 | MSrelapse$DISTWKBL==4 | MSrelapse$DISTWKBL==5 )]<-0
  MSrelapse$DISTWKBL[which(MSrelapse$DISTWKBL==6)]<-1
  MSrelapse$DISTWKBL<-as.factor(MSrelapse$DISTWKBL)
  ###### transformations of continues variables
  #T25FWABL
  MSrelapse$T25FWABL<-log(MSrelapse$T25FWABL+1)
  #NHPTMBL
  MSrelapse$NHPTMBL<-log(MSrelapse$NHPTMBL)
  #NHPTDHBL
  MSrelapse$NHPTDHBL<-log(MSrelapse$NHPTDHBL+1)
  #NHPTDHBL
  MSrelapse$NHPTNHBL<-log(MSrelapse$NHPTNHBL+1)
  #T2VOLBL
  MSrelapse$T2VOLBL<-log(MSrelapse$T2VOLBL+1)
  #T1VOLBL
  MSrelapse$T1VOLBL<-log(MSrelapse$T1VOLBL+1)
  #T25FWP1
  MSrelapse$T25FWP1<-log(MSrelapse$T25FWP1+1)
  #NHPTMP1
  MSrelapse$NHPTMP1<-log(MSrelapse$NHPTMP1+1)
  #NHPTDHP1
  MSrelapse$NHPTDHP1<-log(MSrelapse$NHPTDHP1+1)
  #NHPTNHP1
  MSrelapse$NHPTNHP1<-log(MSrelapse$NHPTNHP1+1)
  #T25FWPC
  MSrelapse$T25FWPC<-log(MSrelapse$T25FWPC+100)
  #NHPTMPC
  MSrelapse$NHPTMPC<-log(MSrelapse$NHPTMPC+100)
  #NHPTDHPC
  MSrelapse$NHPTDHPC<-log(MSrelapse$NHPTDHPC+100)
  #NHPTNHPC
  MSrelapse$NHPTNHPC<-log(MSrelapse$NHPTNHPC+100)
  #PASATPC
  MSrelapse$PASATPC<-log(MSrelapse$PASATPC+101)
  
  ##### remove highly correlated variables (after checking with Spearman correlation)
  todrop<- c("PASATP1", "NHPTDHBL", "NHPTNHBL", "NHPTMP1", "NHPTDHP1", "NHPTNHP1", "T1VOLBL", "VFT100BL","VFT125BL","T25FWP1", "NHPTDHPC", "NHPTNHPC","DIAGYRS")
  MSrelapse <- MSrelapse[,!(names(MSrelapse) %in% todrop)]
  #remove Sentinel study - Not included in AD data - Combination of therapies 
  MSrelapse<-MSrelapse[which(MSrelapse$STUDYID!="SENTINEL"),]
  return(MSrelapse)
}
