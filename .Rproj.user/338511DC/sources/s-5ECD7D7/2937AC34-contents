###################################################################################################################
####################  Function that transforms and selects the needed variables #####################################
###################################################################################################################



FinalDataSMSC.fun=function(datapath){

    ## Read the IPD SMSC data

    # path and loading the data
    SMSCpath=paste(datapath,"/SMSC_phase1_cycles.xlsx",sep="")
    SMSCdata <- read_excel(SMSCpath)
    #keep only needed variables, based on the literature and pre-existing prognostic models on relapses for RRMS patients
    tokeep<-c("unique.visit.id" , "patient.id", "age", "gender","edss", "disease.duration","treatment.naive.prior.visit", "months.since.last.relapse", "nr.relapses.2y.prior.study",
              "nr.Gd.enhanced.lesions","treatment.during.cycle","treatment.time.during.cycle.months","relapse.2y.after.study")
    SMSCdata <- SMSCdata[,(names(SMSCdata) %in% tokeep)]

    ######## Needed transformations for CONTINUOUS variables to approximate the normal distribution #####################

    # disease duration in transformed to log(disease duration +10)
    SMSCdata$disease.duration<-log(SMSCdata$disease.duration+10)
    # months since last relapse relapse is transformed to log(months since last relapse +10)
    SMSCdata$months.since.last.relapse<-log(SMSCdata$months.since.last.relapse+10)

    ######## Needed transformations for CATEGORICAL variables #####################

    # GD enhanced lessions, instead of 0, 1, ..., 18 now we will have 0 and >0
    SMSCdata$nr.Gd.enhanced.lesions[which(SMSCdata$nr.Gd.enhanced.lesions>0)]<-1
    SMSCdata$nr.Gd.enhanced.lesions<-as.factor(SMSCdata$nr.Gd.enhanced.lesions)

    # Number of relapses 2 years prior to study, instead of 0, 1, ..., 6 now we will have 0 and 1 and >1
    SMSCdata$nr.relapses.2y.prior.study[which(SMSCdata$nr.relapses.2y.prior.study>1)]<-2
    SMSCdata$nr.relapses.2y.prior.study<-as.factor(SMSCdata$nr.relapses.2y.prior.study)

    ##make a numeric outcome
    SMSCdata$outcome<-NA
    SMSCdata$outcome[which(SMSCdata$relapse.2y.after.study=="Yes")]<-1
    SMSCdata$outcome[which(SMSCdata$relapse.2y.after.study=="No")]<-0
    SMSCdata$outcome<-as.factor(SMSCdata$outcome)
    #create a variable that indicates the cycle of the patient (i.e. 1st cycle, 2nd cycle, 3rd cycle, etc.)
    SMSCdata$cycle <- ave(SMSCdata$unique.visit.id, SMSCdata$patient.id,  FUN = seq_along)
    #create a column with patient ids as numeric
    SMSCdata$NumericID <- as.numeric(factor(SMSCdata$patient.id,
                                             levels=unique(SMSCdata$patient.id)))


    return(SMSCdata)

}






