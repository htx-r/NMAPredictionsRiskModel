library(MCMCglmm)


SMSCdata_glmer<-SMSCdataC
SMSCdata_glmer$outcome<-as.integer(factor(SMSCdata_glmer$relapse.2y.after.study))-1
SMSCdata_glmer$nr.Gd.enhanced.lesions<-as.numeric(SMSCdata_glmer$nr.Gd.enhanced.lesions)-1
SMSCdata_glmer$treatment.naive.prior.visit<-as.numeric(SMSCdata_glmer$treatment.naive.prior.visit)
SMSCdata_glmer$gender<-as.numeric(as.factor(SMSCdata_glmer$gender))-1
SMSCdata_glmer$nr.relapses.2y.prior.study<-as.numeric(SMSCdata_glmer$nr.relapses.2y.prior.study)-1


MCMCmodel <- MCMCglmm(as.integer(outcome) ~  age + disease.duration + edss + nr.Gd.enhanced.lesions + nr.relapses.2y.prior.study+
                        months.since.last.relapse + treatment.naive.prior.visit + gender + treatment.time.during.cycle.months , 
                      random = ~idh(age+disease.duration+edss+nr.Gd.enhanced.lesions+nr.relapses.2y.prior.study+months.since.last.relapse
                                    +treatment.naive.prior.visit+gender+treatment.time.during.cycle.months):patient.id +
                    + patient.id, data = SMSCdata_glmer, verbose = F, family="categorical")


 prior.m2b.1 = list(R = list(V = 1, fix = 1), G = list(G1 = list(V = 1, nu = 0.002)))
 m2b.1 <- MCMCglmm(Pupated ~ 1, random = ~FSfamily,
                     family = "categorical", data = PlodiaRB, prior = prior.m2b.1,
                     verbose = FALSE)
 