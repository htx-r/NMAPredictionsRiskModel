###################################################################################################################
####################  Function that checks the result of Frequentist Generalized Mixed-Effects Model ###########
###################################################################################################################


Frequentist_glmm.fun=function(dataset){

#treatment is not included ata that moment as for complete cases it has only one category (Yes)
  #also categorical variable were not included as random effects, as the glmer did no run with an error

  glmer_out<-glmer(outcome ~ (1|patient.id) + age + disease.duration + edss + nr.Gd.enhanced.lesions + nr.relapses.2y.prior.study +
                     months.since.last.relapse + treatment.naive.prior.visit + gender + treatment.time.during.cycle.months +
                     (-1+age|patient.id)+(-1+treatment.time.during.cycle.months|patient.id)+(-1+edss|patient.id)+(-1+months.since.last.relapse|patient.id),
                   data=dataset, family = binomial)


  glm_out<-glm((outcome)~ age + disease.duration + edss + nr.Gd.enhanced.lesions + nr.relapses.2y.prior.study +
                 months.since.last.relapse + treatment.naive.prior.visit + gender + treatment.time.during.cycle.months,
               data=dataset, family = binomial)

  return(list(glm_out=glm_out,glmer_out=glmer_out))

}



