#######################################################################################################
###        Selects variables via LASSO based on half RCTs ################################
###       with 100 cross validations (100 different half-datasets)###########################
#############            1 YEAR RELAPSES           ###################################
################################################################################################################

dataset=MSrelapse
library(glmnet)
library(Hmisc)
####################random half RCTs from studies###############################
for (i in 1:100) {
Advance<-dataset[which(dataset$STUDYID=="ADVANCE"),]
Advance.risk<-Advance[sample(nrow(Advance), nrow(Advance)/2),]
todrop<-c("STUDYID","USUBJID","RELAPSE2year")
Advance.risk<-Advance.risk[ , !(names(Advance.risk) %in% todrop)]

Define<-dataset[which(dataset$STUDYID=="DEFINE"),]
Define.risk<-Define[sample(nrow(Define), nrow(Define)/2),]
Define.risk<-Define.risk[ , !(names(Define.risk) %in% todrop)]

Confirm<-dataset[which(dataset$STUDYID=="CONFIRM"),]
Confirm.risk<-Confirm[sample(nrow(Confirm), nrow(Confirm)/2),]
Confirm.risk<-Confirm.risk[ , !(names(Confirm.risk) %in% todrop)]

Affirm<-dataset[which(dataset$STUDYID=="AFFIRM"),]
Affirm.risk<-Affirm[sample(nrow(Affirm), nrow(Affirm)/2),]
Affirm.risk<-Affirm.risk[ , !(names(Affirm.risk) %in% todrop)]

Mscrg<-dataset[which(dataset$STUDYID=="MSCRG"),]
Mscrg.risk<-Mscrg[sample(nrow(Mscrg), nrow(Mscrg)/2),]
Mscrg.risk<-Mscrg.risk[ , !(names(Mscrg.risk) %in% todrop)]
##all half studies together
mrg<-rbind(Advance.risk,Define.risk,Confirm.risk,Affirm.risk,Mscrg.risk)
#####################LASSO preparation####################

###blinded to treatment so drop variable TRT01A
todrop<-c("TRT01A")
mrg.both<-mrg[ , !(names(mrg) %in% todrop)]
### delete NA values (LASSO requierement)
mrg.both<-na.omit(mrg.both)
#### model matrix needed for LASSO
half.matrix<-model.matrix(mrg.both$RELAPSE1year~.,data=mrg.both)
half.matrix<-na.omit(half.matrix)
#################################LASSO################################
######################################################################
##10 cross validations
cv.fit.half<-cv.glmnet(x=half.matrix,y=mrg.both$RELAPSE1year,family="binomial")
### LASSO coefficients
cv.coef.half<-coef(cv.fit.half,s="lambda.1se")
####RESULTS
### non zero coefficients lead to selected variables
cv.pf.em.half<-rownames(cv.coef.half)[as.numeric(cv.coef.half)!=0]

print(cv.pf.em.half)

}

