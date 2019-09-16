##################### Script ro Graph the results from NMR IPD model ############


###model for the logit of RELAPSE for placebo patients only
##only placebo arms are needed
PlaceboData<-RiskData[which(RiskData$TRT01A==4),]

modelFORplacebo<-function(){
  ###likelihood
  for (i in 1:Np){
    outcome[i]~dbern(pplacebo[i])
    ###formula
    logit(pplacebo[i])<-logitpplacebo
  }

  logitpplacebo~dnorm(0,0.01)
}

jagsdataFORplacebo <- list(
  Np=nrow(PlaceboData),
  outcome=PlaceboData$RELAPSE2year
)

#run the model
IPDFORplaceboJAGSmodel <- jags.parallel(data = jagsdataFORplacebo,inits=NULL,parameters.to.save = c('logitpplacebo','pplacebo'),model.file = modelFORplacebo,
                                        n.chains=2,n.iter = 100000,n.burnin = 1000,DIC=F,n.thin = 10)
print(IPDFORplaceboJAGSmodel)


##prepare the predictions
Risknew<-seq(0.01,0.99,0.01)

Risknew<-as.data.frame(Risknew)
logit<-function(x) {log(x/(1-x))}
expit<-function(x) {exp(x)/(1+exp(x))}

logitRisknew<-NA
logitRisknew<-as.data.frame(logitRisknew)
for (i in 1:99) {
  logitRisknew[i,1]<-logit(Risknew[i,1])
}


logitRisknew<-as.data.frame(logitRisknew)


logitp<-NA
logitp<-as.data.frame(logitp)

### I did not substract the mean(logit) as it is equal to 0 - meanRisk=0.5 & meanlogitrisk=0 (logit(0.5)=0)
for (i in 1:99){
  logitp[i,1]<-0.117-0.917+1.219*(logitRisknew[i,1])+0.131*(logitRisknew[i,1])
  logitp[i,2]<-0.117-0.699+1.219*(logitRisknew[i,1])-0.082*(logitRisknew[i,1])
  logitp[i,3]<-0.117-1.337+1.219*(logitRisknew[i,1])-0.390*(logitRisknew[i,1])
  logitp[i,4]<-0.117+1.219*(logitRisknew[i,1])
}

colnames(logitp)<-c("Treatment1","Treatment2","Treatment3","Treatment4")

p<-NA
p<-as.data.frame(p)
p<-expit(logitp)

######### Create the Data for the graph

##For Dimethyl fumarate - Risk & propability to relapse
DF<-cbind(Risknew,p[,1])
DF$Treatment<-1
colnames(DF)<-c("Risknew","prelapse","Treatment")
##For Glatiramer acetate - Risk & propability to relapse
GA<-cbind(Risknew,p[,2])
GA$Treatment<-2
colnames(GA)<-c("Risknew","prelapse","Treatment")
##For Natalizumab - Risk & propability to relapse
N<-cbind(Risknew,p[,3])
N$Treatment<-3
colnames(N)<-c("Risknew","prelapse","Treatment")
##For Placebo - Risk & propability to relapse
Pl<-cbind(Risknew,p[,4])
Pl$Treatment<-4
colnames(Pl)<-c("Risknew","prelapse","Treatment")
##merge data for all the treatments
Graphdata<-rbind(DF,GA,N,Pl)

##### graph

#library(ggplot2)
#library(ggpubr)
#install.packages("gridExtra")
#library(gridExtra)
# Basic line plot with points
ggplot(data=Graphdata, aes(x=Risknew, y=prelapse, group=Treatment)) +
  geom_line()+
  geom_point()

IPDplot<-ggplot(Graphdata, aes(x=Risknew, y=prelapse, group=Treatment)) +
  geom_line(aes(color=Treatment))+
  geom_point(aes(color=Treatment))
IPDplot

#remove no needed items
rm(DF)
rm(GA)
rm(N)
rm(Pl)
rm(logitp)
rm(logitRisknew)
rm(p)
rm(PlaceboData)
rm(Risknew)
rm(IPDFORplaceboJAGSmodel)
