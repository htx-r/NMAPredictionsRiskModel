
###model for the logit of RELAPSE for placebo patients only
##only placebo arms are needed
K<-RiskData[which(RiskData$TRT01A==4),]

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
  Np=nrow(K),
  outcome=K$RELAPSE2year
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

GraphData<-cbind(Risknew,p)

m<-cbind(Risknew,p[,1])
m$Treatment<-1
colnames(m)<-c("Risknew","prelapse","Treatment")

n<-cbind(Risknew,p[,2])
n$Treatment<-2
colnames(n)<-c("Risknew","prelapse","Treatment")
k<-cbind(Risknew,p[,3])
k$Treatment<-3
colnames(k)<-c("Risknew","prelapse","Treatment")
l<-cbind(Risknew,p[,4])
l$Treatment<-4
colnames(l)<-c("Risknew","prelapse","Treatment")
Graphdata<-rbind(m,n,k,l)

##### graph

library(ggplot2)
library(ggpubr)
#install.packages("gridExtra")
library(gridExtra)
# Basic line plot with points
ggplot(data=Graphdata, aes(x=Risknew, y=prelapse, group=Treatment)) +
  geom_line()+
  geom_point()

IPDplot<-ggplot(Graphdata, aes(x=Risknew, y=prelapse, group=Treatment)) +
  geom_line(aes(color=Treatment))+
  geom_point(aes(color=Treatment))
IPDplot

# Use custom color palettes
# p+scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
# Use brewer color palettes
#p+scale_color_brewer(palette="Dark2")
# Use grey scale
#p + scale_color_grey() + theme_classic()


