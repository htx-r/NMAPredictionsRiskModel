################# Script for Plot of predicted Risk #########################

###predicted risk= expit(logit p risk)
expit<-function(x) {exp(x)/(1+exp(x))}
p<-IPDNMRJAGSmodel$BUGSoutput$mean$logitp
p<-expit(IPDNMRJAGSmodel$BUGSoutput$mean$logitp)

####preperation for the graph
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

###Graph for IPD
ggplot(data=Graphdata, aes(x=Risknew, y=prelapse, group=Treatment)) +
  geom_line()+
  geom_point()

IPDplot<-ggplot(Graphdata, aes(x=Risknew, y=prelapse, group=Treatment)) +
  geom_line(aes(color=Treatment))+
  geom_point(aes(color=Treatment))
IPDplot

