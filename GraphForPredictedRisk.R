################# Script for Plot of predicted Risk #########################

###predicted risk= expit(logit p risk)
expit<-function(x) {exp(x)/(1+exp(x))}
p<-IPDNMRJAGSmodel$BUGSoutput$mean$logitp
p<-expit(IPDNMRJAGSmodel$BUGSoutput$mean$logitp)

####preperation for the graph
###For Avonex - Risk & propability to relapse
AV<-cbind(Risknew,p[,1])
AV$Treatment<-1
colnames(AV)<-c("Risknew","prelapse","Treatment")
##For Dymethyl fumarate - Risk & propability to relapse
DF<-cbind(Risknew,p[,2])
DF$Treatment<-2
colnames(DF)<-c("Risknew","prelapse","Treatment")
##For Glatiramer acetate - Risk & propability to relapse
GA<-cbind(Risknew,p[,3])
GA$Treatment<-3
colnames(GA)<-c("Risknew","prelapse","Treatment")
##For Natalizumab - Risk & propability to relapse
N<-cbind(Risknew,p[,4])
N$Treatment<-4
colnames(N)<-c("Risknew","prelapse","Treatment")
##For Placebo - Risk & propability to relapse
Pl<-cbind(Risknew,p[,5])
Pl$Treatment<-5
colnames(Pl)<-c("Risknew","prelapse","Treatment")
##merge data for all the treatments
Graphdata<-rbind(AV,DF,GA,N,Pl)

###Graph for IPD
ggplot(data=Graphdata, aes(x=Risknew, y=prelapse, group=Treatment)) +
  geom_line()+
  geom_point()

IPDplot<-ggplot(Graphdata, aes(x=Risknew, y=prelapse, group=Treatment)) +
  geom_line(aes(color=Treatment))+
  geom_point(aes(color=Treatment))
IPDplot


