################# Script for Plot of predicted Risk #########################

########################################## FOR LASSO MODEL   ######################################
###predicted risk= expit(logit p risk)
expit<-function(x) {exp(x)/(1+exp(x))}
p<-IPDNMRJAGSmodelLASSO$BUGSoutput$mean$logitp
p<-expit(IPDNMRJAGSmodelLASSO$BUGSoutput$mean$logitp)

####preperation for the graph
###For Dymethyl fumarate- Risk & propability to relapse
DF<-cbind(Risknew,p[,1])
DF$Treatment<-1
colnames(DF)<-c("Risknew","prelapse","Treatment")
##For Glatimarer acetate  - Risk & propability to relapse
GA<-cbind(Risknew,p[,2])
GA$Treatment<-2
colnames(GA)<-c("Risknew","prelapse","Treatment")
##For Natalizumab acetate - Risk & propability to relapse
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

IPDplotLASSO<-ggplot(Graphdata, aes(x=Risknew, y=prelapse, group=Treatment)) +
  geom_line(aes(color=Treatment))+
  geom_point(aes(color=Treatment))
IPDplotLASSO

######################################################################################################

################################# For Fabio's model #################################################

###predicted risk= expit(logit p risk)
expit<-function(x) {exp(x)/(1+exp(x))}
p<-IPDNMRJAGSmodelFabio$BUGSoutput$mean$logitp
p<-expit(IPDNMRJAGSmodelFabio$BUGSoutput$mean$logitp)

####preperation for the graph
###For Dymethyl fumarate- Risk & propability to relapse
DF<-cbind(Risknew,p[,1])
DF$Treatment<-1
colnames(DF)<-c("Risknew","prelapse","Treatment")
##For Glatimarer acetate  - Risk & propability to relapse
GA<-cbind(Risknew,p[,2])
GA$Treatment<-2
colnames(GA)<-c("Risknew","prelapse","Treatment")
##For Natalizumab acetate - Risk & propability to relapse
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

IPDplotFabio<-ggplot(Graphdata, aes(x=Risknew, y=prelapse, group=Treatment)) +
  geom_line(aes(color=Treatment))+
  geom_point(aes(color=Treatment))
IPDplotFabio


