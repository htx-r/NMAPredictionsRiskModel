############################ Script for reproduce paper's figures and tables ######################################################



############################   TABLES   ##############################################################################
#### Table 2
x<-LASSOModel$lassomodel$coefficients
y<-FabioModel$fabiomodel$coefficients
x<-as.data.frame(x)
colnames(x)<-c("Coefficients")
rownames(x)<-c("Intercept", "Age", "Baseline Weight", "Baseline EDSS","No. of relapses 1 year prior to study", "Prior MS Treatment group", "REGION: India","REGION: North America","REGION: ROW", "REGION: Western Europe", "Baseline Gd volume", "Baseline SF-36 PCS", "Baseline Actual Distance Walked >500")
x<-round(x,3)
y<-c("Intercept", "Age", "Baseline Weight", "Baseline EDSS","No. of relapses 1 year prior to study", "Prior MS Treatment group", "REGION: India","REGION: North America","REGION: ROW", "REGION: Western Europe", "Baseline Gd volume", "Baseline SF-36 PCS", "Baseline Actual Distance Walked >500")
y<-as.data.frame(y)
Table_LASSOModel<-cbind(y,x)
rownames(Table_LASSOModel) <- c()
colnames(Table_LASSOModel)<-c("Variables","Coefficients")
Table_LASSOModel

y<-FabioModel$fabiomodel$coefficients
y<-as.data.frame(y)
colnames(y)<-c("Coefficients")
y<-round(y,3)
z<-c("Intercept", "Age", "Sex (male vs female)", "Baseline EDSS","Years Since Onset of Symptoms","Ethnicity (white vs other)","No. of relapses 1 year prior to study","Months since pre-study relapse", "Prior MS Treatment group (yes vs no)","T25FW","9HPT","PASAT-3","VFT 2.5%", "Baseline SF-36 MCS", "Baseline SF-36 PCS")
z<-as.data.frame(z)
Table_PreSpecifiedModel<-cbind(z,y)
colnames(Table_PreSpecifiedModel)<-c("Variables","Coefficients")
rownames(Table_PreSpecifiedModel) <- c()
Table_PreSpecifiedModel


#### Table 3
gamma_LASSO<-IPDNMRJAGSmodelLASSO$BUGSoutput$summary[9,]
deltaDF_LASSO<-IPDNMRJAGSmodelLASSO$BUGSoutput$summary[5,]
deltaGA_LASSO<-IPDNMRJAGSmodelLASSO$BUGSoutput$summary[6,]
deltaN_LASSO<-IPDNMRJAGSmodelLASSO$BUGSoutput$summary[7,]
gammaDF_LASSO<-IPDNMRJAGSmodelLASSO$BUGSoutput$summary[10,]
gammaGA_LASSO<-IPDNMRJAGSmodelLASSO$BUGSoutput$summary[11,]
gammaN_LASSO<-IPDNMRJAGSmodelLASSO$BUGSoutput$summary[12,]


LASSOIPDNMR_Table<-rbind(gamma_LASSO,deltaDF_LASSO,deltaGA_LASSO,deltaN_LASSO,gammaDF_LASSO,gammaGA_LASSO,gammaN_LASSO)
LASSOIPDNMR_Table<-as.data.frame(LASSOIPDNMR_Table)
todrop<-c(2,4,5,6,8,9)
LASSOIPDNMR_Table<-LASSOIPDNMR_Table[,-todrop]
LASSOIPDNMR_Table<-round(LASSOIPDNMR_Table,2)


LASSOIPDNMR_Table$CredibleIntervals<-NA
for(i in 1:7){
  LASSOIPDNMR_Table[i,4]<-paste(LASSOIPDNMR_Table[i,1],"(",LASSOIPDNMR_Table[i,2],",", LASSOIPDNMR_Table[i,3], ")")
}
todrop<-c(1,2,3)
LASSOIPDNMR_Table<-LASSOIPDNMR_Table[,-todrop]
LASSOIPDNMR_Table<-as.data.frame(LASSOIPDNMR_Table)

gamma_Prespecified<-IPDNMRJAGSmodelFabio$BUGSoutput$summary[9,]
deltaDF_Prespecified<-IPDNMRJAGSmodelFabio$BUGSoutput$summary[5,]
deltaGA_Prespecified<-IPDNMRJAGSmodelFabio$BUGSoutput$summary[6,]
deltaN_Prespecified<-IPDNMRJAGSmodelFabio$BUGSoutput$summary[7,]
gammaDF_Prespecified<-IPDNMRJAGSmodelFabio$BUGSoutput$summary[10,]
gammaGA_Prespecified<-IPDNMRJAGSmodelFabio$BUGSoutput$summary[11,]
gammaN_Prespecified<-IPDNMRJAGSmodelFabio$BUGSoutput$summary[12,]


PrespecifiedIPDNMR_Table<-rbind(gamma_Prespecified,deltaDF_Prespecified,deltaGA_Prespecified,deltaN_Prespecified,gammaDF_Prespecified,gammaGA_Prespecified,gammaN_Prespecified)
PrespecifiedIPDNMR_Table<-as.data.frame(PrespecifiedIPDNMR_Table)
todrop<-c(2,4,5,6,8,9)
PrespecifiedIPDNMR_Table<-PrespecifiedIPDNMR_Table[,-todrop]

PrespecifiedIPDNMR_Table<-round(PrespecifiedIPDNMR_Table,2)
PrespecifiedIPDNMR_Table$CredibleIntervals<-NA
for(i in 1:7){
PrespecifiedIPDNMR_Table[i,4]<-paste(PrespecifiedIPDNMR_Table[i,1],"(",PrespecifiedIPDNMR_Table[i,2],",", PrespecifiedIPDNMR_Table[i,3], ")")
}

todrop<-c(1,2,3)
PrespecifiedIPDNMR_Table<-PrespecifiedIPDNMR_Table[,-todrop]
PrespecifiedIPDNMR_Table<-as.data.frame(PrespecifiedIPDNMR_Table)

IPDNMR_Table<-cbind(LASSOIPDNMR_Table,PrespecifiedIPDNMR_Table)
colnames(IPDNMR_Table)<-c("LASSO model Mean (95% Cr. Interval)", "Pre-specified model Mean (95% Cr. Interval)" )
rownames(IPDNMR_Table)<-c("??0","??_DF","??_GA","??_N", "??_DF","??_GA","??_N")

#### Table 4
###absolute benefits
LASSOtable
Fabiotable
###ORs
LASSOtableOR
FabiotableOR


############################   FIGURES   ##############################################################################
#### Figure 1

par(mfrow=c(1,2))
#### LASSO model
# Visual assessment of calibration by risk groups
# create 10 risk groups
expit<-function(x) {exp(x)/(1+exp(x))}
predprob<-predict(RiskModel)
predprob<-expit(predprob)
groups <- cut(predprob,breaks=quantile(predprob, prob = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)),labels=c(1:10),include.lowest=TRUE)
RELAPSE2year<-RiskData1$RELAPSE2year
# average the observed and expected probabilities of patients in each risk group
gpdata <- cbind(RiskData1,groups,predprob)
obs <- ddply(gpdata,~groups,summarise,mean=mean(as.numeric(RELAPSE2year)))[,2]-1
exp <- ddply(gpdata,~groups,summarise,mean=mean(predprob))
obsn <- table(RiskData1$RELAPSE2year,groups)[2,]


# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# graph a simple calibration plot over 10 risk groups
par(pty="s")
plot(obs~exp[,2],xlim=c(0,1),ylim=c(0,1),ylab="Observed",xlab="Expected")
title("A  LASSO model", adj=0, cex.main=0.8)
lines(c(0,1),c(0,1),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]))
}
h <- hist(predprob, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(RELAPSE2year))~predprob,span=1))-1
lines_data <- data.frame(predprob,obs_all)
lines_data2 <- lines_data[order(predprob),]
lines(lines_data2[,1],lines_data2[,2])

#### Pre-specified model


dataset<-MSrelapse
X<-dataset
X<-na.omit(X)
RiskModel1<-FabioModel$fabiomodel

# Visual assessment of calibration by risk groups
# create 10 risk groups
expit<-function(x) {exp(x)/(1+exp(x))}
predprob<-predict(RiskModel1)
predprob<-expit(predprob)
groups <- cut(predprob,breaks=quantile(predprob, prob = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)),labels=c(1:10),include.lowest=TRUE)
RELAPSE2year<-RiskData1$RELAPSE2year
# average the observed and expected probabilities of patients in each risk group
gpdata <- cbind(RiskData1,groups,predprob)
obs <- ddply(gpdata,~groups,summarise,mean=mean(as.numeric(RELAPSE2year)))[,2]-1
exp <- ddply(gpdata,~groups,summarise,mean=mean(predprob))
obsn <- table(RiskData1$RELAPSE2year,groups)[2,]


# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# graph a simple calibration plot over 10 risk groups
par(pty="s")
plot(obs~exp[,2],xlim=c(0,1),ylim=c(0,1),ylab="Observed",xlab="Expected")
title("B   Pre-specified model", adj=0, cex.main=0.8)
lines(c(0,1),c(0,1),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]))
}
h <- hist(predprob, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(RELAPSE2year))~predprob,span=1))-1
lines_data <- data.frame(predprob,obs_all)
lines_data2 <- lines_data[order(predprob),]
lines(lines_data2[,1],lines_data2[,2])


#### Figure 2

ggarrange(PrognosticRiskLASSO,PrognosticRiskFabio,ncol = 1,nrow = 2,labels = c("A  LASSO model","B  Pre-specified model"), hjust=-0.2,font.label = list(size = 11))


#### Figure 3

ggarrange(IPDplotLASSO,IPDplotFabio,ncol = 1,nrow = 2,labels = c("A","B"),font.label = list(size = 11))

################################## Appendix figures ############################################


### ORs plot
ggarrange(IPDplotLASSO_OR,IPDplotFabio_OR,ncol = 1,nrow = 2,labels = c("A","B"),font.label = list(size = 11))

### Flow-chart Appendix

library(DiagrammeR)


grViz("digraph flowchart {
      # node definitions with substituted label text
      node [fontname = Helvetica, shape = rectangle]        
      tab1 [label = '@@1']
      tab2 [label = '@@2']
      tab3 [label = '@@3']
      tab4 [label = '@@4']
      tab5 [label = '@@5']
      # edge definitions with the node IDs
      tab1 -> tab2 -> tab3 ;
      tab3->tab4;
      tab3->tab5;
      
      }
      
      [1]: 'Number of prognostic factors, np=57'
      [2]: 'Number of prognostic factors with missing data less than 50%, np=53'
      [3]: 'Number of prognostic factors correlated less than 70%, np=31 '
      [4]: 'Number of prognostic factors in LASSO model, np=12 '
      [5]: 'Number of prognostic factors in pre-specified model, np=14 '
      
      ")
