##################SCRIPT for the calibration plot ##########
Calibrationplots.fun<-function(model){
if (model=="LASSOModel"){
dataset<-MSrelapse
X<-dataset
X<-na.omit(X)
RiskModel1<-LASSOModel$lassomodel
RiskData1<-X

#### LASSO model
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
plot(obs~exp[,2],xlim=c(0,1),ylim=c(0,1),col="red",ylab="Observed",xlab="Expected")
lines(c(0,1),c(0,1),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col="green")
}
h <- hist(predprob, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(RELAPSE2year))~predprob,span=1))-1
lines_data <- data.frame(predprob,obs_all)
lines_data2 <- lines_data[order(predprob),]
lines(lines_data2[,1],lines_data2[,2],col="blue")
legend(0.0,0.9,c("Risk groups","Reference line","95% CI","Loess"),col=c("red","black","green","blue"),lty=c(0,2,1,1),pch=c(1,NA,NA,NA),bty="n")
}


if (model=="FabioModel"){
  dataset<-MSrelapse
  X<-dataset
  X<-na.omit(X)
  RiskModel1<-FabioModel$fabiomodel
  RiskData1<-X
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
  plot(obs~exp[,2],xlim=c(0,1),ylim=c(0,1),col="red",ylab="Observed",xlab="Expected")
  lines(c(0,1),c(0,1),lty=2)
  for(i in 1:10){
    lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col="green")
  }
  h <- hist(predprob, breaks=50, plot=FALSE)
  for(i in 1:length(h$mids)){
    lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],1-((h$counts[i]/max(h$counts))/10)))
  }

  # Add a loess smoother to the plot
  obs_all <- predict(loess((as.numeric(RELAPSE2year))~predprob,span=1))-1
  lines_data <- data.frame(predprob,obs_all)
  lines_data2 <- lines_data[order(predprob),]
  lines(lines_data2[,1],lines_data2[,2],col="blue")
  legend(0.0,0.9,c("Risk groups","Reference line","95% CI","Loess"),col=c("red","black","green","blue"),lty=c(0,2,1,1),pch=c(1,NA,NA,NA),bty="n")
}

rm(dataset)
rm(exp)


rm(lines_data)
rm(lines_data2)

rm(RiskData1)
rm(RiskModel1)
rm(X)
rm(groups)
rm(i)
rm(lci)
rm(obs)
rm(obsn)
rm(obs_all)
rm(predprob)
rm(RELAPSE2year)

rm(uci)
rm(gpdata)
rm(h)
}
