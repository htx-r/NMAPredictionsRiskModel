
# Visual assessment of calibration by risk groups

# create 10 risk groups
RiskDataOmitNa<-na.omit(RiskData)
groups <- cut(RiskDataOmitNa$Risk,breaks=quantile(RiskDataOmitNa$Risk, prob = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)),labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of patients in each risk group 
gpdata <- cbind(RiskDataOmitNa,groups)
obs <- (ddply(gpdata,~groups,summarise,mean=mean(as.numeric(DAY30)))[,2])-1
exp <- ddply(gpdata,~groups,summarise,mean=mean(pred_prob))
obsn <- table(DAY30,groups)[1,] 

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
h <- hist(pred_prob, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(DAY30)-1)~pred_prob,span=1))
lines_data <- data.frame(pred_prob,obs_all)
lines_data2 <- lines_data[order(pred_prob),] 
lines(lines_data2[,1],lines_data2[,2],col="blue")
legend(0.0,0.9,c("Risk groups","Reference line","95% CI","Loess"),col=c("red","black","green","blue"),lty=c(0,2,1,1),pch=c(1,NA,NA,NA),bty="n")
