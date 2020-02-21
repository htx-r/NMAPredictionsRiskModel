##########################################################################################
########################### CHECK THE RESULTS OF IPD and AD WITH NETMETA ##################
################################################################################



dataAD<-cbind(c(1,1,2,2,2,3,3,4,4,5,5),c(1,4,1,2,4,3,4,4,2,4,2),c(343,172,309,161,157,574,281,25,25,126,125),c(100,93,92,53,73,173,157,19,11,97,89))
dataAD<-as.data.frame(dataAD)
colnames(dataAD)<-c("St","Dr","N_ran","N_rel")

TestPair <- pairwise(treat=Dr, event=N_rel, n=N_ran, data=dataAD, sm="OR", studlab=St, allstudies = TRUE)

net1 <- netmeta(TE, seTE, treat1, treat2, studlab, data = TestPair, sm = "OR", comb.random=FALSE, comb.fixed=TRUE, prediction=TRUE, ref=4)

summary(net1, digits = 2)

forest(net1,ref=4,fontsize=10)
