######################### SCRIPT PLOTTING RISK & LOGITRISK DISTRIBUTIONS  ###########


#library(ggpubr)
#install.packages("gridExtra")
#library(gridExtra)

####col.names of dataset
dataset=RiskData
names(dataset)[names(dataset) == "TRT01A"] <- "Treatment"
dataset$Treatment<-recode(dataset$Treatment,"1='Dimethyl fumarate';2='Galtiramer acetate'; 3='Natalizumab';4='Placebo'")
names(dataset)[names(dataset) == "logitRiskLASSO"] <- "LogitRiskLASSO"
names(dataset)[names(dataset) == "logitRiskFabio"] <- "LogitRiskFabio"
##### A. The density of the risk score in the whole dataset
#a <- ggplot(dataset, aes(x = LogitRisk))
# y axis scale = ..density.. (default behaviour)
#a + geom_density(fill = "lightgray") +
#  geom_vline(aes(xintercept = mean(LogitRisk)),
#             linetype = "dashed", size = 0.6,color = "#FC4E07")

#################################################LASSO model########################################
#################################################################################################

# A. The distribution in the whole dataset

#logitrisk
summary(dataset$LogitRiskLASSO) ### mean=-0.5611, median=-0.5991
LogitRiskLASSO<-ggdensity(dataset, x = "LogitRiskLASSO",
          fill = "#0073C2FF", color = "#0073C2FF",
          add = "mean", rug = TRUE)
#risk
summary(dataset$RiskLASSO)
RiskDistLASSO<-ggdensity(dataset, x = "RiskLASSO",
                    fill = "#0073C2FF", color = "#0073C2FF",
                    add = "mean", rug = TRUE)


##### B. The distribution per arm for each study

##data per study
S<-dataset[which(dataset$STUDYID==1),]
K<-dataset[which(dataset$STUDYID==2),]
L<-dataset[which(dataset$STUDYID==3),]


#####FOR LOGIT
a<-ggdensity(S, x = "LogitRiskLASSO",
          add = "mean", rug = TRUE, xlim=c(-4,4),
          color = "Treatment", fill = "Treatment",
          palette = c("blue", "red"), xlab = "Logit Risk score for DEFINE study")


b<-ggdensity(K, x = "LogitRiskLASSO",
              add = "mean", rug = TRUE,xlim=c(-4,4),
              color = "Treatment", fill = "Treatment",
              palette = c("blue", "yellow","red"), xlab = "Logit Risk score for CONFIRM study")

c<-ggdensity(L, x = "LogitRiskLASSO",
          add = "mean", rug = TRUE,xlim=c(-4,4),
          color = "Treatment", fill = "Treatment",
          palette = c("grey", "red"), xlab = "Logit Risk score for AFFIRM study")

RandomizationLogitRiskLASSO<-ggarrange(a,b,c,
          labels = c("A", "B", "C"),
          ncol = 2, nrow = 2)
###FOR RISK
a<-ggdensity(S, x = "RiskLASSO",
             add = "mean", rug = TRUE, xlim=c(0,1),
             color = "Treatment", fill = "Treatment",
             palette = c("blue", "red"), xlab = "Risk score for DEFINE study")


b<-ggdensity(K, x = "RiskLASSO",
             add = "mean", rug = TRUE,xlim=c(0,1),
             color = "Treatment", fill = "Treatment",
             palette = c("blue", "yellow","red"), xlab = "Risk score for CONFIRM study")

c<-ggdensity(L, x = "RiskLASSO",
             add = "mean", rug = TRUE,xlim=c(0,1),
             color = "Treatment", fill = "Treatment",
             palette = c("grey", "red"), xlab = "Risk score for AFFIRM study")

RandomizationRiskLASSO<-ggarrange(a, b, c,
          labels = c("A", "B", "C"),
          ncol = 2, nrow = 2)



#C. The distributionin those with true relapse and true non-relapse in the entire dataset
dataset$RELAPSE2year<-as.factor(dataset$RELAPSE2year)
#logitRisk
PrognosticLogitRiskLASSO<-ggdensity(dataset, x = "LogitRiskLASSO",
          add = "mean", rug = TRUE,
          color = "RELAPSE2year", fill = "RELAPSE2year",
          palette = c("blue", "yellow"), xlab = "Logit Risk score as a prognostic factor")

#Risk
PrognosticRiskLASSO<-ggdensity(dataset, x = "RiskLASSO",
          add = "mean", rug = TRUE,
          color = "RELAPSE2year", fill = "RELAPSE2year",
          palette = c("blue", "yellow"), xlab = "Risk score as a prognostic factor")
##Risk as prognostic factor
summary(dataset$RiskLASSO[dataset$RELAPSE2year==1])
summary(dataset$RiskLASSO[dataset$RELAPSE2year==0])
t.test(dataset$RiskLASSO[dataset$RELAPSE2year==1])
t.test(dataset$RiskLASSO[dataset$RELAPSE2year==0])
#D. The distribution in those with true relapse and true non-relapse for each arm in each study

##logit risk
S$RELAPSE2year<-as.factor(S$RELAPSE2year)
d<-ggdensity(S, x = "LogitRiskLASSO",merge=T,
          add = "mean", rug = TRUE,xlim=c(-2,2),
          color = "RELAPSE2year", fill = "Treatment",
          palette = c("blue", "yellow"), xlab = "Logit Risk score in DEFINE study")
t.test(S$LogitRiskLASSO[S$Treatment=="Dimethyl fumarate" & S$RELAPSE2year==0])
t.test(S$LogitRiskLASSO[S$Treatment=="Dimethyl fumarate" & S$RELAPSE2year==1])
t.test(S$LogitRiskLASSO[S$Treatment=="Placebo" & S$RELAPSE2year==0])
t.test(S$LogitRiskLASSO[S$Treatment=="Placebo" & S$RELAPSE2year==1])

summary(S$LogitRiskLASSO[S$Treatment=="Dimethyl fumarate" & S$RELAPSE2year==0])
summary(S$LogitRiskLASSO[S$Treatment=="Dimethyl fumarate" & S$RELAPSE2year==1])
summary(S$LogitRiskLASSO[S$Treatment=="Placebo" & S$RELAPSE2year==0])
summary(S$LogitRiskLASSO[S$Treatment=="Placebo" & S$RELAPSE2year==1])


K$RELAPSE2year<-as.factor(K$RELAPSE2year)
e<-ggdensity(K, x = "LogitRiskLASSO",merge=T,
          add = "mean", rug = TRUE,xlim=c(-2,2),
          color = "RELAPSE2year", fill = "Treatment",
          palette = c("blue", "yellow", "red"), xlab = "Logit Risk score in CONFIRM study")
t.test(K$LogitRiskLASSO[K$Treatment=="Dimethyl fumarate" & K$RELAPSE2year==0])
t.test(K$LogitRiskLASSO[K$Treatment=="Dimethyl fumarate" & K$RELAPSE2year==1])
t.test(K$LogitRiskLASSO[K$Treatment=="Galtiramer acetate" & K$RELAPSE2year==0])
t.test(K$LogitRiskLASSO[K$Treatment=="Galtiramer acetate" & K$RELAPSE2year==1])
t.test(K$LogitRiskLASSO[K$Treatment=="Placebo" & K$RELAPSE2year==0])
t.test(K$LogitRiskLASSO[K$Treatment=="Placebo" & K$RELAPSE2year==1])

summary(K$LogitRiskLASSO[K$Treatment=="Dimethyl fumarate" & K$RELAPSE2year==0])
summary(K$LogitRiskLASSO[K$Treatment=="Dimethyl fumarate" & K$RELAPSE2year==1])
summary(K$LogitRiskLASSO[K$Treatment=="Galtiramer acetate" & K$RELAPSE2year==0])
summary(K$LogitRiskLASSO[K$Treatment=="Galtiramer acetate" & K$RELAPSE2year==1])
summary(K$LogitRiskLASSO[K$Treatment=="Placebo" & K$RELAPSE2year==0])
summary(K$LogitRiskLASSO[K$Treatment=="Placebo" & K$RELAPSE2year==1])
L$RELAPSE2year<-as.factor(L$RELAPSE2year)
f<-ggdensity(L, x = "LogitRiskLASSO",merge=T,
          add = "mean", rug = TRUE, xlim=c(-2,2),
          color = "RELAPSE2year", fill = "Treatment",
          palette = c("blue", "yellow"), xlab = "Logit Risk score in AFFIRM study")
t.test(L$LogitRiskLASSO[L$Treatment=="Natalizumab" & L$RELAPSE2year==0])
t.test(L$LogitRiskLASSO[L$Treatment=="Natalizumab" & L$RELAPSE2year==1])
t.test(L$LogitRiskLASSO[L$Treatment=="Placebo" & L$RELAPSE2year==0])
t.test(L$LogitRiskLASSO[L$Treatment=="Placebo" & L$RELAPSE2year==1])

summary(L$LogitRiskLASSO[L$Treatment=="Natalizumab" & L$RELAPSE2year==0])
summary(L$LogitRiskLASSO[L$Treatment=="Natalizumab" & L$RELAPSE2year==1])
summary(L$LogitRiskLASSO[L$Treatment=="Placebo" & L$RELAPSE2year==0])
summary(L$LogitRiskLASSO[L$Treatment=="Placebo" & L$RELAPSE2year==1])


EffectModLogitRiskLASSO<-ggarrange(d,e,f, labels = c("A","B","C"))


###Risk
S$RELAPSE2year<-as.factor(S$RELAPSE2year)
d<-ggdensity(S, x = "RiskLASSO",merge=T,
             add = "mean", rug = TRUE,xlim=c(0,1),
             color = "RELAPSE2year", fill = "Treatment",
             palette = c("blue","red"), xlab = "Risk score in DEFINE study")
t.test(S$RiskLASSO[S$Treatment=="Dimethyl fumarate" & S$RELAPSE2year==0])
t.test(S$RiskLASSO[S$Treatment=="Dimethyl fumarate" & S$RELAPSE2year==1])
t.test(S$RiskLASSO[S$Treatment=="Placebo" & S$RELAPSE2year==0])
t.test(S$RiskLASSO[S$Treatment=="Placebo" & S$RELAPSE2year==1])

summary(S$RiskLASSO[S$Treatment=="Dimethyl fumarate" & S$RELAPSE2year==0])
summary(S$RiskLASSO[S$Treatment=="Dimethyl fumarate" & S$RELAPSE2year==1])
summary(S$RiskLASSO[S$Treatment=="Placebo" & S$RELAPSE2year==0])
summary(S$RiskLASSO[S$Treatment=="Placebo" & S$RELAPSE2year==1])


K$RELAPSE2year<-as.factor(K$RELAPSE2year)
e<-ggdensity(K, x = "RiskLASSO",merge=T,
             add = "mean", rug = TRUE,xlim=c(0,1),
             color = "RELAPSE2year", fill = "Treatment",
             palette = c("blue", "yellow", "red"), xlab = "Risk score in CONFIRM study")
t.test(K$RiskLASSO[K$Treatment=="Dimethyl fumarate" & K$RELAPSE2year==0])
t.test(K$RiskLASSO[K$Treatment=="Dimethyl fumarate" & K$RELAPSE2year==1])
t.test(K$RiskLASSO[K$Treatment=="Galtiramer acetate" & K$RELAPSE2year==0])
t.test(K$RiskLASSO[K$Treatment=="Galtiramer acetate" & K$RELAPSE2year==1])
t.test(K$RiskLASSO[K$Treatment=="Placebo" & K$RELAPSE2year==0])
t.test(K$RiskLASSO[K$Treatment=="Placebo" & K$RELAPSE2year==1])

summary(K$RiskLASSO[K$Treatment=="Dimethyl fumarate" & K$RELAPSE2year==0])
summary(K$RiskLASSO[K$Treatment=="Dimethyl fumarate" & K$RELAPSE2year==1])
summary(K$RiskLASSO[K$Treatment=="Galtiramer acetate" & K$RELAPSE2year==0])
summary(K$RiskLASSO[K$Treatment=="Galtiramer acetate" & K$RELAPSE2year==1])
summary(K$RiskLASSO[K$Treatment=="Placebo" & K$RELAPSE2year==0])
summary(K$RiskLASSO[K$Treatment=="Placebo" & K$RELAPSE2year==1])
L$RELAPSE2year<-as.factor(L$RELAPSE2year)
f<-ggdensity(L, x = "RiskLASSO",merge=T,
             add = "mean", rug = TRUE, xlim=c(0,1),
             color = "RELAPSE2year", fill = "Treatment",
             palette = c("gray", "red"), xlab = "Risk score in AFFIRM study")
t.test(L$RiskLASSO[L$Treatment=="Natalizumab" & L$RELAPSE2year==0])
t.test(L$RiskLASSO[L$Treatment=="Natalizumab" & L$RELAPSE2year==1])
t.test(L$RiskLASSO[L$Treatment=="Placebo" & L$RELAPSE2year==0])
t.test(L$RiskLASSO[L$Treatment=="Placebo" & L$RELAPSE2year==1])

summary(L$RiskLASSO[L$Treatment=="Natalizumab" & L$RELAPSE2year==0])
summary(L$RiskLASSO[L$Treatment=="Natalizumab" & L$RELAPSE2year==1])
summary(L$RiskLASSO[L$Treatment=="Placebo" & L$RELAPSE2year==0])
summary(L$RiskLASSO[L$Treatment=="Placebo" & L$RELAPSE2year==1])


EffectModRiskLASSO<-ggarrange(d,e,f,labels = c("A","B","C"))



################################################# Fabio's model ########################################
#################################################################################################

# A. The distribution in the whole dataset

#logitrisk
summary(dataset$LogitRiskFabio) ### mean=-0.5611, median=-0.5991
LogitRiskLASSOFabio<-ggdensity(dataset, x = "LogitRiskFabio",
                          fill = "#0073C2FF", color = "#0073C2FF",
                          add = "mean", rug = TRUE)
#risk
summary(dataset$RiskFabio)
RiskDistFabio<-ggdensity(dataset, x = "RiskFabio",
                         fill = "#0073C2FF", color = "#0073C2FF",
                         add = "mean", rug = TRUE)


##### B. The distribution per arm for each study

##data per study
S<-dataset[which(dataset$STUDYID==1),]
K<-dataset[which(dataset$STUDYID==2),]
L<-dataset[which(dataset$STUDYID==3),]


#####FOR LOGIT
a<-ggdensity(S, x = "LogitRiskFabio",
             add = "mean", rug = TRUE, xlim=c(-4,4),
             color = "Treatment", fill = "Treatment",
             palette = c("blue", "red"), xlab = "Logit Risk score for DEFINE study")


b<-ggdensity(K, x = "LogitRiskFabio",
             add = "mean", rug = TRUE,xlim=c(-4,4),
             color = "Treatment", fill = "Treatment",
             palette = c("blue", "yellow","red"), xlab = "Logit Risk score for CONFIRM study")

c<-ggdensity(L, x = "LogitRiskFabio",
             add = "mean", rug = TRUE,xlim=c(-4,4),
             color = "Treatment", fill = "Treatment",
             palette = c("grey", "red"), xlab = "Logit Risk score for AFFIRM study")

RandomizationLogitRiskFabio<-ggarrange(a,b,c,
                                       labels = c("A", "B", "C"),
                                       ncol = 2, nrow = 2)
###FOR RISK
a<-ggdensity(S, x = "RiskFabio",
             add = "mean", rug = TRUE, xlim=c(0,1),
             color = "Treatment", fill = "Treatment",
             palette = c("blue", "red"), xlab = "Risk score for DEFINE study")


b<-ggdensity(K, x = "RiskFabio",
             add = "mean", rug = TRUE,xlim=c(0,1),
             color = "Treatment", fill = "Treatment",
             palette = c("blue", "yellow","red"), xlab = "Risk score for CONFIRM study")

c<-ggdensity(L, x = "RiskFabio",
             add = "mean", rug = TRUE,xlim=c(0,1),
             color = "Treatment", fill = "Treatment",
             palette = c("grey", "red"), xlab = "Risk score for AFFIRM study")

RandomizationRiskFabio<-ggarrange(a, b, c,
                                  labels = c("A", "B", "C"),
                                  ncol = 2, nrow = 2)



#C. The distributionin those with true relapse and true non-relapse in the entire dataset
dataset$RELAPSE2year<-as.factor(dataset$RELAPSE2year)
#logitRisk
PrognosticLogitRiskFabio<-ggdensity(dataset, x = "LogitRiskFabio",
                                    add = "mean", rug = TRUE,
                                    color = "RELAPSE2year", fill = "RELAPSE2year",
                                    palette = c("blue", "yellow"), xlab = "Logit Risk score as a prognostic factor")

#Risk
PrognosticRiskFabio<-ggdensity(dataset, x = "RiskFabio",
                               add = "mean", rug = TRUE,
                               color = "RELAPSE2year", fill = "RELAPSE2year",
                               palette = c("blue", "yellow"), xlab = "Risk score as a prognostic factor")
##Risk as prognostic factor
summary(dataset$RiskFabio[dataset$RELAPSE2year==1])
summary(dataset$RiskFabio[dataset$RELAPSE2year==0])
t.test(dataset$RiskFabio[dataset$RELAPSE2year==1])
t.test(dataset$RiskFabio[dataset$RELAPSE2year==0])
#D. The distribution in those with true relapse and true non-relapse for each arm in each study

##logit risk
S$RELAPSE2year<-as.factor(S$RELAPSE2year)
d<-ggdensity(S, x = "LogitRiskFabio",merge=T,
             add = "mean", rug = TRUE,xlim=c(-2,2),
             color = "RELAPSE2year", fill = "Treatment",
             palette = c("blue", "yellow"), xlab = "Logit Risk score in DEFINE study")
t.test(S$LogitRiskFabio[S$Treatment=="Dimethyl fumarate" & S$RELAPSE2year==0])
t.test(S$LogitRiskFabio[S$Treatment=="Dimethyl fumarate" & S$RELAPSE2year==1])
t.test(S$LogitRiskFabio[S$Treatment=="Placebo" & S$RELAPSE2year==0])
t.test(S$LogitRiskFabio[S$Treatment=="Placebo" & S$RELAPSE2year==1])

summary(S$LogitRiskFabio[S$Treatment=="Dimethyl fumarate" & S$RELAPSE2year==0])
summary(S$LogitRiskFabio[S$Treatment=="Dimethyl fumarate" & S$RELAPSE2year==1])
summary(S$LogitRiskFabio[S$Treatment=="Placebo" & S$RELAPSE2year==0])
summary(S$LogitRiskFabio[S$Treatment=="Placebo" & S$RELAPSE2year==1])


K$RELAPSE2year<-as.factor(K$RELAPSE2year)
e<-ggdensity(K, x = "LogitRiskFabio",merge=T,
             add = "mean", rug = TRUE,xlim=c(-2,2),
             color = "RELAPSE2year", fill = "Treatment",
             palette = c("blue", "yellow", "red"), xlab = "Logit Risk score in CONFIRM study")
t.test(K$LogitRiskFabio[K$Treatment=="Dimethyl fumarate" & K$RELAPSE2year==0])
t.test(K$LogitRiskFabio[K$Treatment=="Dimethyl fumarate" & K$RELAPSE2year==1])
t.test(K$LogitRiskFabio[K$Treatment=="Galtiramer acetate" & K$RELAPSE2year==0])
t.test(K$LogitRiskFabio[K$Treatment=="Galtiramer acetate" & K$RELAPSE2year==1])
t.test(K$LogitRiskFabio[K$Treatment=="Placebo" & K$RELAPSE2year==0])
t.test(K$LogitRiskFabio[K$Treatment=="Placebo" & K$RELAPSE2year==1])

summary(K$LogitRiskFabio[K$Treatment=="Dimethyl fumarate" & K$RELAPSE2year==0])
summary(K$LogitRiskFabio[K$Treatment=="Dimethyl fumarate" & K$RELAPSE2year==1])
summary(K$LogitRiskFabio[K$Treatment=="Galtiramer acetate" & K$RELAPSE2year==0])
summary(K$LogitRiskFabio[K$Treatment=="Galtiramer acetate" & K$RELAPSE2year==1])
summary(K$LogitRiskFabio[K$Treatment=="Placebo" & K$RELAPSE2year==0])
summary(K$LogitRiskFabio[K$Treatment=="Placebo" & K$RELAPSE2year==1])
L$RELAPSE2year<-as.factor(L$RELAPSE2year)
f<-ggdensity(L, x = "LogitRiskFabio",merge=T,
             add = "mean", rug = TRUE, xlim=c(-2,2),
             color = "RELAPSE2year", fill = "Treatment",
             palette = c("blue", "yellow"), xlab = "Logit Risk score in AFFIRM study")
t.test(L$LogitRiskFabio[L$Treatment=="Natalizumab" & L$RELAPSE2year==0])
t.test(L$LogitRiskFabio[L$Treatment=="Natalizumab" & L$RELAPSE2year==1])
t.test(L$LogitRiskFabio[L$Treatment=="Placebo" & L$RELAPSE2year==0])
t.test(L$LogitRiskFabio[L$Treatment=="Placebo" & L$RELAPSE2year==1])

summary(L$LogitRiskFabio[L$Treatment=="Natalizumab" & L$RELAPSE2year==0])
summary(L$LogitRiskFabio[L$Treatment=="Natalizumab" & L$RELAPSE2year==1])
summary(L$LogitRiskFabio[L$Treatment=="Placebo" & L$RELAPSE2year==0])
summary(L$LogitRiskFabio[L$Treatment=="Placebo" & L$RELAPSE2year==1])


EffectModLogitRiskFabio<-ggarrange(d,e,f, labels = c("A","B","C"))


###Risk
S$RELAPSE2year<-as.factor(S$RELAPSE2year)
d<-ggdensity(S, x = "RiskFabio",merge=T,
             add = "mean", rug = TRUE,xlim=c(0,1),
             color = "RELAPSE2year", fill = "Treatment",
             palette = c("blue","red"), xlab = "Risk score in DEFINE study")
t.test(S$RiskFabio[S$Treatment=="Dimethyl fumarate" & S$RELAPSE2year==0])
t.test(S$RiskFabio[S$Treatment=="Dimethyl fumarate" & S$RELAPSE2year==1])
t.test(S$RiskFabio[S$Treatment=="Placebo" & S$RELAPSE2year==0])
t.test(S$RiskFabio[S$Treatment=="Placebo" & S$RELAPSE2year==1])

summary(S$RiskFabio[S$Treatment=="Dimethyl fumarate" & S$RELAPSE2year==0])
summary(S$RiskFabio[S$Treatment=="Dimethyl fumarate" & S$RELAPSE2year==1])
summary(S$RiskFabio[S$Treatment=="Placebo" & S$RELAPSE2year==0])
summary(S$RiskFabio[S$Treatment=="Placebo" & S$RELAPSE2year==1])


K$RELAPSE2year<-as.factor(K$RELAPSE2year)
e<-ggdensity(K, x = "RiskFabio",merge=T,
             add = "mean", rug = TRUE,xlim=c(0,1),
             color = "RELAPSE2year", fill = "Treatment",
             palette = c("blue", "yellow", "red"), xlab = "Risk score in CONFIRM study")
t.test(K$RiskFabio[K$Treatment=="Dimethyl fumarate" & K$RELAPSE2year==0])
t.test(K$RiskFabio[K$Treatment=="Dimethyl fumarate" & K$RELAPSE2year==1])
t.test(K$RiskFabio[K$Treatment=="Galtiramer acetate" & K$RELAPSE2year==0])
t.test(K$RiskFabio[K$Treatment=="Galtiramer acetate" & K$RELAPSE2year==1])
t.test(K$RiskFabio[K$Treatment=="Placebo" & K$RELAPSE2year==0])
t.test(K$RiskFabio[K$Treatment=="Placebo" & K$RELAPSE2year==1])

summary(K$RiskFabio[K$Treatment=="Dimethyl fumarate" & K$RELAPSE2year==0])
summary(K$RiskFabio[K$Treatment=="Dimethyl fumarate" & K$RELAPSE2year==1])
summary(K$RiskFabio[K$Treatment=="Galtiramer acetate" & K$RELAPSE2year==0])
summary(K$RiskFabio[K$Treatment=="Galtiramer acetate" & K$RELAPSE2year==1])
summary(K$RiskFabio[K$Treatment=="Placebo" & K$RELAPSE2year==0])
summary(K$RiskFabio[K$Treatment=="Placebo" & K$RELAPSE2year==1])
L$RELAPSE2year<-as.factor(L$RELAPSE2year)
f<-ggdensity(L, x = "RiskFabio",merge=T,
             add = "mean", rug = TRUE, xlim=c(0,1),
             color = "RELAPSE2year", fill = "Treatment",
             palette = c("gray", "red"), xlab = "Risk score in AFFIRM study")
t.test(L$RiskFabio[L$Treatment=="Natalizumab" & L$RELAPSE2year==0])
t.test(L$RiskFabio[L$Treatment=="Natalizumab" & L$RELAPSE2year==1])
t.test(L$RiskFabio[L$Treatment=="Placebo" & L$RELAPSE2year==0])
t.test(L$RiskFabio[L$Treatment=="Placebo" & L$RELAPSE2year==1])

summary(L$RiskFabio[L$Treatment=="Natalizumab" & L$RELAPSE2year==0])
summary(L$RiskFabio[L$Treatment=="Natalizumab" & L$RELAPSE2year==1])
summary(L$RiskFabio[L$Treatment=="Placebo" & L$RELAPSE2year==0])
summary(L$RiskFabio[L$Treatment=="Placebo" & L$RELAPSE2year==1])


EffectModRiskFabio<-ggarrange(d,e,f,labels = c("A","B","C"))

boxplot(RiskData$RiskLASSO,RiskData$RiskFabio,col=c("red","blue"), names = c("LASSO","Fabio"))
#remove no needed items
rm(a)
rm(b)
rm(c)
rm(d)
rm(e)
rm(f)
rm(K)
rm(L)
rm(S)

