dataset=RiskData
library(ggpubr)
#install.packages("gridExtra")
library(gridExtra)

names(dataset)[names(dataset) == "TRT01A"] <- "Treatment"
dataset$Treatment<-recode(dataset$Treatment,"1='Dimethyl fumerate'; 2='Galtiramer acetate'; 3='Natalizumab';4='Placebo'")
names(dataset)[names(dataset) == "logitRisk"] <- "LogitRisk"
##### A. The density of the risk score in the whole dataset
a <- ggplot(dataset, aes(x = LogitRisk))
# y axis scale = ..density.. (default behaviour)
a + geom_density(fill = "lightgray") +
  geom_vline(aes(xintercept = mean(LogitRisk)),
             linetype = "dashed", size = 0.6,color = "#FC4E07")

summary(dataset$LogitRisk) ### mean=-0.5611, median=-0.5991



##### B. The distribution of Risk per arm for each study

##data per study
S<-dataset[which(dataset$STUDYID==1),]
K<-dataset[which(dataset$STUDYID==2),]
L<-dataset[which(dataset$STUDYID==3),]

#install.packages("ggpubr")

# Basic density plot with mean line and marginal rug
ggdensity(dataset, x = "LogitRisk",
          fill = "#0073C2FF", color = "#0073C2FF",
          add = "mean", rug = TRUE)


summary(dataset$LogitRisk) ### mean=-0.5611, median=-0.5991

RiskDist<-ggdensity(dataset, x = "Risk",
          fill = "#0073C2FF", color = "#0073C2FF",
          add = "mean", rug = TRUE)


summary(dataset$Risk) ### mean=0.37106, median=0.35456


# Change outline and fill colors by groups ("Treatment")
# Use a custom palette
#####FOR LOGIT
a<-ggdensity(S, x = "LogitRisk",
          add = "mean", rug = TRUE, xlim=c(-4,4),
          color = "Treatment", fill = "Treatment",
          palette = c("#0073C2FF", "#FC4E07"), xlab = "Logit Risk score for DEFINE study")


b<-ggdensity(K, x = "LogitRisk",
              add = "mean", rug = TRUE,xlim=c(-4,4),
              color = "Treatment", fill = "Treatment",
              palette = c("blue", "red","green"), xlab = "Logit Risk score for CONFIRM study")

c<-ggdensity(L, x = "LogitRisk",
          add = "mean", rug = TRUE,xlim=c(-4,4),
          color = "Treatment", fill = "Treatment",
          palette = c("blue", "yellow"), xlab = "Logit Risk score for AFFIRM study")
ggarrange(a, b, c,
          labels = c("A", "B", "C"),
          ncol = 2, nrow = 2)
###FOR RISK
a<-ggdensity(S, x = "Risk",
             add = "mean", rug = TRUE, xlim=c(0,1),
             color = "Treatment", fill = "Treatment",
             palette = c("#0073C2FF", "#FC4E07"), xlab = "Risk score for DEFINE study")


b<-ggdensity(K, x = "Risk",
             add = "mean", rug = TRUE,xlim=c(0,1),
             color = "Treatment", fill = "Treatment",
             palette = c("blue", "red","green"), xlab = "Risk score for CONFIRM study")

c<-ggdensity(L, x = "Risk",
             add = "mean", rug = TRUE,xlim=c(0,1),
             color = "Treatment", fill = "Treatment",
             palette = c("blue", "yellow"), xlab = "Risk score for AFFIRM study")
RandomizationRisk<-ggarrange(a, b, c,
          labels = c("A", "B", "C"),
          ncol = 2, nrow = 2)



#####The distribution of the risk in those with true relapse and true non-relapse in the entire dataset
dataset$RELAPSE2year<-as.factor(dataset$RELAPSE2year)

ggdensity(dataset, x = "LogitRisk",
          add = "mean", rug = TRUE,
          color = "RELAPSE2year", fill = "RELAPSE2year",
          palette = c("blue", "red"), xlab = "Logit Risk score as a prognostic factor")


PrognosticRisk<-ggdensity(dataset, x = "Risk",
          add = "mean", rug = TRUE,
          color = "RELAPSE2year", fill = "RELAPSE2year",
          palette = c("blue", "red"), xlab = "Risk score as a prognostic factor")


#########
#######logit risk

S$RELAPSE2year<-as.factor(S$RELAPSE2year)
d<-ggdensity(S, x = "LogitRisk",merge=T,
          add = "mean", rug = TRUE,xlim=c(-2,2),
          color = "RELAPSE2year", fill = "Treatment",
          palette = c("blue", "red"), xlab = "Logit Risk score in DEFINE study")
t.test(S$LogitRisk[S$Treatment=="Dimethyl fumerate" & S$RELAPSE2year==0])
t.test(S$LogitRisk[S$Treatment=="Dimethyl fumerate" & S$RELAPSE2year==1])
t.test(S$LogitRisk[S$Treatment=="Placebo" & S$RELAPSE2year==0])
t.test(S$LogitRisk[S$Treatment=="Placebo" & S$RELAPSE2year==1])

summary(S$LogitRisk[S$Treatment=="Dimethyl fumerate" & S$RELAPSE2year==0])
summary(S$LogitRisk[S$Treatment=="Dimethyl fumerate" & S$RELAPSE2year==1])
summary(S$LogitRisk[S$Treatment=="Placebo" & S$RELAPSE2year==0])
summary(S$LogitRisk[S$Treatment=="Placebo" & S$RELAPSE2year==1])


K$RELAPSE2year<-as.factor(K$RELAPSE2year)
e<-ggdensity(K, x = "LogitRisk",merge=T,
          add = "mean", rug = TRUE,xlim=c(-2,2),
          color = "RELAPSE2year", fill = "Treatment",
          palette = c("blue", "red", "yellow"), xlab = "Logit Risk score in CONFIRM study")
t.test(K$LogitRisk[K$Treatment=="Dimethyl fumerate" & K$RELAPSE2year==0])
t.test(K$LogitRisk[K$Treatment=="Dimethyl fumerate" & K$RELAPSE2year==1])
t.test(K$LogitRisk[K$Treatment=="Galtiramer acetate" & K$RELAPSE2year==0])
t.test(K$LogitRisk[K$Treatment=="Galtiramer acetate" & K$RELAPSE2year==1])
t.test(K$LogitRisk[K$Treatment=="Placebo" & K$RELAPSE2year==0])
t.test(K$LogitRisk[K$Treatment=="Placebo" & K$RELAPSE2year==1])

summary(K$LogitRisk[K$Treatment=="Dimethyl fumerate" & K$RELAPSE2year==0])
summary(K$LogitRisk[K$Treatment=="Dimethyl fumerate" & K$RELAPSE2year==1])
summary(K$LogitRisk[K$Treatment=="Galtiramer acetate" & K$RELAPSE2year==0])
summary(K$LogitRisk[K$Treatment=="Galtiramer acetate" & K$RELAPSE2year==1])
summary(K$LogitRisk[K$Treatment=="Placebo" & K$RELAPSE2year==0])
summary(K$LogitRisk[K$Treatment=="Placebo" & K$RELAPSE2year==1])
L$RELAPSE2year<-as.factor(L$RELAPSE2year)
f<-ggdensity(L, x = "LogitRisk",merge=T,
          add = "mean", rug = TRUE, xlim=c(-2,2),
          color = "RELAPSE2year", fill = "Treatment",
          palette = c("blue", "red"), xlab = "Logit Risk score in AFFIRM study")
t.test(L$LogitRisk[L$Treatment=="Natalizumab" & L$RELAPSE2year==0])
t.test(L$LogitRisk[L$Treatment=="Natalizumab" & L$RELAPSE2year==1])
t.test(L$LogitRisk[L$Treatment=="Placebo" & L$RELAPSE2year==0])
t.test(L$LogitRisk[L$Treatment=="Placebo" & L$RELAPSE2year==1])

summary(L$LogitRisk[L$Treatment=="Natalizumab" & L$RELAPSE2year==0])
summary(L$LogitRisk[L$Treatment=="Natalizumab" & L$RELAPSE2year==1])
summary(L$LogitRisk[L$Treatment=="Placebo" & L$RELAPSE2year==0])
summary(L$LogitRisk[L$Treatment=="Placebo" & L$RELAPSE2year==1])


ggarrange(d,e,f, labels = c("A","B","C"))


###Risk
S$RELAPSE2year<-as.factor(S$RELAPSE2year)
d<-ggdensity(S, x = "Risk",merge=T,
             add = "mean", rug = TRUE,xlim=c(0,1),
             color = "RELAPSE2year", fill = "Treatment",
             palette = c("blue", "red"), xlab = "Risk score in DEFINE study")
t.test(S$Risk[S$Treatment=="Dimethyl fumerate" & S$RELAPSE2year==0])
t.test(S$Risk[S$Treatment=="Dimethyl fumerate" & S$RELAPSE2year==1])
t.test(S$Risk[S$Treatment=="Placebo" & S$RELAPSE2year==0])
t.test(S$Risk[S$Treatment=="Placebo" & S$RELAPSE2year==1])

summary(S$Risk[S$Treatment=="Dimethyl fumerate" & S$RELAPSE2year==0])
summary(S$Risk[S$Treatment=="Dimethyl fumerate" & S$RELAPSE2year==1])
summary(S$Risk[S$Treatment=="Placebo" & S$RELAPSE2year==0])
summary(S$Risk[S$Treatment=="Placebo" & S$RELAPSE2year==1])


K$RELAPSE2year<-as.factor(K$RELAPSE2year)
e<-ggdensity(K, x = "Risk",merge=T,
             add = "mean", rug = TRUE,xlim=c(0,1),
             color = "RELAPSE2year", fill = "Treatment",
             palette = c("blue", "red", "yellow"), xlab = "Risk score in CONFIRM study")
t.test(K$Risk[K$Treatment=="Dimethyl fumerate" & K$RELAPSE2year==0])
t.test(K$Risk[K$Treatment=="Dimethyl fumerate" & K$RELAPSE2year==1])
t.test(K$Risk[K$Treatment=="Galtiramer acetate" & K$RELAPSE2year==0])
t.test(K$Risk[K$Treatment=="Galtiramer acetate" & K$RELAPSE2year==1])
t.test(K$Risk[K$Treatment=="Placebo" & K$RELAPSE2year==0])
t.test(K$Risk[K$Treatment=="Placebo" & K$RELAPSE2year==1])

summary(K$Risk[K$Treatment=="Dimethyl fumerate" & K$RELAPSE2year==0])
summary(K$Risk[K$Treatment=="Dimethyl fumerate" & K$RELAPSE2year==1])
summary(K$Risk[K$Treatment=="Galtiramer acetate" & K$RELAPSE2year==0])
summary(K$Risk[K$Treatment=="Galtiramer acetate" & K$RELAPSE2year==1])
summary(K$Risk[K$Treatment=="Placebo" & K$RELAPSE2year==0])
summary(K$Risk[K$Treatment=="Placebo" & K$RELAPSE2year==1])
L$RELAPSE2year<-as.factor(L$RELAPSE2year)
f<-ggdensity(L, x = "Risk",merge=T,
             add = "mean", rug = TRUE, xlim=c(0,1),
             color = "RELAPSE2year", fill = "Treatment",
             palette = c("blue", "red"), xlab = "Risk score in AFFIRM study")
t.test(L$Risk[L$Treatment=="Natalizumab" & L$RELAPSE2year==0])
t.test(L$Risk[L$Treatment=="Natalizumab" & L$RELAPSE2year==1])
t.test(L$Risk[L$Treatment=="Placebo" & L$RELAPSE2year==0])
t.test(L$Risk[L$Treatment=="Placebo" & L$RELAPSE2year==1])

summary(L$Risk[L$Treatment=="Natalizumab" & L$RELAPSE2year==0])
summary(L$Risk[L$Treatment=="Natalizumab" & L$RELAPSE2year==1])
summary(L$Risk[L$Treatment=="Placebo" & L$RELAPSE2year==0])
summary(L$Risk[L$Treatment=="Placebo" & L$RELAPSE2year==1])


EffectModRISK<-ggarrange(d,e,f, labels = c("A","B","C"))
RiskDist
RandomizationRisk
PrognosticRisk
EffectModRISK
