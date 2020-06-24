#For Begging Frequencies & Begging Rates Computations Only: 4 Gaussian models in total
#################################################################################################

# set working directory
setwd("C:/Users/Miki/Google Drive/Thesis!!!/Plots & R Script")

# read my full table for the purpose of analyses (full_ana). 
#The dofference between this table and my original full table is that this table has no P culumns or comment or missing values written as NA (just empty instead)!
con<-read.delim("Condensed.txt", header=T, sep="\t",na.strings=c("", "NA"))
str(con) 
half<-read.delim("HalfYearTable.txt", header=T, sep="\t",na.strings=c("", "NA"))
str(half)
zoocon<-read.delim("Zoo condensed.txt", header=T, sep="\t",na.strings=c("", "NA"))
str(zoocon) 
conzoowild<-read.delim("Condensed_zoowildcompar.txt", header=T, sep="\t",na.strings=c("", "NA"))
str(conzoowild)

#########################################################################################################
#all the packages I need
library(Hmisc)
library(corrplot)
library(PerformanceAnalytics)
library(devtools)
library(ggpubr)
library(coin)
library(car)
library(pgirmess)
library(psych)
library(wesanderson)
library(lme4)
library(simr) 
library(influence.ME)
library(MuMIn)
library(ggplot2)
library(GGally)
library(reshape2)
library(compiler)
library(parallel)
library(boot)
library(lattice)
library(sciplot)
library(plotrix)
library(viridis)
library(hrbrthemes)
library(tidyverse)
library(RColorBrewer)
library(magrittr)
library(ggeffects)
library(effects)
library(sjmisc)
library(splines)


#there are four Gaussian models with three being begging frequency and one being begging rate
###########################################Following the lecture slides on Gaussian models by Roger Mundry on February 5, 2020#################################################################
#model 1: Begging Frequency ~ Age + Sex + Complexity +Site (conzoowild)(vi)
# handle the missing data 
conzoowild=conzoowild[!is.na(conzoowild$AgeClass),]
conzoowild=conzoowild[!is.na(conzoowild$Sex),]
conzoowild=conzoowild[!is.na(conzoowild$Complexity),]
conzoowild=conzoowild[!is.na(conzoowild$Site),]
conzoowild=conzoowild[!is.na(conzoowild$FoodItem),]
conzoowild=conzoowild[!is.na(conzoowild$BeggingFrequency),]
str(conzoowild)
conzoowild$BeggerAge=as.factor(paste(conzoowild$Begger,round(conzoowild$AgeClass)))
#check the disribution of the DM
hist(conzoowild$BeggingFrequency)
describe(conzoowild$BeggingFrequency)
#very skewed so let's log transform it
hist(log(conzoowild$BeggingFrequency))
describe(log(conzoowild$BeggingFrequency))
#this is okay ish, keep it like this
conzoowild$log.BeggingFrequency=log(conzoowild$BeggingFrequency)
hist(sqrt(conzoowild$AgeClass))
describe(sqrt(conzoowild$AgeClass))
#since not that different, keep the original
hist(conzoowild$Complexity)
describe(conzoowild$Complexity)
#keep it as it is
#z transform all the covariates
conzoowild$z.AgeClass=as.vector(scale(conzoowild$AgeClass))
conzoowild$z.Complexity=as.vector(scale(conzoowild$Complexity))
describe(conzoowild$z.AgeClass)
describe(conzoowild$AgeClass)
#just keep the original ones for now, use those z-transformed ones in case the model doesn't converge.
#random intercepts are mother/begger beggerage and fooditem
full=lmer(log.BeggingFrequency~z.AgeClass+z.Complexity+Sex+Site+(1|BeggerAge)+(1|FoodItem),data=conzoowild,REML=F)
full.1=lmer(log.begfreq~AgeClass+Complexity+Sex+Site+(1|FoodItem),data=conzoowild,REML=F)
#collineraty
xx=lm(log.begfreq~z.AgeClass+z.Complexity+Sex+Site,data=conzoowild)
vif(xx)
#there is no issue
summary(full)
#let's see all the effects --> to get LRT, df, and p values.
drop1(full, test="Chisq")
drop1(full.1, test="Chisq")
null.1<-lmer(log.begfreq~(1|FoodItem), data=conzoowild,REML=F)
anova(null.1, full.1, test="Chisq")
summary(full.1)
#maybe keep full.1?
#to get estimates and SE
round(summary(full.1)$coefficients,3)
round(summary(full)$coefficients,3)
#to get sd. for random effects
summary(full.1)$varcor
summary(full)$varcor
#to get CIs
confint.merMod(object=full.1)
confint.merMod(object=full)
#it does work!
#plot the significant results: Site
plot(allEffects(full), ylab="log(Probability of Begging Frequency)", selection=4)

############################################################################################################
#model 2 & 3: begging frequency ~ Age + Sex + Complexity + Rarity (iv)& begging rate (vii)~ same effects (con) 
# handle the missing data 
con=con[!is.na(con$BeggingFrequency),]
con=con[!is.na(con$BeggingRate),]
con=con[!is.na(con$AgeClass),]
con=con[!is.na(con$Sex),]
con=con[!is.na(con$Complexity),]
con=con[!is.na(con$Rarity),]
con=con[!is.na(con$FoodItem),]
str(con)
con$BeggerAge=as.factor(paste(con$Begger,round(con$AgeClass)))
#check the disribution of the DM
hist(con$BeggingFrequency)
describe(con$BeggingFrequency)
hist(con$BeggingRate)
describe(con$BeggingRate)
#both very skewed so let's log transform it
hist(log(con$BeggingFrequency))
describe(log(con$BeggingFrequency))
hist(log(con$BeggingRate))
describe(log(con$BeggingRate))
#inspect the distribution of all the covariates 
hist(con$AgeClass)
describe(con$AgeClass)
describe(sqrt(con$AgeClass))
con$sqrt.Age=sqrt(con$AgeClass)
hist(con$Complexity)
describe(con$Complexity)
describe(con$Rarity)
describe(log(con$Rarity))
#keep the log one!
con$log.Rarity=log(con$Rarity)
#this is okay ish, keep it like these!
con$log.BeggingFrequency=log(con$BeggingFrequency)
con$log.BeggingRate=log(con$BeggingRate)
full.bf=lmer(log.BeggingFrequency~sqrt.Age+Complexity+Sex+log.Rarity+(1|BeggerAge)+(1|FoodItem), data=con,REML=F)
full.br=lmer(log.BeggingRate~sqrt.Age+Complexity+Sex+log.Rarity+(1|BeggerAge)+(1|FoodItem), data=con,REML=F)
#both converged, so there is no need to z.transformation.
#collineraty
xx1=lm(log.BeggingFrequency~sqrt.Age+Complexity+Sex+log.Rarity,data=con)
vif(xx1)
x2=lm(log.BeggingRate~sqrt.Age+Complexity+Sex+log.Rarity,data=con)
vif(x2)
#there is no issue!
#to get LRT, df, and p values
drop1(full.bf, test="Chisq")
drop1(full.br, test="Chisq")
#to get estimates and SE
round(summary(full.bf)$coefficients,3)
round(summary(full.br)$coefficients,3)
#to get sd. for random effects
summary(full.bf)$varcor
summary(full.br)$varcor
#to get CIs
confint.merMod(object=full.bf)
confint.merMod(object=full.br)
#it does work!
#plot the significant results: Complexity and Rarity
plot(allEffects(full.br), ylab="log(Probability of Begging Rate)", selection=2)
plot(allEffects(full.br), ylab="log(Probability of Begging Rate)", selection=4)

##############################################################################################################################
#model 4: begging frequency (~ age + Sex + Complexity + Rarity + Desirability (zoocon)(v)
zoocon=zoocon[!is.na(zoocon$Rarity),]
str(zoocon)
#check the disribution of the DM
hist(zoocon$BeggingFreq)
describe(zoocon$BeggingFreq)
#both very skewed so let's log transform it
hist(log(zoocon$BeggingFreq))
describe(log(zoocon$BeggingFreq))
#this is okay ish, keep it like these!
zoocon$log.BeggingFrequency=log(zoocon$BeggingFreq)
hist(zoocon$AgeClass)
describe(zoocon$AgeClass)) #just keep it although not ideal
hist(zoocon$Complexity)
describe(zoocon$Complexity)
hist(zoocon$Rarity)
describe(zoocon$Rarity)
hist(zoocon$Desirability)
describe(zoocon$Desirability)
full=lmer(log.BeggingFrequency~AgeClass+Complexity+Sex+Rarity+Desirability+(1|Species_Item_Simple), data=zoocon,REML=F)
#collineraty
xx=lm(log.BeggingFrequency~AgeClass+Complexity+Sex+Rarity+Desirability,data=zoocon)
vif(xx)
#there is no issue!
#get LRT, df, and p values.
drop1(full, test="Chisq")
#to get estimates and SE
round(summary(full)$coefficients,3)
#to get sd. for random effects
summary(full)$varcor
#to get CIs
confint.merMod(object=full)
#it does work!
#plot the significant results: Complexity and Rarity
plot(allEffects(full), ylab="log(Probability of Begging Frequency)", selection=1)
plot(allEffects(full), ylab="log(Probability of Begging Frequency)", selection=4)

#one thing keep in mind though: there are only three different age classes so it's prelimitary with the age effect found...







