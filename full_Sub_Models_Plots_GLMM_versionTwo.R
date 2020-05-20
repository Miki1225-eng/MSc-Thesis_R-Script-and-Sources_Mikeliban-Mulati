#The Role of the Mother in Feeding Skill Acquisition in Immature Sumatran Orangutans_RScript_Miki_31.03
#Supervisors: Dr. Caroline Schuppli; Dr. Robert Hepach 
#########################################################################################################

# set working directory
setwd("C:/Users/Miki/Google Drive/Thesis!!!/Plots & R Script")

# read my full table. 
wild<-read.delim("Subset_wild_ExactAge.txt", header=T, sep="\t",na.strings=c("", "NA"))
str(wild)
wild$AgeClass<-as.factor(wild$AgeClass)
library(psych)
describe(wild$Age)



#IVs and DM
#IVs for the wild data: Exact age of the beggar, sex of the beggar, complexity of the food item, rarity of the food item. 
#DM: success (yes/no).
#random intercepts: 
#Beggar and Date (suggested by Miki)
#Beggar, Mother, time window (half-year based)--age category and it's called 'AgeClass' in the model (suggested by Caroline). 


# handle the missing data 
wild=wild[!is.na(wild$Success),]
wild=wild[!is.na(wild$ProcessingSteps_Complexity),]
wild=wild[!is.na(wild$PopFreq_Rarity),]
wild=wild[!is.na(wild$FollowHours),]
wild=wild[!is.na(wild$AgeClass),]
wild=wild[!is.na(wild$ExactAge),]

#########################################################################

#all the packages I need
library(car)
library(wesanderson)
library(lme4)
library(simr) 
library(influence.ME)
library(MuMIn)
library(ggplot2)
library(Hmisc)
library(GGally)
library(reshape2)
library(compiler)
library(parallel)
library(boot)
library(lattice)
library(sciplot)
mycols<-colors()[c(10,19,72,143,373,419,516,642)]
mycols5=col=c("lightblue","mistyrose","lavender","beige","lightgreen")
mycols2=col=c("lavender","mistyrose")


#wild data:between every variable and the DM:
afs4<- melt(wild[,c("Success","ProcessingSteps_Complexity", "PopFreq_Rarity","ExactAge")],
  id.vars="Success")
ggplot(afs4,aes(factor(Success),y=value,fill=factor(Success))) +
  geom_boxplot() + 
  facet_wrap(~variable,ncol = 5,scales="free_y")
#patterns seem not very contrasting except for complexity variable.
RelationtoB<-spineplot(x=wild$RelationtoB,y=wild$Success,xlab="Relation of the target to the beggar",ylab="Sccess",breaks=lims,
          main="Success by different relationships")
#mother is the most tolerent one and overall the most common target.

#target’s relation to the beggar: because captive ones are with their relatives for most of the time a day in close distance so I expect them to beg from different individuals more often than the wild ones).
library(viridis)
library(hrbrthemes)
library(tidyverse)
library(RColorBrewer)
relatiwild<- table(wild$TargetRelation)
Desire<-pie3D(relatiwild,theta=pi/4,col=c("lightgreen","beige"),explode=0.1,labels=names(relatiwild),main="Target's relation to the beggar:Wild") 
#put these two plots side by side to compare.

##################################################################################################################################
#here are some of the very important steps before running the models:
# Important concepts:
# Linear predictor.
# Link function.
# Error distribution.
# Likelihood and maximum likelihood.
# -> As a consequence we calculate the p-values not by testing the beta-coefficient against '0' but through so-called likelihood ratio tests.
# Deviance (-2LL) and comparing deviances.
# these points above from glm binomial example script by Dr. Hepach (Script R Training seminar May 29th 2019).   


#Model structure suggested by Caroline
full_wild<- glmer(Success ~ ExactAge + Sex + ProcessingSteps_Complexity + PopFreq_Rarity +(1|Beggar)+(1|Mother) + (1|AgeClass) ,data =wild, family = binomial, control = glmerControl(optimizer = "bobyqa"),nAGQ = 1)
#I'm not trying to get inferences from the random intercepts, the singular fit issue can be therefore ignored (Roger Mundry, 2020). 
#print the mod results without correlations among fixed effects
print(full_wild,corr = FALSE)
summary(full_wild)
#significance from complexity only. 
# The fitted values we get from running fitted() are the values for the linear predictor and the fitted values related to the dependent measure are already calculated based on the above equation.
hist(fitted(full_wild))
#looks fine!
plot(wild$Sex,fitted(full_wild), pch=19)
plot(wild$ExactAge,fitted(full_wild), pch=19)
plot(wild$ProcessingSteps_Complexity,fitted(full_wild), pch=19)
plot(wild$PopFreq_Rarity,fitted(full_wild),pch=19)

# To calculate p-values, we proceed as follows:
# Calculate an omnibus model comparison test.  
# Note that this step is not always necessary, especially if there are a priori predictions with 
regards to predictor variables. Here we simply want to acknowledge the concept of an omnibus test.
# The null model essentialy predicts the same value for each focal. 
nullModelwild= glmer(Success ~ 1+(1|Beggar)+(1|Mother)+(1|AgeClass), family = binomial, data = wild)
# Through comparing the two models we can address the following important question: Do the predictors have a combined incluence on the dependent measure.
anova(nullModelwild,full_wild,test="Chisq")
# here we see an improve in model fit!
# There is one function that allows me to do both model comparisons at once:
drop1(full_wild, test="Chisq")
#drop1--> to get each of the predictors at once

#significance in Age and Complexity!
# Get the the 95% confidence interval for the beta-coefficients:
confint(full_wild)
cbind(coefficients(full_wild),confint(full_wild))

#very large effect from rarity but very little from age, which is strange since from the plot afs4 there seems no difference in success in rarity. 
lattice::dotplot(ranef(full_wild,which = "Beggar", condVar = TRUE), scales = list(y = list(alternating = 0)))
lattice::dotplot(ranef(full_wild,which = "Mother", condVar = TRUE), scales = list(y = list(alternating = 0)))
lattice::dotplot(ranef(full_wild,which = "AgeClass", condVar = TRUE), scales = list(y = list(alternating = 0)))

#Model structure suggested by Miki
full_wild<- glmer(Success ~ Age + Sex+RelationtoB+ ProcessingSteps_Complexity+PopFreq_Rarity_Onetofive+(1|Begger)+(1|Date),data = wild,family = binomial,nAGQ = 1)
# print the mod results without correlations among fixed effects
print(full_wild,corr = FALSE)
anova(full_wild)
summary(full_wild)
#significance from age, sex, and complexity.
lattice::dotplot(ranef(full_wild,which = "Begger", condVar = TRUE), scales = list(y = list(alternating = 0)))
lattice::dotplot(ranef(full_wild,which = "Date", condVar = TRUE), scales = list(y = list(alternating = 0)))

#get confidence intervals (CIs) and rough estimates using the SEs
seswild<-sqrt(diag(vcov(full_wild)))
# table of estimates with 95% CI
tabwild<- cbind(Est= fixef(full_wild), LL = fixef(full_wild) - 1.96 * seswild, UL = fixef(full_wild) + 1.96 *seswild)
exp(tabwild)


#####################################################################################################################################
#plots
ggpredict(full_wild,"ExactAge")
ggpredict(full_wild,"Sex")
ggpredict(full_wild,"ProcessingSteps_Complexity")
ggpredict(full_wild,"PopFreq_Rarity")
iv1wild<- ggpredict(full_wild,"ExactAge")
plot(iv1wild)
#to plot the raw data:ExactAge and Success
ggboxplot(wild,x="Success",y="ExactAge",color ="Success",palette = "jco")
#seems more yes to the older ones but not that different.
iv2wild<- ggpredict(full_wild,"Sex")
plot(iv2wild)
plot(wild$Sex, wild$Success, xlab="Sex", ylab="Success")
#even a bit larger difference but not significantly.
iv3wild<- ggpredict(full_wild,"ProcessingSteps_Complexity")
plot(iv3wild)
#this supports my hypothesis!
ggboxplot(wild,x="Success",y="ProcessingSteps_Complexity",color ="Success",palette = "jco")
iv4wild<- ggpredict(full_wild,"PopFreq_Rarity")
plot(iv4wild)
#this supports my hypothesis but not significant enough.
ggboxplot(wild,x="Success",y="PopFreq_Rarity",color ="Success",palette = "jco")



