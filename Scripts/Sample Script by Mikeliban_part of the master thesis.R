# Title: Sample Script by Mikeliban_part of the master thesis.R
# Author: Mulati Mikeliban
# Date: 24. June 2020
# Title of master thesis: "The Role of the Mother in Feeding Skill Acquisition in Immature Sumatran Orangutans"
# Master thesis supervisors: Dr. Caroline Schuppli; Dr. Robert Hepach 
# Script contains a part of the full analysis from the master thesis: generalized linear mixed model with binomial outcome using wild sample only

#******************************************
#*                                        *
#*  Setup                                 *
#*                                        *
#******************************************
# set working directory
setwd("")

# load needed packages
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
library(effects)
library(Matrix)
library(ggeffects)
library(stargazer)
library(psych)

# read in dataset
wild<-read.delim("Subset_wild_ExactAge.txt", header=T, sep="\t", na.strings=c("", "NA"))

#******************************************
#*                                        *
#*  Preparatory Steps                     *
#*                                        *
#******************************************
str(wild)

# focus on the begging target being the mother only
wild_mo<-subset(wild,RelationtoB=="mother")
str(wild_mo)

# handle missing data [only for running the model!]
wild_mo=wild_mo[!is.na(wild_mo$Success),]
wild_mo=wild_mo[!is.na(wild_mo$ExactAge),]
wild_mo=wild_mo[!is.na(wild_mo$FoodItem),]
wild_mo=wild_mo[!is.na(wild_mo$Mother),]
wild_mo=wild_mo[!is.na(wild_mo$Complexity),]
wild_mo=wild_mo[!is.na(wild_mo$Rarity),]
wild_mo=wild_mo[!is.na(wild_mo$Begger.id),]
str(wild_mo)

# get distribution of all the covariates
ggpairs(wild_mo[,c("ExactAge","FollowHours","Complexity","Rarity")])
#-> rarity and age variables seem very skewed. 

# random intercepts in the model: focal ID ("Begger.id"), "Mother", food item ID "Food item", & BeggerAge (i.e. per subject at different age in days)
# check if data is balanced in each of those
table(wild_mo$Begger.id)
table(wild_mo$Mother)
table(wild_mo$FoodItem)

# check the frequency of occurrence of success per chimp ID and food item
table(wild_mo$Success,wild_mo$Begger.id)
table(wild_mo$Success,wild_mo$FoodItem)

# check the distributions of the covariates: ExactAge, AgeClass, Complexity, Rarity 
hist(wild_mo$ExactAge)
hist(wild_mo$AgeClass)
hist(wild_mo_mo$Complexity)
hist(wild_mo$Rarity)
# only the complexity is okay-ish.

# inspect cross-tabulation of food items and begger id: could there be a special taste of an individual?
table(wild_mo$FoodItem,wild_mo$Begger.id)
table(wild_mo$Begger.id,wild_mo$ExactAge)
table(wild_mo$Begger.id,wild_mo$Mother)
# -> they are both partly crossed


#******************************************
#*                                        *
#*  GLM-Model                             *
#*                                        *
#******************************************
# fit the model
contr=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000))
wild_mo$BeggerAge=as.factor(paste(wild_mo$Begger.id,round(wild_mo$ExactAge)))
hist(wild_mo$ExactAge)
describe(wild_mo$ExactAge)
hist(sqrt(wild_mo$ExactAge))
describe(sqrt(wild_mo$ExactAge))
# nice, keep the 'sqrt.Age' because skewness is nearly 0.

hist(wild_mo$Complexity)
describe(wild_mo$Complexity)
# keep it as it is because skewness is nearly 0.

hist(wild_mo$Rarity)
describe(wild_mo$Rarity)
hist(sqrt(wild_mo$Rarity))
describe(sqrt(wild_mo$Rarity))

# keep it as 'sqrt.Rarity'because skewness is nearly 0.
wild_mo$sqrt.Rarity=sqrt(wild_mo$Rarity)
wild_mo$sqrt.Age=sqrt(wild_mo$ExactAge)
interaction<- glmer(Success ~ sqrt.Age + Sex + sqrt.Age * Sex + Complexity + sqrt.Rarity+I(sqrt.Age^2) + 
                   (1|BeggerAge) + (1|Mother/Begger.id) + 
                   (0 + sqrt.Age | Mother/Begger.id) +(1|FoodItem) + 
                   (0 + Complexity | FoodItem)+(0 + sqrt.Rarity| FoodItem),family=binomial(link="logit"),
                   data=wild_mo,control=contr)

# singular fit issue
summary(interaction)$varcor 
ll.old=logLik(interaction)
round(ll.old, 3); logLik(interaction)

# exclude random slopes for ExactAge
interaction.s<- glmer(Success ~ ExactAge + Sex + ExactAge*Sex+Complexity + Rarity +(1|BeggerAge)+ (1|FoodItem)+ (0+Complexity|FoodItem)+(0+Rarity|FoodItem),family=binomial(link="logit"), data=wild_mo, control=contr)

# it worked! Let's see if they are actually very different
summary(interaction.s)$varcor  
logLik(interaction); logLik(interaction.s)
round(summary(interaction)$coefficients, 3)
round(summary(interaction.s)$coefficients, 3)
# they are exactly the same so just keep the original one since it's more informative.

# collinearity
wild_mo$Success=as.numeric(wild_mo$Success)
co=lm(Success ~ ExactAge + Sex +Complexity + Rarity,data=wild_mo)
vif(co)
plot(x=wild_mo$Sex, y=wild_mo$ExactAge)
# since there is overlap and rather small variation of age within sexes, it makes sense to include them as an interaction

# results: full vs. null model
wild_mo$Success=as.factor(wild_mo$Success)
null=glmer(Success ~ 1+(1|BeggerAge)+(1|Mother/Begger.id)+(0+ExactAge|Mother/Begger.id) + (1|FoodItem)+ (0+Complexity|FoodItem)+(0+Rarity|FoodItem),family=binomial(link="logit"), data=wild_mo, control=contr)
anova(null, interaction, test="Chisq")
# clear effect of the test predictors on the probability of being successful (chi square=19.614, df=6, P=0.003243)

# to get LRT, df, and p values of each of the fixed effects:
drop1(interaction, test="Chisq")

# get estimates and SE
round(summary(interaction)$coefficients,3)

# get sd. for random effects:
summary(interaction)$varcor

# get CIs
confint.merMod(object=interaction)

# plot the significant results: Age, complexity and the interaction between age and sex. 
plot(allEffects(interaction), ylab="Probability of Begging Success", selection=1)
plot(allEffects(interaction), ylab="Probability of Begging Success", selection=3)
plot(allEffects(interaction), ylab="Probability of Begging Success", selection=4)
library(emmeans)
(mylist <- list(sqrt.Age=seq(0,4,by=0.1),Sex=c("female","male")))
emmip(interaction,Sex ~ sqrt.Age, at=mylist,CIs=TRUE)

### End of file ###