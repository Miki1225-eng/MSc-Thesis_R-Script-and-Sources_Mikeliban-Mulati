#The Role of the Mother in Feeding Skill Acquisition in Immature Sumatran Orangutans_RScript_Miki_31.03
#Supervisors: Dr. Caroline Schuppli; Dr. Robert Hepach 
#########################################################################################################

# set working directory
setwd("C:/Users/Miki/Google Drive/Thesis!!!/Plots & R Script")

#######################################################################################################
#begging success models

# read my full table. 
wild<-read.delim("Subset_wild_ExactAge.txt", header=T, sep="\t", na.strings=c("", "NA"))
str(wild)

# handle the missing data [only for running the model!]
wild_mo<-subset(wild,RelationtoB=="mother")
str(wild_mo)
wild_mo=wild_mo[!is.na(wild_mo$Success),]
wild_mo=wild_mo[!is.na(wild_mo$ExactAge),]
wild_mo=wild_mo[!is.na(wild_mo$FoodItem),]
wild_mo=wild_mo[!is.na(wild_mo$Mother),]
wild_mo=wild_mo[!is.na(wild_mo$Complexity),]
wild_mo=wild_mo[!is.na(wild_mo$Rarity),]
wild_mo=wild_mo[!is.na(wild_mo$Begger.id),]
str(wild_mo)
#########################################################################
#all the packages I will need
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

ggpairs(wild[,c("ExactAge","FollowHours","Complexity","Rarity")])

#################Following the lecture on GLMM (binomial) slides by Roger Mundry (February 5, 2020) and source codes as part of course materials##############
#random effects in the model: focal ID ("Beggar"), "Mother", & food item ID "Food item" 
#preparatory steps
#are the data balanced in each of those?
table(wild$Begger.id)
table(wild$Mother)
table(wild$FoodItem)
#all not really balanced..
#check the frequency of occurrence of success per chimp ID and food item
table(wild$Success,wild$Begger.id)
table(wild$Success,wild$FoodItem)
#check the distributions of the covariates: ExactAge, AgeClass, Complexity, Rarity 
hist(wild$ExactAge)
hist(sqrt(wild$ExactAge))
hist(wild$AgeClass)
hist(wild$Complexity)
hist(wild$Rarity)
#only the complexity is okay ish.
#inspect cross-tabulation of food items and begger id: could there be a special taste of an individual?
table(wild$FoodItem,wild$Begger.id)
table(wild$Begger.id,wild$ExactAge)
table(wild$Begger.id,wild$Mother)
#they are both partly crossed
###################################################################################
#fit the model
contr=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000))
wild_mo$BeggerAge=as.factor(paste(wild_mo$Begger.id,round(wild_mo$ExactAge)))
hist(wild_mo$ExactAge)
describe(wild_mo$ExactAge)
hist(sqrt(wild_mo$ExactAge))
describe(sqrt(wild_mo$ExactAge))
#nice, keep the 'sqrt.Age'
hist(wild_mo$Complexity)
describe(wild_mo$Complexity)
#keep it as it is.
hist(wild_mo$Rarity)
describe(wild_mo$Rarity)
hist(sqrt(wild_mo$Rarity))
describe(sqrt(wild_mo$Rarity))
#not super ideal but much better! Keep it as 'sqrt.Rarity'!
wild_mo$sqrt.Rarity=sqrt(wild_mo$Rarity)
wild_mo$sqrt.Age=sqrt(wild_mo$ExactAge)
interaction<- glmer(Success ~ sqrt.Age + Sex + sqrt.Age * Sex + Complexity + sqrt.Rarity+I(sqrt.Age^2) + 
           (1|BeggerAge) + (1|Mother/Begger.id) + 
           (0 + sqrt.Age|Mother/Begger.id) +(1|FoodItem) + 
           (0 + Complexity | FoodItem)+(0 + sqrt.Rarity| FoodItem),family=binomial(link="logit"),
           data=wild_mo,control=contr)
#to get LRT, df, and p results:
drop1(interaction, test="Chisq")
#to get estimates and SE
round(summary(interaction)$coefficients,3)
#to get sd. for random effects:
summary(interaction)$varcor
#to get CIs
confint.merMod(object=interaction)
#it just didn't work..
#plot the significant results: Age, complexity and the interaction between age and sex. 
plot(allEffects(interaction), ylab="Probability of Begging Success", selection=1)
plot(allEffects(interaction), ylab="Probability of Begging Success", selection=3)
plot(allEffects(interaction), ylab="Probability of Begging Success", selection=4)
library(emmeans)
(mylist <- list(sqrt.Age=seq(0,4,by=0.1),Sex=c("female","male")))
emmip(interaction,Sex ~ sqrt.Age, at=mylist,CIs=TRUE)
###########################more info about the wild data below######################################
#singular fit issue
summary(interaction)$varcor 
ll.old=logLik(interaction)
round(ll.old, 3); logLik(interaction)
#exclude random slopes for ExactAge
interaction.s<- glmer(Success ~ ExactAge + Sex + ExactAge*Sex+Complexity + Rarity + I(ExactAge^2)+(1|BeggerAge)+ (1|FoodItem)+ (0+Complexity|FoodItem)+(0+Rarity|FoodItem),family=binomial(link="logit"), data=wild, control=contr)
#well, it worked! Let's see if they are actually very different
summary(interaction.s)$varcor  
logLik(interaction); logLik(interaction.s)
round(summary(interaction)$coefficients, 3)
round(summary(interaction.s)$coefficients, 3)
#well, they are exactly the same so just keep the original since it shows more info maybe?
ranef.diagn.plot(interaction)
#those distributions look really fine!
#collinearity
wild_mo$Success=as.numeric(wild_mo$Success)
co=lm(Success ~ sqrt.Age + Sex +Complexity + sqrt.Rarity , data=wild_mo)
vif(co)
plot(x=wild$Sex, y=wild$ExactAge)
#since there is overlap and rather small variation of age within sexes, it makes sense to include them as an interaction
#results: full vs. null model
wild_mo$Success=as.factor(wild_mo$Success)
null=glmer(Success ~ 1+(1|BeggerAge)+(1|Mother/Begger.id)+(0+ExactAge|Mother/Begger.id) + (1|FoodItem)+ (0+Complexity|FoodItem)+(0+Rarity|FoodItem),family=binomial(link="logit"), data=wild, control=contr)
anova(null, interaction, test="Chisq")
#clear effect of the test predictors on the probability of being successful (?2=19.614, df =6, P=0.003243)
###################################################################################################################
#Results
#individual predictors
round(summary(interaction)$coefficients,3)
tests=as.data.frame(drop1(interaction,test="Chisq"))
round(tests, 3)
#intercepts
exp(fixef(interaction)["(Intercept)"])/(1+exp(fixef(interaction)["(Intercept)"]))
#random effects
summary(interaction)$varcor
#get confidence intervals (CIs) and rough estimates using the SEs
source("boot_glmm.r")
boot.res=boot.glmm.pred(model.res=interaction,excl.warnings=T,nboots=1000, para=T)
round(boot.res$ci.estimates, 3)
#didn't work...
confint(interaction)
cbind(coefficients(interaction),confint(interaction))
#didn't work either. 
############################################################################################################################################################
#full model
full_ana<-read.delim("Full_ExactAge_SimpleFood_WildZoo_AnalysesPurpose_NoPColumns_noNA_Miki.txt", header=T, sep="\t",na.strings=c("", "NA"))
str(full_ana) 
full_mo<-subset(full_ana,RelationtoB=="mother")
str(full_mo)
full_mo=full_mo[!is.na(full_mo$Success),]
full_mo=full_mo[!is.na(full_mo$ExactAge),]
full_mo=full_mo[!is.na(full_mo$FoodItem),]
full_mo=full_mo[!is.na(full_mo$Mother),]
full_mo=full_mo[!is.na(full_mo$Complexity),]
full_mo=full_mo[!is.na(full_mo$Begger.id),]
str(full_mo)
contr=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000))
full_mo$BeggerAge=as.factor(paste(full_mo$Begger.id,round(full_mo$ExactAge)))
describe(full_mo$Complexity)
#looking good with the distribution and everything else
hist(sqrt(full_mo$ExactAge))
full_mo$sqrt.Age=sqrt(full_mo$ExactAge)
describe(full_mo$sqrt.Age)
#much better with the distribution, keep it as 'sqrt.Age'
#z-tramsform all the variates 
full_mo$z.Complexity=as.vector(scale(full_mo$Complexity))
full_mo$z.Age=as.vector(scale(full_mo$ExactAge))
full<-glmer(Success ~ z.Age + Site + Sex + z.Age * Sex + z.Complexity + I(z.Age^2) + 
           I(z.Age^2) * Sex + (1 | BeggerAge) + (1 | Mother/Begger.id) + 
           (0 + z.Age | Mother/Begger.id) + (0 + I(z.Age^2)| Mother/Begger.id)+(1 | FoodItem) + 
           (0 + z.Complexity | FoodItem),family=binomial(link="logit"),data=full_mo,control=contr)
full.x<-glmer(Success ~ sqrt.Age + Site + Sex + sqrt.Age * Sex + Complexity + I(sqrt.Age^2) + 
             (1 | BeggerAge) + (1 | Mother/Begger.id) + (0 + sqrt.Age | Mother/Begger.id)+
             (1 | FoodItem) + (0 + Complexity | FoodItem),family=binomial(link="logit"),data=full_mo,control=contr)
#to get LRT, df, and p results:
drop1(full, test="Chisq")
drop1(full.x, test="Chisq")
#to get estimates and SE
round(summary(full.x)$coefficients,3)
#to get sd. for random effects:
summary(full.x)$varcor
#Age and Site are significant!
plot(allEffects(full.x), ylab="Probability of Begging Success", selection=1) #the site variable
plot(allEffects(full.x), ylab="Probability of Begging Success", selection=3) #the age variable
#overall, full.x seems better to me..
################################################################################################################################################################
#zoo model
zooage<-read.delim("Subset_zoo_ExactAge.txt", header=T, sep="\t",na.strings=c("", "NA"))
str(zooage) 
zoo_mo<-subset(zooage,RelationtoB=="mother")
str(zoo_mo)
zoo_mo=zoo_mo[!is.na(zoo_mo$Success),]
zoo_mo=zoo_mo[!is.na(zoo_mo$ExactAge),]
zoo_mo=zoo_mo[!is.na(zoo_mo$Sex),]
zoo_mo=zoo_mo[!is.na(zoo_mo$Complexity),]
zoo_mo=zoo_mo[!is.na(zoo_mo$Rarity),]
zoo_mo=zoo_mo[!is.na(zoo_mo$FoodItem),]
zoo_mo=zoo_mo[!is.na(zoo_mo$Desirability),]
str(zoo_mo)
hist(zoo_mo$ExactAge)
describe(zoo_mo$ExactAge)
hist(log(zoo_mo$ExactAge))
describe(log(zoo_mo$ExactAge))
#nice, just keep eaxct age as it is
hist(zoo_mo$Complexity)
describe(sqrt(zoo_mo$Complexity))
#keep it as sqrt.Complexity
hist(zoo_mo$Rarity)
describe(zoo_mo$Rarity)
hist(sqrt(wild_mo$Rarity))
describe(log(zoo_mo$Rarity))
hist(zoo_mo$Desirability)
describe(zoo_mo$Desirability)
#keep the rarity and desirability as they are
zoo_mo$sqrt.Complexity=sqrt(zoo_mo$Complexity)
zoo<- glmer(Success ~ ExactAge+Sex+Complexity+Rarity+Desirability+(1|FoodItem)+
            (0+Complexity|FoodItem)+(0+Rarity|FoodItem)+(0+Desirability|FoodItem),
            family=binomial,data=zoo_mo,control=contr)
#to get LRT, df, and p results:
drop1(zoo,test="Chisq")
#to get estimates and SE
round(summary(zoo)$coefficients,3)
#to get sd. for random effects:
summary(zoo)$varcor
#no significant result.

################################assisting section###############################################
#To see if age has a linear or quadratic effect on success you can compare two models with the anova function:1. Success ~ Age + (1|Individual) and 2. Success~ Age + I(Age^2) + (1|Individual) -> if Nr 2. is preferred you can include age like this in the full model. 
is.ExactAge.linear.v1<- glmer(Success ~ ExactAge + (1|Mother/Begger.id),family=binomial(link="logit"), data=wild_mo, control=contr)
is.ExactAge.linear.v2<- glmer(Success ~ ExactAge + I(ExactAge^2)+(1|Mother/Begger.id),family=binomial(link="logit"), data=wild_mo, control=contr)
anova(is.ExactAge.linear.v1,is.ExactAge.linear.v2,test="Chisq")
#is.ExactAge.linear.v2 is preferred! 














