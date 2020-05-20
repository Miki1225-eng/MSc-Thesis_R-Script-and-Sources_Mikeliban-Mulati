############################################################################
#Part 4: comparison of wild and zoo data sets

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
library(sjmisc)
library(splines)
library(effsize)
library(BBmisc)
library(moments)

#For the zoo vs wild comparison, I would look at the following models: 
#a) BeggingFreq ~ Site + AgeClass + (1|Beggar)  (Gaussian, data:conzoowild)
#b) Success ~ Site + ExactAge + (1|AgeClass) +(1|Beggar) (Binomial,data:full_ana)

# set working directory
setwd("C:/Users/Miki/Google Drive/Thesis!!!/Plots & R Script")
conzoowild<-read.delim("Condensed_zoowildcompar.txt", header=T, sep="\t",na.strings=c("", "NA"))
str(conzoowild)
full_ana<-read.delim("Full_ExactAge_SimpleFood_WildZoo_AnalysesPurpose_NoPColumns_noNA_Miki.txt", header=T, sep="\t",na.strings=c("", "NA"))
str(full_ana) 
zooage<-read.delim("Subset_zoo_ExactAge.txt", header=T, sep="\t",na.strings=c("", "NA"))
str(zooage)  
zoocon<-read.delim("Zoo condensed.txt", header=T, sep="\t",na.strings=c("", "NA"))
str(zoocon) 

#model a)
# Check the balance of the design. How many subjects were tested per condition/age cell and how often?
table(conzoowild$Beggar,conzoowild$AgeClass) 
table(conzoowild$AgeClass,conzoowild$Site)
conzoowild=conzoowild[!is.na(conzoowild$BeggingFreq),]
conzoowild=conzoowild[!is.na(conzoowild$AgeClass),]
conzoowild$AgeClass=as.numeric(conzoowild$AgeClass)

# Check assumptions for a (multiple) linear regression
# Independence of observations.
# Distribution of 'BeggingFreq'.
hist(log(conzoowild$BeggingFreq))

# -> Looks ok.
# Model structure: BeggingFreq ~ Site + AgeClass.
full_conwz = lm(BeggingFreq~ Site + AgeClass, data = conzoowild)
summary(full_conwz)

# Distribution of model residuals.
hist(log(residuals(full_conwz)))
# -> Looks ok.

# Variance of model residuals.
plot(fitted(full_conwz), residuals(full_conwz), pch=19, cex=1.5)
#-> Looks good overall but there are fewer large fitted values and model fitted appears to be increase with the increasing values of the fitted variable (i.e., 'funnel-pattern'). Check that the model results (test statisics and probability values) are similar with a transformed dependent measure.

# Check for influential cases.
# Collinearity.
vif(full_conwz)
# All vifs (1.024017) smaller than 10 suggesting no issue of collinearity.

# Influential cases.
# Running the model with and without each factor and looking at the change in regression slopes. 
# CAUTION. The function dfbetas returns standardized values whereas the function dfbeta returns unstandardized values.
df.res = dfbetas(full_conwz)
df.res = data.frame(df.res)
max(abs(df.res$Site))
max(abs(df.res$AgeClass))
# 'Site' doesn't look good coz its maximum values is beyond 1.

# To check unstandardized dfbetas a direct comparison of betas with and without the respective data point included. See the following table:
table.dfbeta = data.frame(round(cbind(coefficients(full_conwz), coefficients(full_conwz)+ t(apply(X=dfbeta(full_conwz), MARGIN=2, FUN=range))), 5))
names(table.dfbeta) <- c('Coefficients from the model', 'Smallest coefficient following casewise exclusion', 'Largest coefficient following casewise exclusion')
table.dfbeta
# We compare the variation withn each row of the predictor. Simialr to interpreting plots, we are looking for odd pattern, in this case large variation.

# Running the model with and without each factor and looking at the change in fitted values. 
# The function dffits returns standardized values.
dffit.res = dffits(full_conwz)
dffit.res = data.frame(dffit.res)
max(dffit.res)
# looks good because the maximum values is below 2.

# Cook's distance provides a standardized measure of influence per case.
max(cooks.distance(full_conwz))
# Looks good; it's lower than 1.

# Generate Plot.
# x-axis: AgeClass, y-axis:BeggingFreq, regression line for each category of Site.
# Let's look at the linear equation using the output from summary(full_conwz)
# Y ~ 0.0193644 +0.2377918*X(Site)-0.0002954 *X(AgeClass).

# To enter the model lines for each level of the categorical predcitor  proceed as follows:
# For Site=zoo
# Intercept:  0.0193644 +0.2377918*1-0.0002954 *0 -->0.2571562
# Slope: -0.0002954
x.zoo = seq(min(conzoowild$AgeClass),max(conzoowild$AgeClass), by=0.01)
y.zoo =  0.2571562-0.0002954*x.zoo
plot(conzoowild$AgeClass,conzoowild$BeggingFreq, pch=19, cex=1)
lines(x.zoo, y.zoo,col="blue",lwd=1.5)

#ESs and CIs! Cohen's D 
cohen.d(conzoowild$BeggingFreq,conzoowild$AgeClass)
#d estimate: -2.344574 (large)!
#95 percent confidence interval:
#    lower     upper 
#-2.517038 -2.172109 
#cohen.d(conzoowild$BeggingFreq,conzoowild$Site)
#d estimate: -2.987189 (large)!
#95 percent confidence interval:
#    lower     upper 
#-3.372149 -2.602228 

#plots
iv1<- ggpredict(full_conwz,"Site")
plot(iv1)
iv2<- ggpredict(full_conwz,"AgeClass")
plot(iv2)
#Begging freq ~ AgeClass & Site -  as a box plot colored with Site
ggboxplot(conzoowild,x="AgeClass",y="BeggingFreq",color ="Site",palette = "jco")
ggboxplot(conzoowild,x="Site",y="BeggingFreq",color ="Site",palette = "jco")+stat_compare_means()


#model b)
full_ana=full_ana[!is.na(full_ana$Success),]
full_ana=full_ana[!is.na(full_ana$ExactAge),]
full_ana=full_ana[!is.na(full_ana$Site),]
str(full_ana)
full_gen<- glm(Success ~ Site + ExactAge ,data = full_ana, family = binomial)
# print the mod results without correlations among fixed effects
summary(full_gen)

# Note that both the command structure and the summary-table are quite comparable to what we would get from using the command lm().
# Some important points to keep in mind:
# The fitted values we get from running fitted() are the values for the linear predictor and the fitted values related to the dependent measure are already calculated based on the above equation.

# Distribution of model residuals.
hist(log(residuals(full_gen)))
# hmm...

# Variance of model residuals.
plot(fitted(full_gen), residuals(full_gen), pch=19, cex=1.5)
#-> Looks good overall but there are fewer large fitted values and model fitted appears to be increase with the increasing values of the fitted variable (i.e., 'funnel-pattern'). Check that the model results (test statisics and probability values) are similar with a transformed dependent measure.

# The fitted values are already bound between 0 and 1:
hist(1/(fitted(full_gen)))
#looking fine.

plot(full_ana$ExactAge,fitted(full_gen), pch=15)
#the upper line shows zoo values,the lower one shows the wild. 

# The estimates, i.e., the beta-coefficients, from the summary table can be interpreted in the 
#same way we interpret them from a lm()-output. With increase the advertisment th likelihood 
#of the album selling well increases. In addition, low attraction levels were associated with
#a reduced likelihood that the album would sell well, at least compared to high levels of 
#attractiveness. Here we are only looking at the estimates!

# To calculate p-values, we proceed as follows:
# Calculate an omnibus model comparison test.  
# Note that this step is not always necessary, especially if there are a priori predictions with 
regards to predictor variables. Here we simply want to acknowledge the concept of an omnibus test.

# The null model essentialy predicts the same value for each 'participant'. 
full_gen_nullModel = glm(Success ~ 1, family = binomial, data = full_ana)

# Through comparing the two models we can address the following important question: 
#Do the predictors have a combined incluence on the dependent measure.
anova(full_gen_nullModel,full_gen, test="Chisq")
# -> We see an improve in model fit!

# There is one function that allows us to do both model comparisons at once:
drop1(full_gen, test="Chisq")
#drop1--> to get each of the predictors at once

# Get the the 95% confidence interval for the beta-coefficients:
confint(full_gen)
cbind(coefficients(full_gen), confint(full_gen)) 
#cbind--> to get the the 95% confidence interval

# Generate Plot.
# x-axis: ExactAge, y-axis: Success, model line for each category of Site.
# Let's look at the linear equation using the output from summary.
# Y ~ -0.22296+0.58400*X(Site)+0.11977*X(ExactAge).

# To enter the model lines for each level of the categorical predcitor proceed as follows:
# For Site:zoo.
# Intercept:  -0.22296+0.58400*1+0.11977*0 -> 0.36104
# Slope: 0.11977
x.zoo = seq(min(full_ana$ExactAge),max(full_ana$ExactAge), by=0.1)
y.zoo =  0.36104+ 0.11977*x.zoo
# And now the important addition of adding the model line!
y.zoo = exp(y.zoo)/(1+exp(y.zoo))
plot(full_ana$Success,full_ana$ExactAge, xlab="Success", ylab="ExactAge",pch=19, cex=1)
lines(x.zoo, y.zoo, col="blue", lwd=1.5)

#plots
#a) Success ~ Site 
plot(full_ana$Site,full_ana$Success,xlab="Site", ylab="Success")
#b) Succuss ~ ExactAge 
ggboxplot(full_ana,x="Success",y="ExactAge")+stat_compare_means()
iv1<- ggpredict(full_gen,"Site")
plot(iv1)
iv2<- ggpredict(full_gen,"ExactAge [all]")
plot(iv2)

########################################################################################################################
#For the analyses on the Zoo data only I would do the following 2 main models: 
#c) Success ~ ExactAge + Sex+ Complexity+ Rarity + Desirability (Binomial glm, data:zooage)
#d) BeggingFrequency ~ AgeClass + Complexity + Rarity + Desirability (lm, data: zoocon)

#model c)
str(zooage)
zooage$AgeClass=as.factor(zooage$AgeClass)
zooage=zooage[!is.na(zooage$Success),]
zooage=zooage[!is.na(zooage$ExactAge),]
zooage=zooage[!is.na(zooage$Sex),]
zooage=zooage[!is.na(zooage$ProcessingSteps_Complexity),]
zooage=zooage[!is.na(zooage$PopFreq_Rarity_Onetofive),]
zooage=zooage[!is.na(zooage$Desirability_OnetoSeven),]
str(zooage)
full_zooage<- glm(Success ~ ExactAge + Sex + ProcessingSteps_Complexity+ PopFreq_Rarity_Onetofive + Desirability_OnetoSeven,data = zooage, family = binomial)
# print the mod results without correlations among fixed effects
summary(full_zooage)

# Note that both the command structure and the summary-table are quite comparable to what we would get from using the command lm().
# Some important points to keep in mind:
# The fitted values we get from running fitted() are the values for the linear predictor and the fitted values related to the dependent measure are already calculated based on the above equation.

# Distribution of model residuals.
hist(sqrt(residuals(full_zooage)))
# looks better than before transformed.

# Variance of model residuals.
plot(fitted(full_zooage), residuals(full_zooage), pch=19, cex=1.5)
#-> Looks good overall but there are fewer large fitted values and model fitted appears to be increase with the increasing values of the fitted variable (i.e., 'funnel-pattern'). Check that the model results (test statisics and probability values) are similar with a transformed dependent measure.

# The fitted values are already bound between 0 and 1:
hist(sqrt(fitted(full_zooage)))
#looking good.

# The estimates, i.e., the beta-coefficients, from the summary table can be interpreted in the 
#same way we interpret them from a lm()-output. With increase the advertisment th likelihood 
#of the album selling well increases. In addition, low attraction levels were associated with
#a reduced likelihood that the album would sell well, at least compared to high levels of 
#attractiveness. Here we are only looking at the estimates!

# To calculate p-values, we proceed as follows:
# Calculate an omnibus model comparison test.  
# Note that this step is not always necessary, especially if there are a priori predictions with 
regards to predictor variables. Here we simply want to acknowledge the concept of an omnibus test.

# The null model essentialy predicts the same value for each 'participant'. 
full_zooage_nullModel = glm(Success ~ 1, family = binomial, data = zooage)

# Through comparing the two models we can address the following important question: 
#Do the predictors have a combined incluence on the dependent measure.
anova(full_zooage_nullModel,full_zooage, test="Chisq")
# no effect overall. 

# There is one function that allows us to do both model comparisons at once:
drop1(full_zooage, test="Chisq")
#drop1--> to get each of the predictors at once

# Get the the 95% confidence interval for the beta-coefficients:
confint(full_zooage)
cbind(coefficients(full_zooage), confint(full_zooage)) 
#cbind--> to get the the 95% confidence interval

# Important information about the model:
# Linear model equation: Y ~ 1.57761 -0.81196*X(ExactAge)+0.50063*X(Sex)-0.17393*X(Complexity)+0.33002*X(Rarity)-0.05816*X(Desirability)
# Note that Y refers to the linear predictor, not the actual dependent measure!

# Power analysis
sample.size = 132
intercept = 1.57761
slope.Sex = 0.50063
slope.Exactage=-0.81196
slope.Complexity=-0.17393
slope.Rarity=0.33002
slope.Desirability=-0.05816
nr.simulations = 100
alpha = 0.05

p.vals=rep(x=0, time= nr.simulations)

set.seed(11)

for(a in 1: nr.simulations){

	# Sample factor
	Sex.sample = sample(x=zooage$Sex,size= sample.size, replace=T)
	levels(Sex.sample) <- c(0,1)
	Sex.sample = as.numeric(as.character(Sex.sample))

	# Generate linear predictors, not the dependent measure!  
	lp.sample = intercept + slope.Sex*Sex.sample+slope.Complexity+slope.Rarity+slope.Desirability+slope.Exactage
	dm.sample=rbinom(n=sample.size, size=1, prob = exp(lp.sample)/(1+exp(lp.sample)))
	
	# Run the model.
	model.sample = glm(dm.sample ~ Sex.sample+lp.sample, family="binomial")
	model.sample = as.data.frame(drop1(model.sample, test="Chisq"))
	
	# Store '1' for all predictors if p < .05 
	if(is.na(model.sample["Sex.sample", "Pr(>Chi)"] <= alpha && model.sample["lp.sample", "Pr(>|t|)"] <= alpha)){p.vals[a] <- 1}

	rm("Sex.sample", "dm.sample", "lp.sample", "model.sample")
}

# Look at the result.
mean(p.vals)
#1
#now all predictors show p < .05, but why there is no significance in drop 1 function?


#plots
#predicted probabilities 
iv1<- ggpredict(full_zooage,"ExactAge [all]")
plot(iv1)
iv2<- ggpredict(full_zooage,"Sex")
plot(iv2)
iv3<- ggpredict(full_zooage,"ProcessingSteps_Complexity")
plot(iv3)
iv4<- ggpredict(full_zooage,"PopFreq_Rarity_Onetofive")
plot(iv4)
iv5<- ggpredict(full_zooage,"Desirability_OnetoSeven")
plot(iv5)


#model d): BeggingFrequency ~ AgeClass + Sex+Complexity + Rarity + Desirability (lm, data: zoocon)
# Check the balance of the design. How many subjects were tested per condition/age cell and how often?
table(zoocon$Beggar,zoocon$AgeClass) 
table(zoocon$AgeClass,zoocon$Sex,zoocon$Complexity,zoocon$Rarity,zoocon$Desirability)
zoocon=zoocon[!is.na(zoocon$BeggingFreq),]
zoocon=zoocon[!is.na(zoocon$AgeClass),]
zoocon=zoocon[!is.na(zoocon$Complexity),]
zoocon=zoocon[!is.na(zoocon$Rarity),]
zoocon=zoocon[!is.na(zoocon$Desirability),]
str(zoocon)
# Check assumptions for a (multiple) linear regression
# Independence of observations.
# Distribution of 'BeggingFreq'.
hist(log(zoocon$BeggingFreq))
# not sure about the distribution. Run a kurtosis skewness test (The kurtosis of any univariate normal distribution is 3)
kurtosis(log(zoocon$BeggingFreq))
#way better after transforming

full_zoocon = lm(BeggingFreq~AgeClass+Sex+Complexity+Rarity+Desirability,data = zoocon)
summary(full_zoocon)

# Omnibus test:
null.model = lm(BeggingFreq~1, data = zoocon)
anova(null.model,full_zoocon,test="Chisq")
#generally no effect.

# Distribution of model residuals.
hist((residuals(full_zoocon)))
# -> Looks ok.

# Variance of model residuals.
plot(fitted(full_zoocon), residuals(full_zoocon), pch=19, cex=1.5)
#-> Looks good overall but there are fewer smaller fitted values and model fitted appears to be increase with the increasing values of the fitted variable (i.e., 'funnel-pattern'). Check that the model results (test statisics and probability values) are similar with a transformed dependent measure.

# Check for influential cases.
# Collinearity.
vif(full_zoocon)
# All vifs (1.024017) smaller than 10 suggesting no issue of collinearity.

#ESs and CIs! Cohen's D 
cohen.d(zoocon$BeggingFreq,zoocon$AgeClass)
#d estimate: -4.100388 (large)
95 percent confidence interval:
    lower     upper 
-4.905446 -3.295330  
cohen.d(zoocon$BeggingFreq,zoocon$Sex)
#d estimate: -0.2672307 (small)
95 percent confidence interval:
     lower      upper 
-0.9318214  0.3973600 
cohen.d(zoocon$BeggingFreq,zoocon$Complexity)
#d estimate: -0.7762471 (medium)
95 percent confidence interval:
     lower      upper 
-1.2502706 -0.3022236 
cohen.d(zoocon$BeggingFreq,zoocon$Rarity)
#d estimate: -3.576563 (large)
95 percent confidence interval:
    lower     upper 
-4.313503 -2.839623
cohen.d(zoocon$BeggingFreq,zoocon$Desirability)
#d estimate: -4.466015 (large)
95 percent confidence interval:
    lower     upper 
-5.320374 -3.611656 

# Important information about the model:
# Linear model equation: Y ~  -0.005864 +0.098674*X(AgeClass) -0.032203*X(Sex)-0.073583 *X(Complexity)+ 0.054946 *X(Rarity)-0.023933*X(Desirability)
# Note that Y refers to the linear predictor, not the actual dependent measure!

# Power analysis
sample.size = 38
intercept = -0.005864
slope.Sex = -0.032203
slope.AgeClass= 0.098674
slope.Complexity= -0.073583
slope.Rarity= 0.054946
slope.Desirability= -0.023933
nr.simulations = 100
alpha = 0.05

p.vals=rep(x=0, time= nr.simulations)

set.seed(11)

for(a in 1: nr.simulations){

	# Sample factor
	Sex.sample = sample(x=zoocon$Sex,size= sample.size, replace=T)
	levels(Sex.sample) <- c(0,1)
	Sex.sample = as.numeric(as.character(Sex.sample))

	# Generate linear predictors, not the dependent measure!  
	lp.sample = intercept + slope.Sex*Sex.sample+slope.Complexity+slope.Rarity+slope.Desirability+slope.AgeClass
	dm.sample=rbinom(n=sample.size, size=1, prob = exp(lp.sample)/(1+exp(lp.sample)))
	
	# Run the model.
	model.sample = glm(dm.sample ~ Sex.sample+lp.sample, family="binomial")
	model.sample = as.data.frame(drop1(model.sample, test="Chisq"))
	
	# Store '1' for all predictors if p < .05 
	if(is.na(model.sample["Sex.sample", "Pr(>Chi)"] <= alpha && model.sample["lp.sample", "Pr(>|t|)"] <= alpha)){p.vals[a] <- 1}

	rm("Sex.sample", "dm.sample", "lp.sample", "model.sample")
}

# Look at the result.
mean(p.vals)
#1
#now all predictors show p < .05, but why there is no significance in drop 1 function?

#plots
#Begging freq ~ Complexity_raw data.
ggplot(zoocon, aes(Complexity,BeggingFreq)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE)
iv1<- ggpredict(full_zoocon,"Complexity")
plot(iv1)
#Begging freq ~ Rarity
ggplot(zoocon, aes(Rarity,BeggingFreq)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE)
iv2<- ggpredict(full_zoocon,"Rarity")
plot(iv2)
#Begging freq ~ Desirability
ggplot(zoocon, aes(Desirability,BeggingFreq)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE)
iv3<- ggpredict(full_zoocon,"Desirability")
plot(iv3)
#Begging freq ~ AgeClass -  as a box plot colored with Sex
ggplot(zoocon, aes(AgeClass,BeggingFreq)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE)
iv4<- ggpredict(full_zoocon,"AgeClass [all]")
plot(iv4)
#Begging freq ~ Sex
ggboxplot(zoocon,x="Sex",y="BeggingFreq",color ="Sex",palette = "jco")
iv5<- ggpredict(full_zoocon,"Sex")
plot(iv5)















