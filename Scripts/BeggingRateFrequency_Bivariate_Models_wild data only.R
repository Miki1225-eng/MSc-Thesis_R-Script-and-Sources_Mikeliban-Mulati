#For Begging Frequencies and Begging Rates Computations: with wild data set only
#################################################################################################

# set working directory
setwd("C:/Users/Miki/Google Drive/Thesis!!!/Plots & R Script")

# read my full table for the purpose of analyses (full_ana). 
#The dofference between this table and my original full table is that this table has no P culumns or comment or missing values written as NA (just empty instead)!
con<-read.delim("Condensed.txt", header=T, sep="\t",na.strings=c("", "NA"))
str(con) 
half<-read.delim("HalfYearTable.txt", header=T, sep="\t",na.strings=c("", "NA"))
str(half)

# handle the missing data 
half=half[!is.na(half$Age),]
half=half[!is.na(half$BeggingFreq),]


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
library(sjmisc)
library(splines)


############################################################################################################

#discover all the correlations on condensed_food-based_Suggested by Miki
#begging rate/begging freq ~ age/sex/complexity/rarity
#begging rate/begging freq ~ age
bfa<-cor(con$Age,con$BeggingFreq,use="pairwise",method="pearson")
bra<-cor(con$Age,con$BeggingRate,use="pairwise",method="pearson") 
#no correlation in any.
#visulization
ggqqplot(con$Age)
#seems ok, can try Shapiro-Wilk test to test the distribution.
shapiro.test(con$Age)
#age is definitely not normally distributed.
shapiro.test(con$BeggingFreq)
shapiro.test(con$BeggingRate) 
#none of them is normally distributed, meaning that no parametric test is applicable.   
#success rate/begging rate/begging freq ~ gender: Mann-Whitney U-Test &  Kruskal-Wallis test
brsw<-wilcox.exact(con$BeggingRate~con$Sex)
brskph<-kruskalmc(con$BeggingRate~con$Sex)
brsk<-kruskal.test(con$BeggingRate~con$Sex)
#there is a significant gender difference in begging rate! 
#visulization:box plot!
bfsw<-wilcox.exact(con$BeggingFreq~con$Sex)
bfskph<-kruskalmc(con$BeggingFreq~con$Sex)
bfsk<-kruskal.test(con$BeggingFreq~con$Sex) 
#there is a significant gender difference in begging frequency!
#visulization:box plot!
sexbf<-ggboxplot(con,x="Sex",y="BeggingFreq",color = "Sex",palette = "jco")+stat_compare_means()
#begging rate/begging freq ~ complexity
bfc<-cor(con$ProcessingSteps_Complexity,con$BeggingFreq,use="pairwise",method="spearman")
brc<-cor(con$ProcessingSteps_Complexity,con$BeggingRate,use="pairwise",method="spearman") 
#no correlation in either.
#success rate/begging rate/begging freq ~ rarity
bfr<-cor(con$PopFreq_Rarity,con$BeggingFreq,use="pairwise",method="pearson")
brr<-cor(con$PopFreq_Rarity,con$BeggingRate,use="pairwise",method="pearson") 
#no correlation in any.


##############################################################################################################################
#discover all the correlations on halfyearly_0.5yearbased.
#begging rate/begging freq ~ age/gender
#begging rate/begging freq ~ age
hbfa<-cor(half$Age,half$BeggingFreq,use="pairwise",method="pearson")
hbra<-cor(half$Age,half$BeggingRate,use="pairwise",method="pearson") 
#a bit corelation between age and begging frequency.
#visulization
ggqqplot(half$Age)
#seems ok, can try Shapiro-Wilk test to test the distribution.
shapiro.test(half$Age)
#a bit skewed.
shapiro.test(half$BeggingFreq)
shapiro.test(half$BeggingRate)
#none of them is normally distributed, meaning that no parametric test is applicable once again.
ggscatter(half,x = "Age",y = "BeggingFreq", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",col="yellow",
          xlab = "Age",ylab = "Begging Frequency")   
#success rate/begging rate/begging freq ~ gender: Mann-Whitney U-Test &  Kruskal-Wallis test
hbrsw<-wilcox.exact(half$BeggingRate~half$Sex)
hbrskph<-kruskalmc(half$BeggingRate~half$Sex)
hbrsk<-kruskal.test(half$BeggingRate~half$Sex)
#there is no significant gender difference in begging rate. 
hbfsw<-wilcox.exact(half$BeggingFreq~half$Sex)
hbfskph<-kruskalmc(half$BeggingFreq~half$Sex)
hbfsk<-kruskal.test(half$BeggingFreq~half$Sex) 
#there is a significant gender difference in begging frequency!
#visulization:box plot!
hsexbf<-ggboxplot(half,x="Sex",y="BeggingFreq",color ="Sex",palette = "jco")+stat_compare_means()


################################################################################################################################
#In terms of Begging frequency, the full model should look like this: 
#Begging Frequency ~ Age + Sex 
#Here we are only interested in the overall effects of AgeClass and Sex.
str(half)

# First steps (similar to 'standard' linear models):
# Check the distribution of the dependent measure and the covariate:

# Check the balance of the design. How many subjects were tested per condition/age cell and how often?
table(half$Focal,half$Sex) 
table(half$Age, half$Sex, half$BeggingFreq,half$NrBeggingEvents) 

# Model structure.
full_begfreq = lmer(BeggingFreq ~ Age + Sex +(1|Mother/Focal), data = half, REML=F)
drop1(full_begfreq, test="Chisq")
summary(full_begfreq)

# Omnibus test:
null.model = lmer(BeggingFreq ~ (1|Mother/Focal), data = half, REML=F)
anova(null.model,full_begfreq,test="Chisq")

# Sometimes models do not converge. Here you change the number of simulations to, e.g., 500000.
contr=lmerControl(optCtrl=list(maxfun=500000))

# Check residuals and residuals vs. fitted vallues.
hist((residuals(full_begfreq)))
qqplot((residuals(full_begfreq)), rnorm(length((residuals(full_begfreq)))), pch=19, col="grey")
qqline((residuals(full_begfreq)), rnorm(length((residuals(full_begfreq)))))
plot(fitted(full_begfreq), residuals(full_begfreq))
#looks a bit funny... try it again with adjusted follow hours.

# Look at the random effects:
ranef(full_begfreq)

# Check for influential cases.
# Collinearity.
vif(lm(BeggingFreq ~ Age + Sex + NrBeggingEvents, data = half))
#we don't put any  interactions between any two predictors in this model. 

# Confidence intervals for the beta coefficients.
confint.merMod(full_begfreq)

# Effect size simialr to R^2
r.squaredGLMM(full_begfreq)

# -> R2m (Effect of fixed effects)
# -> R2c (Effect of fixed and random effects) which is way larger in the model than R2m. 

# Power analysis

power.res = powerSim(full_begfreq, test=compare(null.model, method = c("lr")), seed = 1, nsim=10)

# Run multiple power analysis simulations (of the simulations). 

power.list2 <- vector('list',10)
power.mean2 = rep(NA, 10)
power.ci2 = matrix(rep(NA, 10*2),ncol=2)

for(a in 1:10){	
power.temp = powerSim(full_begfreq, test=compare(null.model, method = c("lr")), seed = 1, nsim=10)
temp = summary(power.temp)
power.list2[a] <- power.temp
power.mean2[a] <- temp$mean
power.ci2[a,] <- as.numeric(temp[4:5])
rm("power.temp", "temp")
}

#plots: begging frequency
#a) Begging frequency ~ Sex -  as a box plot.
ggboxplot(half,x="Sex",y="BeggingFreq",color ="Sex",palette = "jco")
iv1wild<- ggpredict(full_begfreq,"Sex")
plot(iv1wild)
#b) Begging frequency ~ Age - as a scatter plot with Age on the x axis and begging frequency on the y axis. 
ggplot(half, aes(Age,BeggingFreq)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE)
iv2wild<- ggpredict(full_begfreq,"Age")
plot(iv2wild)

########################################################################################################################
#For Begging rate, using the GLMM approach, the full model would look like this: 
#Begging rate ~ Age + Sex + Rarity + Complexity +(1|Baggar) + (1|Mother) with a Gaussian family distribution. 
#since there are rarity and complexity variables, I'd use the condensed table.
con$Age=as.numeric(con$Age)

# First steps (similar to 'standard' linear models):
# Check the distribution of the dependent measure and the covariate:

# Check the balance of the design. How many subjects were tested per condition/age cell and how often?
table(con$Beggar,con$Sex) 
table(con$Age,con$Sex,con$ProcessingSteps_Complexity,con$PopFreq_Rarity,con$NrBeggingEvents) 
str(con)
# Model structure.
full_begrate = lmer(BeggingRate~ Age + Sex +ProcessingSteps_Complexity+ PopFreq_Rarity +(1|Mother/Beggar), data = con, REML=F)
drop1(full_begrate, test="Chisq")
summary(full_begrate)

# Omnibus test:
null.model = lmer(BeggingRate ~ (1|Mother/Beggar), data = con, REML=F)
anova(null.model,full_begrate,test="Chisq")

# Sometimes models do not converge. Here you change the number of simulations to, e.g., 500000.
contr=lmerControl(optCtrl=list(maxfun=500000))

# Check residuals and residuals vs. fitted vallues.
hist(log(residuals(full_begrate)))
plot((fitted(full_begrate)), (residuals(full_begrate)))
#after transforming, not so problematic. 

# Look at the random effects:
ranef(full_begrate)

# Check for influential cases.
# Collinearity.
vif(lm(BeggingRate~ Age + Sex +ProcessingSteps_Complexity+ PopFreq_Rarity, data = con))
#we don't put any  interactions between any two predictors in this model. 

# Confidence intervals for the beta coefficients.
# The following takes time!
confint.merMod(full_begrate)

# Effect size simialr to R^2
r.squaredGLMM(full_begrate)

# -> R2m (Effect of fixed effects)
# -> R2c (Effect of fixed and random effects) 
#these have the same value. 

#plots: begging rate
#a) Begging Rate ~ Sex -  as a box plot
ggboxplot(con,x="Sex",y="BeggingRate",color ="Sex",palette = "jco")
iv1wild<- ggpredict(full_begrate,"Sex")
plot(iv1wild)
#b) Begging Rate ~ Age - as a scatter plot with Age on the x axis and begging rate on the y axis. 
ggplot(con, aes(Age,BeggingRate)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE)
iv2wild<- ggpredict(full_begrate,"Age")
plot(iv2wild)
#c) Begging Rate ~ Complexity - as a scatter plot with Complexity on the x axis and begging rate on the y axis. 
ggplot(con, aes(ProcessingSteps_Complexity,BeggingRate)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE)
iv3wild<- ggpredict(full_begrate,"ProcessingSteps_Complexity")
plot(iv3wild)
#d) Begging Rate ~ Rarity - as a scatter plot with Rarity on the x axis and begging rate on the y axis. 
ggplot(con, aes(PopFreq_Rarity,BeggingRate)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE)
iv4wild<- ggpredict(full_begrate,"PopFreq_Rarity")
plot(iv4wild)











