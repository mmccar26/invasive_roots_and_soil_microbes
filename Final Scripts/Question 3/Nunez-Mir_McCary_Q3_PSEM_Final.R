###Nunez-Mir and McCary, invasive root traits analysis
##Question 3 analysis

######## Piecewise Structural Equation Modeling#################
##Load relevant R packages
library(corrplot)
library(tidySEM)
library(piecewiseSEM)
library(lme4)
library(car)
library(olsrr)
library(tidyverse)
library(ggplot2)

##Import data
data<-read.csv("question3.csv", header = TRUE)
##data is already scaled

###********************************************************
###Check for variance inflation for main root traits (possible vif concerns)
###********************************************************

#check correlations
# Correlation plot of root traits
c <- cor(data[16:23])
corrplot(c, method="number")
##highest correlation is 0.65, no suggestive multicollinearity violations
##no vif concerns for the main regressions

###******************************************************************
###Piecewise SEM analyses
##******************************************************************

#### chosen model##############
model1 <- psem(lm(traitsPC1 ~ InvFraction, data = data),
               lm(traitsPC2 ~ InvFraction, data = data),
               lm(AvgOfpOxalate ~ traitsPC2 + traitsPC1, data = data),
               lm(AvgOfphH2o ~ traitsPC2 + traitsPC1, data = data),
               lm(AvgOfctonRatio ~ traitsPC2 + traitsPC1 + InvFraction, data = data),
               lm(Ratio ~  AvgOfpOxalate + AvgOfphH2o + AvgOfctonRatio + InvFraction + traitsPC2, data = data),
               AvgOfctonRatio %~~% AvgOfphH2o,
               AvgOfctonRatio %~~% AvgOfpOxalate
)
summary(model1, standardize = "scale")
fisherC(model1)
plot(model1)
rsquared(model1)


#### chosen model############## I know we made a blood oath but this is the chosen model with Ecosystem as random effect. Same outcome
model2 <- psem(lmer(traitsPC1 ~ InvFraction+ (1|Ecosystem), data = data),
               lmer(traitsPC2 ~ InvFraction+ (1|Ecosystem), data = data),
               lmer(AvgOfpOxalate ~ traitsPC2 + traitsPC1+ (1|Ecosystem), data = data),
               lm(AvgOfphH2o ~ traitsPC2 + traitsPC1, data = data),
               lmer(AvgOfctonRatio ~ traitsPC2 + traitsPC1 + InvFraction+ (1|Ecosystem), data = data),
               lm(Ratio ~  AvgOfpOxalate + AvgOfphH2o + AvgOfctonRatio + InvFraction, data = data),
               AvgOfctonRatio %~~% AvgOfphH2o,
               AvgOfctonRatio %~~% AvgOfpOxalate,
               AvgOfphH2o %~~% AvgOfpOxalate
)
summary(model2, standardize = "scale")
fisherC(model2)
plot(model2)




