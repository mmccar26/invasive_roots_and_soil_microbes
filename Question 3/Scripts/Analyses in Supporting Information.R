
library(phylolm)
library(ggpubr)
library(adephylo)
library(ape)
library(phylobase)
library(lme4)
library(lmerTest)
library(corrplot)
library(MASS)

########### Sensitivity Analyses (using data with lower and upper bounds of imputated estimates) ############

######### Lower bound data ##########

data = read.csv("Traits comparisons_Data lower 95.12.21.23.csv",header=T,row.names = 1)
data$nativeStatusCode <- factor(data$nativeStatusCode, levels = c("N","I"))

## Checking normality

hist(data$rootCcontent)
hist(data$rootNcontent)#
hist(data$rootCNratio)#
hist(data$mycorrhizaIntensity)#
hist(data$rootDiameter)#
hist(data$SRL)#
hist(data$RMF)#
hist(data$rootingDepth)#
hist(data$rootTissueDensity)#

## Transformations

b <- boxcox(lm(data$rootNcontent ~ 1))
lambda <- b$x[which.max(b$y)]
data$rootN_t = (data$rootNcontent ^ lambda - 1) / lambda
hist(data$rootN_t)

b <- boxcox(lm(data$rootCNratio ~ 1))
lambda <- b$x[which.max(b$y)]
data$rootCN_t = (data$rootCNratio ^ lambda - 1) / lambda
hist(data$rootCN_t)

data$rootMI_t = log(data$mycorrhizaIntensity+1)
hist(data$rootMI_t)

b <- boxcox(lm(data$rootDiameter ~ 1))
lambda <- b$x[which.max(b$y)]
data$diam_t = (data$rootDiameter ^ lambda - 1) / lambda
hist(data$diam_t)

b <- boxcox(lm(data$SRL ~ 1))
lambda <- b$x[which.max(b$y)]
data$srl_t = (data$SRL ^ lambda - 1) / lambda
hist(data$srl_t)

b <- boxcox(lm(data$RMF ~ 1))
lambda <- b$x[which.max(b$y)]
data$rmf_t = (data$RMF ^ lambda - 1) / lambda
hist(data$rmf_t)

b <- boxcox(lm(data$rootingDepth ~ 1))
lambda <- b$x[which.max(b$y)]
data$rd_t = (data$rootingDepth ^ lambda - 1) / lambda
hist(data$rd_t)

b <- boxcox(lm(data$rootTissueDensity ~ 1))
lambda <- b$x[which.max(b$y)]
data$rtd_t = (data$rootTissueDensity ^ lambda - 1) / lambda
hist(data$rtd_t)


# Root C

rootclm <- lm(rootCcontent~nativeStatusCode,data=data)
summary(rootclm)

rootcphylolm <- phylolm(rootCcontent~nativeStatusCode,data=data,phy=myTree)
summary(rootcphylolm)

plot(data$rootCcontent,rootcphylolm$residuals)

# root N

rootNlm <- lm(rootNcontent~nativeStatusCode,data=data)
summary(rootNlm)

rootNphylolm <- phylolm(rootN_t~nativeStatusCode,data=data,phy=myTree)
summary(rootNphylolm)

plot(data$rootN_t,rootNphylolm$residuals)#

# root CN

rootCNlm <- lm(rootCNratio~nativeStatusCode,data=data)
summary(rootCNlm)

rootCNphylolm <- phylolm(rootCN_t~nativeStatusCode,data=data,phy=myTree)
summary(rootCNphylolm)

plot(data$rootCN_t,rootCNphylolm$residuals)#


# mycorrhizae colonization Intensity**

MIlm <- lm(mycorrhizaIntensity ~ nativeStatusCode,data=data)
summary(MIlm)

rootMIphylolm <- phylolm(rootMI_t~nativeStatusCode,data=data,phy=myTree)
summary(rootMIphylolm)

plot(data$rootMI_t,rootMIphylolm$residuals)


# Root Diameter*

diamlm <- lm(rootDiameter ~ nativeStatusCode,data=data)
summary(diamlm)

diamphylolm <- phylolm(diam_t~nativeStatusCode,data=data,phy=myTree)
summary(diamphylolm)

plot(data$diam_t,diamphylolm$residuals)

# SRL***

srllm <- lm(SRL ~ nativeStatusCode,data=data)
summary(srllm)

srlphylolm <- phylolm(srl_t~nativeStatusCode,data=data,phy=myTree)
summary(srlphylolm)

plot(data$srl_t,srlphylolm$residuals)


# RMF**

rmflm <- lm(RMF ~ nativeStatusCode,data=data)
summary(rmflm)

rmfphylolm <- phylolm(rmf_t~nativeStatusCode,data=data,phy=myTree)
summary(rmfphylolm)

plot(data$rmf_t,rmfphylolm$residuals)


# Rooting Depth

rdlm <- lm(rootingDepth~ nativeStatusCode,data=data)
summary(rdlm)

rdphylolm <- phylolm(rd_t~nativeStatusCode,data=data,phy=myTree)
summary(rdphylolm)

plot(data$rd_t,rdphylolm$residuals)


# Root Tissue Density

rtdlm <- lm(rootTissueDensity ~ nativeStatusCode,data=data)
summary(rtdlm)

rtdphylolm <- phylolm(rtd_t ~ nativeStatusCode,data=data,phy=myTree)
summary(rtdphylolm)

plot(data$rtd_t,rtdphylolm$residuals)


## Plot results

f1=ggboxplot(data, x = "nativeStatusCode", y = "rootCcontent",
             color = "nativeStatusCode", palette = "jco")

f2=ggboxplot(data, x = "nativeStatusCode", y = "rootNcontent",
             color = "nativeStatusCode", palette = "jco")

f3=ggboxplot(data, x = "nativeStatusCode", y = "rootCN_t",
             color = "nativeStatusCode", palette = "jco")

f4=ggboxplot(data, x = "nativeStatusCode", y = "rootMI_t",
             color = "nativeStatusCode", palette = "jco")

f5=ggboxplot(data, x = "nativeStatusCode", y = "diam_t",
             color = "nativeStatusCode", palette = "jco")

f6=ggboxplot(data, x = "nativeStatusCode", y = "srl_t",
             color = "nativeStatusCode", palette = "jco")

f7=ggboxplot(data, x = "nativeStatusCode", y = "rmf_t",
             color = "nativeStatusCode", palette = "jco")

f8=ggboxplot(data, x = "nativeStatusCode", y = "rd_t",
             color = "nativeStatusCode", palette = "jco")

f9=ggboxplot(data, x = "nativeStatusCode", y = "rtd_t",
             color = "nativeStatusCode", palette = "jco")

ggarrange(f1,f2,f3,f4,f5,f6,f7,f8,f9,common.legend = TRUE, legend = "bottom")


######### Upper bound data ###########

data = read.csv("Traits comparisons_Data upper 95.12.21.23.csv",header=T,stringsAsFactors = T, row.names = 1)
data$nativeStatusCode <- factor(data$nativeStatusCode, levels = c("N","I"))

## Checking normality

hist(data$rootCcontent)#
hist(data$rootNcontent)#
hist(data$rootCNratio)#
hist(data$mycorrhizaIntensity)#
hist(data$rootDiameter)#
hist(data$SRL)#
hist(data$RMF)
hist(data$rootingDepth)#
hist(data$rootTissueDensity)#

## Transformations

b <- boxcox(lm(data$rootCcontent ~ 1))
lambda <- b$x[which.max(b$y)]
data$rootC_t = (data$rootCcontent ^ lambda - 1) / lambda
hist(data$rootC_t)

b <- boxcox(lm(data$rootNcontent ~ 1))
lambda <- b$x[which.max(b$y)]
data$rootN_t = (data$rootNcontent ^ lambda - 1) / lambda
hist(data$rootN_t)

b <- boxcox(lm(data$rootCNratio ~ 1))
lambda <- b$x[which.max(b$y)]
data$rootCN_t = (data$rootCNratio ^ lambda - 1) / lambda
hist(data$rootCN_t)

data$rootMI_t = log(data$mycorrhizaIntensity+1)
hist(data$rootMI_t)

b <- boxcox(lm(data$rootDiameter ~ 1))
lambda <- b$x[which.max(b$y)]
data$diam_t = (data$rootDiameter ^ lambda - 1) / lambda
hist(data$diam_t)

b <- boxcox(lm(data$SRL ~ 1))
lambda <- b$x[which.max(b$y)]
data$srl_t = (data$SRL ^ lambda - 1) / lambda
hist(data$srl_t)

b <- boxcox(lm(data$RMF ~ 1))
lambda <- b$x[which.max(b$y)]
data$rmf_t = (data$RMF ^ lambda - 1) / lambda
hist(data$rmf_t)

b <- boxcox(lm(data$rootingDepth ~ 1))
lambda <- b$x[which.max(b$y)]
data$rd_t = (data$rootingDepth ^ lambda - 1) / lambda
hist(data$rd_t)

b <- boxcox(lm(data$rootTissueDensity ~ 1))
lambda <- b$x[which.max(b$y)]
data$rtd_t = (data$rootTissueDensity ^ lambda - 1) / lambda
hist(data$rtd_t)

# Root C

rootclm <- lm(rootCcontent~nativeStatusCode,data=data)
summary(rootclm)

rootcphylolm <- phylolm(rootC_t~nativeStatusCode,data=data,phy=myTree)
summary(rootcphylolm)

plot(data$rootCcontent,rootcphylolm$residuals)

# root N

rootNlm <- lm(rootNcontent~nativeStatusCode,data=data)
summary(rootNlm)

rootNphylolm <- phylolm(rootN_t~nativeStatusCode,data=data,phy=myTree)
summary(rootNphylolm)

plot(data$rootN_t,rootNphylolm$residuals)#

# root CN*

rootCNlm <- lm(rootCNratio~nativeStatusCode,data=data)
summary(rootCNlm)

rootCNphylolm <- phylolm(rootCN_t~nativeStatusCode,data=data,phy=myTree)
summary(rootCNphylolm)

plot(data$rootCN_t,rootCNphylolm$residuals)#


# mycorrhizae colonization Intensity

MIlm <- lm(mycorrhizaIntensity ~ nativeStatusCode,data=data)
summary(MIlm)

rootMIphylolm <- phylolm(rootMI_t~nativeStatusCode,data=data,phy=myTree)
summary(rootMIphylolm)

plot(data$rootMI_t,rootMIphylolm$residuals)


# Root Diameter*

diamlm <- lm(rootDiameter ~ nativeStatusCode,data=data)
summary(diamlm)

diamphylolm <- phylolm(diam_t~nativeStatusCode,data=data,phy=myTree)
summary(diamphylolm)

plot(data$diam_t,diamphylolm$residuals)

# SRL***

srllm <- lm(SRL ~ nativeStatusCode,data=data)
summary(srllm)

srlphylolm <- phylolm(srl_t~nativeStatusCode,data=data,phy=myTree)
summary(srlphylolm)

plot(data$srl_t,srlphylolm$residuals)


# RMF***

rmflm <- lm(RMF ~ nativeStatusCode,data=data)
summary(rmflm)

rmfphylolm <- phylolm(RMF~nativeStatusCode,data=data,phy=myTree)
summary(rmfphylolm)

plot(data$rmf_t,rmfphylolm$residuals)


# Rooting Depth

rdlm <- lm(rootingDepth~ nativeStatusCode,data=data)
summary(rdlm)

rdphylolm <- phylolm(rd_t~nativeStatusCode,data=data,phy=myTree)
summary(rdphylolm)

plot(data$rd_t,rdphylolm$residuals)


# Root Tissue Density

rtdlm <- lm(rootTissueDensity ~ nativeStatusCode,data=data)
summary(rtdlm)

rtdphylolm <- phylolm(rtd_t ~ nativeStatusCode,data=data,phy=myTree)
summary(rtdphylolm)

plot(data$rtd_t,rtdphylolm$residuals)


## Plot results

f1=ggboxplot(data, x = "nativeStatusCode", y = "rootCcontent",
             color = "nativeStatusCode", palette = "jco")

f2=ggboxplot(data, x = "nativeStatusCode", y = "rootNcontent",
             color = "nativeStatusCode", palette = "jco")

f3=ggboxplot(data, x = "nativeStatusCode", y = "rootCN_t",
             color = "nativeStatusCode", palette = "jco")

f4=ggboxplot(data, x = "nativeStatusCode", y = "rootMI_t",
             color = "nativeStatusCode", palette = "jco")

f5=ggboxplot(data, x = "nativeStatusCode", y = "diam_t",
             color = "nativeStatusCode", palette = "jco")

f6=ggboxplot(data, x = "nativeStatusCode", y = "srl_t",
             color = "nativeStatusCode", palette = "jco")

f7=ggboxplot(data, x = "nativeStatusCode", y = "rmf_t",
             color = "nativeStatusCode", palette = "jco")

f8=ggboxplot(data, x = "nativeStatusCode", y = "rd_t",
             color = "nativeStatusCode", palette = "jco")

f9=ggboxplot(data, x = "nativeStatusCode", y = "rtd_t",
             color = "nativeStatusCode", palette = "jco")

ggarrange(f1,f2,f3,f4,f5,f6,f7,f8,f9,common.legend = TRUE, legend = "bottom")



########### Complete Cases Analyses (using original data prior to imputation) ############

data = read.csv("Traits comparisons_Data_completecases.csv",header=T,stringsAsFactors = T, row.names = 1)
data$nativeStatusCode <- factor(data$nativeStatusCode, levels = c("N","I"))

## Checking normality

hist(data$rootCcontent)
hist(data$rootNcontent)#
hist(data$rootCNratio)#
hist(data$mycorrhizaIntensity)#
hist(data$rootDiameter)#
hist(data$SRL)#
hist(data$RMF)
hist(data$rootingDepth)#
hist(data$rootTissueDensity)#

## Transformations

b <- boxcox(lm(data$rootNcontent ~ 1))
lambda <- b$x[which.max(b$y)]
data$rootN_t = (data$rootNcontent ^ lambda - 1) / lambda
hist(data$rootN_t)

b <- boxcox(lm(data$rootCNratio ~ 1))
lambda <- b$x[which.max(b$y)]
data$rootCN_t = (data$rootCNratio ^ lambda - 1) / lambda
hist(data$rootCN_t)

b <- boxcox(lm(data$mycorrhizaIntensity+1 ~ 1))
lambda <- b$x[which.max(b$y)]
data$rootMI_t = (data$mycorrhizaIntensity+1 ^ lambda - 1) / lambda
hist(data$rootMI_t)

b <- boxcox(lm(data$rootDiameter ~ 1))
lambda <- b$x[which.max(b$y)]
data$diam_t = (data$rootDiameter ^ lambda - 1) / lambda
hist(data$diam_t)

b <- boxcox(lm(data$SRL ~ 1))
lambda <- b$x[which.max(b$y)]
data$srl_t = (data$SRL ^ lambda - 1) / lambda
hist(data$srl_t)

b <- boxcox(lm(data$rootingDepth ~ 1))
lambda <- b$x[which.max(b$y)]
data$rd_t = (data$rootingDepth ^ lambda - 1) / lambda
hist(data$rd_t)

b <- boxcox(lm(data$rootTissueDensity ~ 1))
lambda <- b$x[which.max(b$y)]
data$rtd_t = (data$rootTissueDensity ^ lambda - 1) / lambda
hist(data$rtd_t)

# Root C

rootclm <- lm(rootCcontent~nativeStatusCode,data=data)
summary(rootclm)

rootcphylolm <- phylolm(rootCcontent~nativeStatusCode,data=data,phy=myTree)
summary(rootcphylolm)

plot(data$rootCcontent,rootcphylolm$residuals)

# root N

rootNlm <- lm(rootNcontent~nativeStatusCode,data=data)
summary(rootNlm)

rootNphylolm <- phylolm(rootN_t~nativeStatusCode,data=data,phy=myTree)
summary(rootNphylolm)

plot(data$rootN_t,rootNphylolm$residuals)#

# root CN*

rootCNlm <- lm(rootCNratio~nativeStatusCode,data=data)
summary(rootCNlm)

rootCNphylolm <- phylolm(rootCN_t~nativeStatusCode,data=data,phy=myTree)
summary(rootCNphylolm)

plot(data$rootCN_t,rootCNphylolm$residuals)#


# mycorrhizae colonization Intensity

MIlm <- lm(mycorrhizaIntensity ~ nativeStatusCode,data=data)
summary(MIlm)

rootMIphylolm <- phylolm(rootMI_t~nativeStatusCode,data=data,phy=myTree)
summary(rootMIphylolm)

plot(data$rootMI_t,rootMIphylolm$residuals)


# Root Diameter

diamlm <- lm(rootDiameter ~ nativeStatusCode,data=data)
summary(diamlm)

diamphylolm <- phylolm(diam_t~nativeStatusCode,data=data,phy=myTree)
summary(diamphylolm)

plot(data$diam_t,diamphylolm$residuals)

# SRL***

srllm <- lm(SRL ~ nativeStatusCode,data=data)
summary(srllm)

srlphylolm <- phylolm(srl_t~nativeStatusCode,data=data,phy=myTree)
summary(srlphylolm)

plot(data$srl_t,srlphylolm$residuals)


# RMF***

rmflm <- lm(RMF ~ nativeStatusCode,data=data)
summary(rmflm)

rmfphylolm <- phylolm(RMF~nativeStatusCode,data=data,phy=myTree)
summary(rmfphylolm)

plot(data$rmf_t,rmfphylolm$residuals)


# Rooting Depth

rdlm <- lm(rootingDepth~ nativeStatusCode,data=data)
summary(rdlm)

rdphylolm <- phylolm(rd_t~nativeStatusCode,data=data,phy=myTree)
summary(rdphylolm)

plot(data$rd_t,rdphylolm$residuals)


# Root Tissue Density

rtdlm <- lm(rootTissueDensity ~ nativeStatusCode,data=data)
summary(rtdlm)

rtdphylolm <- phylolm(rtd_t ~ nativeStatusCode,data=data,phy=myTree)
summary(rtdphylolm)

plot(data$rtd_t,rtdphylolm$residuals)


## Plot results

f1=ggboxplot(data, x = "nativeStatusCode", y = "rootCcontent",
             color = "nativeStatusCode", palette = "jco")

f2=ggboxplot(data, x = "nativeStatusCode", y = "rootNcontent",
             color = "nativeStatusCode", palette = "jco")

f3=ggboxplot(data, x = "nativeStatusCode", y = "rootCN_t",
             color = "nativeStatusCode", palette = "jco")

f4=ggboxplot(data, x = "nativeStatusCode", y = "rootMI_t",
             color = "nativeStatusCode", palette = "jco")

f5=ggboxplot(data, x = "nativeStatusCode", y = "diam_t",
             color = "nativeStatusCode", palette = "jco")

f6=ggboxplot(data, x = "nativeStatusCode", y = "srl_t",
             color = "nativeStatusCode", palette = "jco")

f7=ggboxplot(data, x = "nativeStatusCode", y = "rmf_t",
             color = "nativeStatusCode", palette = "jco")

f8=ggboxplot(data, x = "nativeStatusCode", y = "rd_t",
             color = "nativeStatusCode", palette = "jco")

f9=ggboxplot(data, x = "nativeStatusCode", y = "rtd_t",
             color = "nativeStatusCode", palette = "jco")

ggarrange(f1,f2,f3,f4,f5,f6,f7,f8,f9,common.legend = TRUE, legend = "bottom")


############ Trait correlation matrix #################

traits=na.omit(cbind(data[3],data[13:20]))
cor.traits=cor(traits)
corrplot(cor.traits)


