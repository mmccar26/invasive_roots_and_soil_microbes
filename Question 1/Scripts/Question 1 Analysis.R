
library(phylolm)
library(ggpubr)
library(adephylo)
library(ape)
library(phylobase)
library(lme4)
library(lmerTest)
library(corrplot)
library(MASS)

################ Uploading imputed traits dataset ##################

data = read.csv("Traits comparisons_Data.12.20.23.csv",header=T,row.names = 1)
data$nativeStatusCode <- factor(data$nativeStatusCode, levels = c("N","I"))

######################### Data preparation #########################

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


############### Performing phylogenetic regressions ################

### Upload phylogenetic tree

myTree <- read.tree("C:\\Users\\gnm\\Dropbox\\Roots&food-webs\\Root traits data\\Data imputation\\RPhylopars Inputs\\output_tree_based_on_plant_megatree.tre")

### Run phylogenetic comparisons

# Root C

rootcphylolm <- phylolm(rootCcontent~nativeStatusCode,data=data,phy=myTree)
summary(rootcphylolm)

plot(data$rootCcontent,rootcphylolm$residuals)

# root N

rootNphylolm <- phylolm(rootN_t~nativeStatusCode,data=data,phy=myTree)
summary(rootNphylolm)

plot(data$rootNcontent,rootNphylolm$residuals)

# root CN

rootCNphylolm <- phylolm(rootCN_t~nativeStatusCode,data=data,phy=myTree)
summary(rootCNphylolm)

plot(data$rootCN_t,rootCNphylolm$residuals)

# mycorrhizae colonization Intensity

rootMIphylolm <- phylolm(mycorrhizaIntensity~nativeStatusCode,data=data,phy=myTree)
summary(rootMIphylolm)

plot(data$rootMI_t,rootMIphylolm$residuals)

# Root Diameter

diamphylolm <- phylolm(diam_t~nativeStatusCode,data=data,phy=myTree)
summary(diamphylolm)

plot(data$diam_t,diamphylolm$residuals)

# SRL

srlphylolm <- phylolm(srl_t~nativeStatusCode,data=data,phy=myTree)
summary(srlphylolm)

plot(data$srl_t,srlphylolm$residuals)

# RMF

rmfphylolm <- phylolm(rmf_t~nativeStatusCode,data=data,phy=myTree)
summary(rmfphylolm)

plot(data$rmf_t,rmfphylolm$residuals)


# Rooting Depth

rdphylolm <- phylolm(rd_t~nativeStatusCode,data=data,phy=myTree)
summary(rdphylolm)

plot(data$rd_t,rdphylolm$residuals)

# Root Tissue Density

rtdphylolm <- phylolm(rtd_t ~ nativeStatusCode,data=data,phy=myTree)
summary(rtdphylolm)

plot(data$rtd_t,rtdphylolm$residuals)


## Plot results

f1=ggboxplot(data, x = "nativeStatusCode", y = "rootCcontent",
             color = "gray", palette = c("#56B4E9","firebrick"))+
  labs(x = NULL, y=NULL, subtitle = "Root C concentration\n(mg/g) n=783",colour = "Origin")+
  theme_bw()+
  theme(plot.subtitle = element_text(size = 10.5),axis.text.x = element_text(size = 10.5),
        axis.text.y = element_text(size = 10.5),legend.position="none")

f2=ggboxplot(data, x = "nativeStatusCode", y = "rootNcontent",
             color = "gray", palette = c("#56B4E9","firebrick"))+
  labs(x = NULL,y=NULL, subtitle = "Root N concentration\n(mg/g) n=785",colour = "Origin")+
  theme_bw()+
  theme(plot.subtitle = element_text(size = 10.5),axis.text.x = element_text(size = 10.5),
        axis.text.y = element_text(size = 10.5),legend.position="none")

f3=ggboxplot(data, x = "nativeStatusCode", y = "rootCNratio",
             color = "gray", palette = c("#56B4E9","firebrick"))+
  labs(x = NULL, y=NULL, subtitle = "C:N ratio\n(mg/mg) n=769",colour = "Origin")+
  theme_bw()+
  theme(plot.subtitle = element_text(size = 10.5),axis.text.x = element_text(size = 10.5),
        axis.text.y = element_text(size = 10.5),legend.position="none")

f4=ggboxplot(data, x = "nativeStatusCode", y = "mycorrhizaIntensity",
             color = "gray", palette = c("#56B4E9","firebrick"))+
  labs(x = NULL, y=NULL, subtitle = "Mycorrhizal colonization\n(%) n=787",colour = "Origin")+
  theme_bw()+
  theme(plot.subtitle = element_text(size = 10.5),axis.text.x = element_text(size = 10.5),
        axis.text.y = element_text(size = 10.5),legend.position="none")

f5=ggboxplot(data, x = "nativeStatusCode", y = "rootDiameter",
             color = "gray", palette = c("#56B4E9","firebrick"))+
  labs(x = NULL, y=NULL, subtitle = "Root diameter\n(mm) n=787",colour = "Origin")+
  theme_bw()+
  theme(plot.subtitle = element_text(size = 10.5),axis.text.x = element_text(size = 10.5),
        axis.text.y = element_text(size = 10.5),legend.position="none")

f6=ggboxplot(data, x = "nativeStatusCode", y = "SRL",
             color = "nativeStatusCode", palette = c("#56B4E9","firebrick"))+ 
  labs(x = NULL, y=NULL, subtitle = "Specific root length\n(m/g) n=786",colour = "Origin")+
  annotate("text", y = max(data$SRL,na.rm=T)*.95, x = 1.5, 
           label = "***", col="red")+
  theme_bw()+
  theme(plot.subtitle = element_text(size = 10.5),axis.text.x = element_text(size = 10.5),
        axis.text.y = element_text(size = 10.5),legend.position="none")

f7=ggboxplot(data, x = "nativeStatusCode", y = "RMF",
             color = "gray", palette = c("#56B4E9","firebrick"))+ 
  labs(x = NULL, y=NULL, subtitle = "Root mass fraction\n(g/g) n=772",colour = "Origin")+
  theme_bw()+
  theme(plot.subtitle = element_text(size = 10.5),axis.text.x = element_text(size = 10.5),
        axis.text.y = element_text(size = 10.5),legend.position="none")

f8=ggboxplot(data, x = "nativeStatusCode", y = "rootingDepth",
             color = "gray", palette = c("#56B4E9","firebrick"))+ 
  labs(x = NULL, y=NULL, subtitle = "Rooting depth\n(cm) n=783",colour = "Origin")+
  theme_bw()+
  theme(plot.subtitle = element_text(size = 10.5),axis.text.x = element_text(size = 10.5),
        axis.text.y = element_text(size = 10.5),legend.position="none")

f9=ggboxplot(data, x = "nativeStatusCode", y = "rootTissueDensity",
             color = "nativeStatusCode", palette = c("#56B4E9","firebrick"))+ 
  labs(x = NULL, y=NULL, subtitle = "Root tissue density\n(g/cubic cm) n=776",colour = "Origin")+ 
  annotate("text", y = max(data$rootTissueDensity,na.rm=T)*.95, x = 1.5, 
           label = "*",color="red")+
  theme_bw()+
  theme(plot.subtitle = element_text(size = 10.5),axis.text.x = element_text(size = 10.5),
        axis.text.y = element_text(size = 10.5),legend.position="none")


jpeg("C:\\Users\\gnm\\Dropbox\\Roots&food-webs\\Manuscript\\NatureEE Submission\\Q1Fig.jpeg", width = 6, height = 6, units = 'in', res = 800)

ggarrange(f6,f9,f1,f2,f3,f4,f5,f7,f8,common.legend = TRUE, legend = "none")

dev.off()