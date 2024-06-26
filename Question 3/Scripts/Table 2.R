library(lme4)
library(lmerTest)
library(vegan)
library(dplyr)
library(adespatial)
library(dotwhisker)


################### Data for longitudinal study ##############################
data2=read.csv("paired_consecutiveyears_traitsandsoilcomposition.csv",header = T)

################### Select relevant columns ##################################

data_multipaired <- cbind(data2[1:13],data2[81],data2[83:91],(data2[,colnames(data2) %in% NameList$PLFA.lipid])) 

#### Version of data without missing values

data_multipaired_nona=na.omit(data_multipaired)

#### Create response matrix with PLFA columns

dependent_vars <- matrix(as.numeric(unlist(data_multipaired_nona[,24:64])),nrow=nrow(data_multipaired_nona))

############## Running PERMANOVA via adonis2 ####################

set.seed(123)

### Run adonis2
manova_plots <- adonis2(dependent_vars ~ soilMoisture+CNratio_s+SRL_s+RD_s+RTD_s,by="margin",data = data_multipaired_nona,na.action=na.omit,method="euclidean")
manova_plots

### Run variance partitioning
var <- varpart(dependent_vars, ~soilMoisture,~CNratio_s,~SRL_s,~RTD_s,data=data_multipaired_nona)
var
plot(var)