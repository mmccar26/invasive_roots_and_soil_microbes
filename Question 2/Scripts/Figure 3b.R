library(lme4)
library(lmerTest)
library(vegan)
library(dplyr)
library(dotwhisker)

########## Create dataset with PLFA columns and plot invasion status ###############

data=read.csv("neon_soil_subset.csv", header=T,stringsAsFactors = T)
NameList=read.csv("Included PLFA Variables.csv", header=T)

### Not invaded plots --> No Invasives & No Naturalized
data_notinvaded = data[is.na(data$I)&is.na(data$NI),]
data_notinvaded = cbind(data_notinvaded,plotStatus=rep("Not Invaded",nrow(data_notinvaded)))

### Moderately invaded plots --> Invasive species occupy less than 50% total cover 
data_modinvaded = data[data$I<50&!is.na(data$I),]
data_modinvaded = cbind(data_modinvaded,plotStatus=rep("Moderately Invaded",nrow(data_modinvaded)))

### Highly invaded plots --> Invasive species occupy 50% total cover or greater
data_highinvaded = data[data$I>=50&!is.na(data$I),]
data_highinvaded = cbind(data_highinvaded,plotStatus=rep("Highly Invaded",nrow(data_highinvaded)))

### Join rows for all invasion statuses in the same data frame
plotstatusdata = rbind(data_highinvaded,data_modinvaded,data_notinvaded)
plotstatusdata$plotStatus = as.factor(plotstatusdata$plotStatus)

### Add PLFA columns
data_pca <- cbind(plotstatusdata[1:6],plotstatusdata[8],log1p(plotstatusdata[,colnames(plotstatusdata) %in% NameList$PLFA.lipid]),plotstatusdata[76:81]) 

########### Multivariate homogeneity of variance with unique plots ########## 

### Removing data from older years for plots with multiple years of data ###
data_pca_sub <- data_pca %>% distinct(plotID, .keep_all = T)

### Bray-Curtis distances between samples
dis <- vegdist(data_pca_sub[,8:48],method="bray")

### Specify groups
groups <- data_pca_sub$plotStatus

### Perform betadisper
bdisper_plots <- betadisper(dis, groups, type = "median", bias.adjust = FALSE)

bdisper_plots

### Test significance of dispersion difference
anova(bdisper_plots)

### Perform Tukey test
tukeyresults=TukeyHSD(bdisper_plots)
tukeyresults

### Plot results

results_df <- cbind.data.frame(term=c("Moderately Invaded-Highly Invaded","Not Invaded-Highly Invaded","Not Invaded-Moderately Invaded"),
                               estimate=tukeyresults$group[1:3,1],conf.low=tukeyresults$group[1:3,2],conf.high=tukeyresults$group[1:3,3])

results_df %>% dotwhisker::dwplot(dot_args = list(color="black",size=2),
                                  whisker_args = list(size=1,color="black"))  + theme_bw() + ggtitle("(b)")+
  theme(plot.title = element_text(family = "Times",size = 12),axis.text.x = element_text(family = "Times",size = 12),axis.text.y = element_text(family = "Times", size = 11),legend.position="none",panel.grid.minor.x=element_blank(),panel.grid.major.x=element_blank()) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2)

### Perform PCA to obtain ordination

pca_plots <- prcomp(data_pca_sub[8:48],scale. = T,center = T)

s <- summary(pca_plots)

### Obtain scores for PCs 1,2,and 3 for 3D ordination

PC1 = pca_plots$x[,1]
PC2 = pca_plots$x[,2]
PC3 = pca_plots$x[,3]

write.csv(cbind(PC1,PC2,PC3,plotStatus=data_pca_sub$plotStatus),"dataForOrdination.csv")
