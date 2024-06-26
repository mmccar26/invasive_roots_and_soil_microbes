##Main results R code
##Question 3

#===load libraries======
library(ggplot2)
library(tidyverse)
library(vegan)
library(ggord)
library(BiodiversityR)
library(ggrepel)

#=====import data======
#relative pathname
plt <- file.path(".", "Data", "CWM_traits_data.csv")

##Load data
cwm<-read_csv(plt)%>%
  mutate(Root_C = rootC_cwm * pcentCover,
         Root_N = rootN_cwm * pcentCover,
         CN_ratio = CNratio_cwm * pcentCover,
         Myco_I = MI_cwm * pcentCover,
         Root_diameter = diameter_cwm * pcentCover,
         SRL = SRL_cwm * pcentCover,
         RMF = RMF_cwm * pcentCover,
         RD = RD_cwm * pcentCover,
         RTD = RTD_cwm * pcentCover)%>%
  filter(pcentCover > 20)

#invasives
invasives<-
  cwm%>%
  filter(Status == "I")%>%
  group_by(plotID, Status)%>%
  summarize_at(vars(Root_C:RTD), mean, na.rm = TRUE)%>%
  na.omit()

#=====import data======
#relative pathname
mic <- file.path(".", "Data", "neon_soil_all.csv")

microbes<-read_csv(mic)%>%
  group_by(plotID)%>%
  summarize_at(vars(c8To0Concentration:soilInWaterpH), mean, na.rm = TRUE)

##combine all data to match labels
mic_invas<- merge(microbes, invasives,  by='plotID')

##multivariate analysis
##Microbial community separation
spe<-
  mic_invas%>%
  select(-c("plotID",
            "c10To0Concentration",
            "c11To0Concentration", 
            "c8To0Concentration",
            "lipid2OH12To0Concentration",
            "trans18To1n9Concentration",
            "lipid10Methyl18To1Concentration",
            "lipid3OH12To0Concentration",
            "c16To1n7Concentration",
            "lipid10Methyl17To1Concentration",
            "trans18To2n912Concentration",
            "c18To1n13Concentration",
            "c18To3n3Concentration",
            "c19To0Concentration",
            "lipid2OH14To0Concentration",
            "c19To1Cis10Concentration",
            "c20To3n3Concentration",
            "lipid3OH14To0Concentration",
            "c17To1n7Concentration",
            "c18To0Concentration",
            "c21To0Concentration",
            "totalLipidConcentration",
            d15N:RTD)
  )

#count NAs
colSums(is.na(spe))

##distance measure for microbes
spe.hel <- decostand(spe, method = "hellinger")

#plant environmental dataset
invas<-
  mic_invas%>%
  select(Root_C:RTD)%>%
  log()

mod1 <- rda(decostand(spe, "hel") ~., invas)
mod1 <- rda(decostand(spe, "hel") ~.-Root_C, invas)
ordistep(mod1, direction = c("back"), Pout = 0.1, Pin = 0.1, step=100, perm.max=10000)
anova(mod1, by= "term", perm = 200) 
RsquareAdj(mod1)

vif.cca(mod1)

##RDA plot
##set the ggplot2 theme
BioR.theme <- theme(
  panel.background = element_blank(),
  panel.border = element_rect(fill = NA, colour = "black", size = 2),
  panel.grid = element_blank(),
  axis.line = element_line("gray25"),
  text = element_text(size = 12),
  axis.text = element_text(size = 10, colour = "gray25"),
  axis.title = element_text(size = 14, colour = "gray25"),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 14),
  legend.key = element_blank())

##plot
#Preset conditions for plot
plot1 <- ordiplot(mod1, choices=c(1,2))##set axes 

sites.long2 <- sites.long(plot1, env.data=mic_invas)

axis.long2 <- axis.long(mod1, choices=c(1, 2))

spec.envfit <- envfit(plot1, env=invas)

spec.data.envfit <- data.frame(r=spec.envfit$vectors$r, p=spec.envfit$vectors$pvals)

#set vectors
spp.scrs <- as.data.frame(scores(spec.envfit, display = "vectors"))
spp.scrs <- cbind(spp.scrs, Environ = rownames(spp.scrs))%>%
  slice(6, 7, 8, 9)

#delete outlier
sites.long2<-sites.long2[-41,]

##ggplot2 ordination
plotgg2 <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab(axis.long2[1, "label"]) +
  ylab(axis.long2[2, "label"]) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL), limits = c(-.8, .8)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL), limits = c(-0.8, .8)) +    
  geom_point(data=sites.long2, 
             aes(x=axis1, y=axis2), 
             size=5, shape = 16, color = "red") +
  geom_segment(data=spp.scrs,
               aes(x = 0, xend = RDA1, y = 0, yend = RDA2),
               arrow = arrow(length = unit(0.25, "cm")), color = "black", size = 1) +
  geom_text_repel(data=spp.scrs,  aes(x = RDA1, y = RDA2, label = Environ),
                  size = 6)+
  BioR.theme +
  coord_fixed(ratio=1)



#####################################
##natives
#=====import data======
#relative pathname
plt <- file.path(".", "Data", "CWM_traits_data.csv")

##Load data
cwm<-read_csv(plt)%>%
  mutate(Root_C = rootC_cwm * pcentCover,
         Root_N = rootN_cwm * pcentCover,
         CN_ratio = CNratio_cwm * pcentCover,
         Myco_I = MI_cwm * pcentCover,
         Root_diameter = diameter_cwm * pcentCover,
         SRL = SRL_cwm * pcentCover,
         RMF = RMF_cwm * pcentCover,
         RD = RD_cwm * pcentCover,
         RTD = RTD_cwm * pcentCover)%>%
  filter(pcentCover > 20)

#native
natives<-
  cwm%>%
  filter(Status == "N")%>%
  sample_n(45)%>%
  group_by(plotID, Status)%>%
  summarize_at(vars(Root_C:RTD), mean, na.rm = TRUE)%>%
  na.omit()

#copy data
write.csv(natives, "C:/Users/mattm/Dropbox/Matt's folder/Rice University/Research/Invasive root traits/Analysis/Invasive_root_traits/Data/natives.csv", row.names=F)

#import derived data
nat.p <- file.path(".", "Data", "natives.csv")
##Load data
natives<-read_csv(nat.p)

#=====import data======
#relative pathname
mic <- file.path(".", "Data", "neon_soil_all.csv")

microbes<-read_csv(mic)%>%
  group_by(plotID)%>%
  summarize_at(vars(c8To0Concentration:soilInWaterpH), mean, na.rm = TRUE)

##combine all data to match labels
mic_natives<- merge(microbes, natives,  by='plotID')

##multivariate analysis
##Microbial community separation
spe<-
  mic_natives%>%
  select(-c("plotID",
            "c10To0Concentration",
            "c11To0Concentration", 
            "c8To0Concentration",
            "lipid2OH12To0Concentration",
            "trans18To1n9Concentration",
            "lipid10Methyl18To1Concentration",
            "lipid3OH12To0Concentration",
            "c16To1n7Concentration",
            "lipid10Methyl17To1Concentration",
            "trans18To2n912Concentration",
            "c18To1n13Concentration",
            "c18To3n3Concentration",
            "c19To0Concentration",
            "lipid2OH14To0Concentration",
            "c19To1Cis10Concentration",
            "c20To3n3Concentration",
            "lipid3OH14To0Concentration",
            "c17To1n7Concentration",
            "c18To0Concentration",
            "c21To0Concentration",
            "totalLipidConcentration",
            d15N:RTD)
  )

#count NAs
colSums(is.na(spe))

##distance measure for microbes
spe.hel <- decostand(spe, method = "hellinger")

#plant environmental dataset
nat<-
  mic_natives%>%
  select(Root_C:RTD)%>%
  log()

mod2 <- rda(decostand(spe, "hel") ~., nat)
mod2 <- rda(decostand(spe, "hel") ~. -Root_C, nat)
ordistep(mod2, direction = c("back"), Pout = 0.1, Pin = 0.1, step=100, perm.max=10000)
anova(mod2, by= "term", perm = 200) 
RsquareAdj(mod2)

vif.cca(mod2)

##RDA plot
##set the ggplot2 theme
BioR.theme <- theme(
  panel.background = element_blank(),
  panel.border = element_rect(fill = NA, colour = "black", size = 2),
  panel.grid = element_blank(),
  axis.line = element_line("gray25"),
  text = element_text(size = 12),
  axis.text = element_text(size = 10, colour = "gray25"),
  axis.title = element_text(size = 14, colour = "gray25"),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 14),
  legend.key = element_blank())

nat1<-
  nat%>%
  select(CN_ratio, Myco_I)

##plot
#Preset conditions for plot
plot1 <- ordiplot(mod2, choices=c(1,2))##set axes 

sites.long2 <- sites.long(plot1, env.data=nat1)

axis.long2 <- axis.long(plot1, choices=c(1, 2))

spec.envfit <- envfit(plot1, env=nat1)

spec.data.envfit <- data.frame(r=spec.envfit$vectors$r, p=spec.envfit$vectors$pvals)

#set vectors, 
spp.scrs <- as.data.frame(scores(spec.envfit, display = "vectors"))
spp.scrs <- cbind(spp.scrs, Environ = rownames(spp.scrs))%>%
 slice(2, 3, 4, 7, 9)

##ggplot2 ordination
native_plot <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab(axis.long2[1, "label"]) +
  ylab(axis.long2[2, "label"]) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL), limits = c(-.8, .8)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL), limits = c(-.8, .8)) +    
  geom_point(data=sites.long2, 
             aes(x=axis1, y=axis2), 
             size=5, shape = 17, color = "#56B4E9") +
  geom_segment(data=spp.scrs,
               aes(x = 0, xend = RDA1, y = 0, yend = RDA2),
               arrow = arrow(length = unit(0.25, "cm")), color = "black", size = 1) +
  geom_text_repel(data=spp.scrs,  aes(x = RDA1, y = RDA2, label = Environ),
                  size = 6)+
  BioR.theme +
  coord_fixed(ratio=1)


