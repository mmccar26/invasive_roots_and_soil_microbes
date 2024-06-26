#===load libraries======
library(ggplot2)
library(tidyverse)
library(vegan)
library(eulerr)

#=====import data======
#relative pathname
plt <- file.path(".", "Data", "CWM_traits_data.csv")

##Load data
cwm<-read_csv(plt)%>%
  mutate(rootC_cwm_pc = rootC_cwm * pcentCover,
         rootN_cwm_pc = rootN_cwm * pcentCover,
         CNratio_cwm_pc = CNratio_cwm * pcentCover,
         MI_cwm_pc = MI_cwm * pcentCover,
         diameter_cwm_pc = diameter_cwm * pcentCover,
         SRL_cwm_pc = SRL_cwm * pcentCover,
         RMF_cwm_pc = RMF_cwm * pcentCover,
         RD_cwm_pc = RD_cwm * pcentCover,
         RTD_cwm_pc = RTD_cwm * pcentCover)%>%
 filter(pcentCover > 20)
 
#natives
natives<-
  cwm%>%
  filter(Status == "N")%>%
  group_by(plotID, Status)%>%
  summarize_at(vars(rootC_cwm_pc:RTD_cwm_pc), mean, na.rm = TRUE)%>%
  na.omit()

#invasives
invasives<-
  cwm%>%
  filter(Status == "I")%>%
  filter(pcentCover > 20)%>%
  group_by(plotID, Status)%>%
  summarize_at(vars(rootC_cwm_pc:RTD_cwm_pc), mean, na.rm = TRUE)%>%
  na.omit()

#=====import data======
#relative pathname
mic <- file.path(".", "Data", "neon_soil_all.csv")

microbes<-read_csv(mic)%>%
  group_by(plotID)%>%
  summarize_at(vars(c8To0Concentration:soilInWaterpH), mean, na.rm = TRUE)

##combine all data to match labels
mic_invas<-merge(microbes, invasives, by = "plotID")

mic_natives_invas<- merge(mic_invas, natives,  by='plotID')

##multivariate analysis
##Microbial community separation
spe<-
  mic_natives_invas%>%
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
            d15N:RTD_cwm_pc.y)
  )

#count NAs
colSums(is.na(spe))

##distance measure for microbes
spe.hel <- decostand(spe, method = "hellinger")

#plant environmental dataset
invas<-
  mic_natives_invas%>%
  select(rootC_cwm_pc.x:RTD_cwm_pc.x)%>%
  log()


mod1 <- rda(decostand(spe, "hel") ~., invas)
ordistep(mod1, direction = c("both"), Pout = 0.1, Pin = 0.1, step=100, perm.max=10000)
anova(mod1, by= "term", perm = 200) 
mod1$anova


#Native
nat<-
  mic_natives_invas%>%
  select(rootC_cwm_pc.y:RTD_cwm_pc.y)%>%
  log()

mod2 <- rda(decostand(spe, "hel") ~., nat)
ordistep(mod2, direction = c("both"), Pout = 0.1, Pin = 0.1, step=100, perm.max=10000)
anova(mod2, by= "term", perm = 200) 

#covariation
spe.part.all <- varpart(spe, nat, invas, transfo = 'hel')

summary(spe.part.all)
plot(spe.part.all, bg=2:4)

spe.part.all$part  # access results!

plot(spe.part.all,
     Xnames = c("Natives", "Invasives"), # name the partitions
     bg = c("#56B4E9", "firebrick"), alpha = 80, # colour the circles
     digits = 2, # only show 2 digits
     cex = 1.5)

aFrac <- rda(spe.hel, nat)
anova(aFrac)
## RsquareAdj gives the same result as component [a] of varpart

RsquareAdj(aFrac)
#significance testing
anova.cca(rda(spe.hel, nat))
anova.cca(rda(spe.hel, invas, nat))

##Variation paritiioning varplot
VennDiag <- euler(c("Natives" = 0.124, "Invasives" = 0.271, "Natives&Invasives" = 0.083))

plot(VennDiag, counts = TRUE, font=5, cex=1, alpha=0.5, labels = identical(legend, FALSE),
     fill=c("#56B4E9", "red1"))
