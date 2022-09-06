###Nunez-Mir and McCary, invasive root traits analysis
###Question 1 Analysis
##9 Sept 2022

######## Fixed-effects and Mixed-effects Modeling#################

library(r2glmm)
library(ggplot2)
library(sjPlot)
library(lme4)
library(lmerTest)
library(MASS)
library(ggpmisc)
library(effects)
library(ggpubr)

data3=read.csv("Data\\question1.csv", header = TRUE)

#### BoxCox Transformations

b <- boxcox(lm(data3$AvgOfphH2o ~ 1))
lambda <- b$x[which.max(b$y)]
soilph_t <- (data3$AvgOfphH2o ^ lambda - 1) / lambda

b <- boxcox(lm(data3$AvgOfnitrogenTot ~ 1))
lambda <- b$x[which.max(b$y)]
soilntot_t <- (data3$AvgOfnitrogenTot ^ lambda - 1) / lambda

b <- boxcox(lm(data3$AvgOfctonRatio ~ 1))
lambda <- b$x[which.max(b$y)]
soilcton_t <- (data3$AvgOfctonRatio ^ lambda - 1) / lambda

b <- boxcox(lm((data3$AvgOfpOxalate+1) ~ 1))
lambda <- b$x[which.max(b$y)]
soilphos_t <- (data3$AvgOfpOxalate ^ lambda - 1) / lambda

b <- boxcox(lm((data3$Allbacteria) ~ 1))
lambda <- b$x[which.max(b$y)]
bac_t <- (data3$Allbacteria ^ lambda - 1) / lambda

b <- boxcox(lm((data3$Fungi) ~ 1))
lambda <- b$x[which.max(b$y)]
fungi_t <- (data3$Fungi ^ lambda - 1) / lambda

b <- boxcox(lm((data3$Ratio) ~ 1))
lambda <- b$x[which.max(b$y)]
ratio_t <- (data3$Ratio ^ lambda - 1) / lambda

b <- boxcox(lm((data3$InvCover) ~ 1))
lambda <- b$x[which.max(b$y)]
invcover_t <- as.vector(scale((data3$InvCover ^ lambda - 1) / lambda))

b <- boxcox(lm((data3$InvFraction) ~ 1))
lambda <- b$x[which.max(b$y)]
invfrac_t <- as.vector(scale((data3$InvFraction ^ lambda - 1) / lambda))

covariates_t=scale(cbind(soilph_t,soilcton_t,soilntot_t,soilphos_t))

d=as.data.frame(cbind(as.numeric(invfrac_t),as.numeric(invcover_t),as.numeric(ratio_t), as.numeric(fungi_t), 
                      as.numeric(bac_t),Site=data3$siteID,Ecosystem=data3$Ecosystem))

###### Regressions

### Fungi

#Testing simple and mixed-effects models with and w/out soil property covariates
#--Simple
lm_fungi=lm(fungi_t~invfrac_t+soilph_t+soilcton_t+soilntot_t+soilphos_t)
summary(lm_fungi)
effx1 <- effect("invfrac_t", lm_fungi, partial.residuals=T)
plot(effx1, smooth.residuals=F)

lm_fungi2=lm(fungi_t~invfrac_t)
summary(lm_fungi2)
plot(invfrac_t,fungi_t)
abline(lm_fungi2)

AIC(lm_fungi,lm_fungi2)

#--Mixed effects
lmr_fungi=lmer(fungi_t~invfrac_t+covariates_t+(1|Site),data=d)
summary(lmr_fungi)

lmr_fungi2=lmer(fungi_t~invfrac_t + (1|Site),data=d)
sm=summary(lmr_fungi2)


AIC(lm_fungi,lm_fungi2,lmr_fungi,lmr_fungi2)

#--Plot best models
shapes <- c(17,8,15,18,19)
names(shapes) <- c("Agriculture", "Evergreen", "Grassland", "Mixed forest", "Wetland")
colors <- c("orange","#26453e","#39FF14", "#FF00FF","blue")

# Simple plot
plotf <- ggplot(d,aes(invfrac_t, fungi_t)) + 
            geom_point(aes(col=Ecosystem, shape=Ecosystem,size=2)) + 
            geom_smooth(method=lm, level=0.95, size=1.5, col="black",lty="dotted")+
            xlab("Invasive fraction")+
            ylab("Fungal biomass")+
            guides(size = "none")+
            scale_shape_manual(values= shapes)+
            scale_color_manual(values= colors)+
            guides(colour = guide_legend(override.aes = list(size=4)))+
            ggtitle("(a) Fungi") +
            stat_poly_eq(label.x = 0.85,
                         label.y = 0.95,aes(group=NA, label = paste(..rr.label..)))+
            stat_fit_glance(method = "lm", 
                  method.args = list(formula = fungi_t~invfrac_t),
                  label.x = 0.85,
                  label.y = 0.9,
                  aes(group=NA,label = paste("italic(P)*\"-value = \"*", 
                                             signif(after_stat(p.value), digits = 1),
                                             sep = "")),parse = TRUE) +
            theme_bw()+ theme( panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                    axis.text=element_text(size=10),
                    axis.title=element_text(size=11), legend.title = element_blank(),legend.text = element_text(size=11))
plotf
# Mixed-effects plot
d$fit <- predict(lmr_fungi2) 
plotm1 <- ggplot(d,aes(invfrac_t, fungi_t, group=Site, col=Site, shape=Ecosystem)) + 
  geom_line(aes(y=fit), size=1.5) +
  geom_point(aes(size=2)) + 
  scale_shape_manual(values= shapes)+
  geom_abline(intercept=sm$coefficients[1],slope=sm$coefficients[2],size=2,lty="dotted") +
  geom_text(label=paste0("Estimate = 0.07"),x=1,y=1.7,col="black")+
  geom_text(label=paste0("italic('P')~'-value = 0.6'"),parse=T,x=1,y=1.4,col="black")+
  xlab("Invasive fraction")+
  ylab("Fungal biomass")+
  ggtitle("(a) Fungi") +
  guides(size = "none")+
  theme_bw() + theme( panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                        axis.text=element_text(size=10),
                        axis.title=element_text(size=11), legend.title = element_blank(),legend.text = element_text(size=11))



### Bacteria

#Testing simple and mixed-effects models with and w/out soil property covariates
#--Simple
lm_bac=lm(bac_t~invfrac_t+covariates_t)
summary(lm_bac)

lm_bac2=lm(bac_t~invfrac_t)
summary(lm_bac2)
plot(invfrac_t,bac_t)
abline(lm_bac2)

AIC(lm_bac,lm_bac2)

#--Mixed effects
lmr_bac=lmer(bac_t~invfrac_t+covariates_t+(1|Site),data=d)
summary(lmr_bac)

lmr_bac2=lmer(bac_t~invfrac_t + (1|Site),data=d)
sm=summary(lmr_bac2)
sm

AIC(lm_bac,lm_bac2,lmr_bac,lmr_bac2)

#--Plot best models
shapes <- c(17,8,15,18,19)
names(shapes) <- c("Agriculture", "Evergreen", "Grassland", "Mixed forest", "Wetland")

# Simple plot
plotb <- ggplot(d,aes(invfrac_t, bac_t)) + 
  geom_point(aes(col=Ecosystem, shape=Ecosystem,size=2)) + 
  geom_smooth(method=lm, level=0.95, size=1.5, col="black",lty="dotted")+
  xlab("Invasive fraction")+
  ylab("Bacterial biomass")+
  guides(size = F)+
  scale_shape_manual(values= shapes)+
  scale_color_manual(values= colors)+
  guides(colour = guide_legend(override.aes = list(size=4)))+
  ggtitle("(b) Bacteria") +
  stat_poly_eq(label.x = 0.85,
               label.y = 0.15,aes(group=NA, label = paste(..rr.label..)))+
  stat_fit_glance(method = "lm", 
                  method.args = list(formula = bac_t~invfrac_t),
                  label.x = 0.85,
                  label.y = 0.1,
                  aes(group=NA,label = paste("italic(P)*\"-value = \"*", 
                                             signif(after_stat(p.value), digits = 1),
                                             sep = "")),parse = TRUE) +
  theme_bw()+ theme( panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=10),
                     axis.title=element_text(size=11), legend.title = element_blank(),legend.text = element_text(size=11))
plotb
# Mixed-effects plot
d$fit <- predict(lmr_bac2) 
plotm2 <- ggplot(d,aes(invfrac_t, bac_t, group=Site, col=Site, shape=Ecosystem)) + 
  geom_line(aes(y=fit), size=1.5) +
  geom_point(aes(size=2)) + 
  scale_shape_manual(values= shapes)+
  geom_abline(intercept=sm$coefficients[1],slope=sm$coefficients[2],size=2,lty="dotted") +
  geom_text(label=paste0("Estimate = 0.08"),x=1,y=2.1,col="black")+
  geom_text(label=paste0("italic('P')~'-value = 0.4'"),parse=T,x=1,y=1.95,col="black")+
  xlab("Invasive fraction")+
  ylab("Bacterial biomass")+
  ggtitle("(b) Bacteria") +
  guides(size = "none")+
  theme_bw()+ theme( panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=10),
                     axis.title=element_text(size=11), legend.title = element_blank(),legend.text = element_text(size=11))




### Fungi:Bacteria Ratio

#Testing simple and mixed-effects models with and w/out soil property covariates
#--Simple
lm_ratio=lm(ratio_t~invfrac_t+covariates_t)
summary(lm_ratio)

lm_ratio2=lm(ratio_t~invfrac_t)
summary(lm_ratio2)
plot(invfrac_t,ratio_t)
abline(lm_ratio2)

AIC(lm_ratio,lm_ratio2)

#--Mixed effects
lmr_ratio=lmer(ratio_t~invfrac_t+covariates_t+(1|Site),data=d)
summary(lmr_ratio)

lmr_ratio2=lmer(ratio_t~invfrac_t + (invfrac_t|Site),data=d)
sm=summary(lmr_ratio2)
sm

AIC(lm_ratio,lm_ratio2,lmr_ratio,lmr_ratio2)

#--Plot best models
shapes <- c(17,8,15,18,19)
names(shapes) <- c("Agriculture", "Evergreen", "Grassland", "Mixed forest", "Wetland")

# Simple plot
plotf2b <- ggplot(d,aes(invfrac_t, ratio_t)) + 
  geom_point(aes(col=Ecosystem, shape=Ecosystem,size=2)) + 
  geom_smooth(method=lm, level=0.95, size=1.5, col="black")+
  xlab("Invasive fraction")+
  ylab("Fungi:bacteria ratio")+
  guides(size = "none")+
  scale_shape_manual(values= shapes)+
  scale_color_manual(values= colors)+
  ggtitle("(c) Fungi:bacteria") +
  guides(colour = guide_legend(override.aes = list(size=4)))+
  stat_poly_eq(label.x = 0.85,
               label.y = 0.95,aes(group=NA, label = paste(..rr.label..)))+
  stat_fit_glance(method = "lm", 
                  method.args = list(formula = ratio_t~invfrac_t),
                  label.x = 0.85,
                  label.y = 0.9,
                  aes(group=NA,label = paste("italic(P)*\"-value = \"*", 
                                             signif(after_stat(p.value), digits = 2),
                                             sep = "")),parse = TRUE) +
  theme_bw()+ theme( panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=10),
                     axis.title=element_text(size=11), legend.title = element_blank(),legend.text = element_text(size=11))
plotf2b

# Mixed-effects plot
d$fit <- predict(lmr_ratio2) 
plotm3 <- ggplot(d,aes(invfrac_t, ratio_t, group=Site, col=Site, shape=Ecosystem)) + 
  geom_line(aes(y=fit), size=1.5) +
  geom_point(aes(size=2)) + 
  scale_shape_manual(values= shapes)+
  geom_abline(intercept=sm$coefficients[1],slope=sm$coefficients[2],size=2,lty="dotted") +
  geom_text(label=paste0("Estimate = -0.19"),x=1,y=0.1,col="black")+
  geom_text(label=paste0("italic('P')~'-value = 0.3'"),parse=T,x=1,y=-0.06,col="black")+
  xlab("Invasive fraction")+
  ylab("Fungi:bacteria ratio")+
  ggtitle("(c) Fungi:bacteria") +
  guides(size = "none")+
  theme_bw()+ theme( panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=10),
                     axis.title=element_text(size=11), legend.title = element_blank(),legend.text = element_text(size=11))



######### Exporting plots for pub
jpeg("C:\\Users\\gnune\\Dropbox\\Roots&food-webs\\Manuscript\\Figures\\Q1Regressions.jpeg", width = 14, height = 5, units = 'in', res = 800)

fig = ggarrange(plotf,plotb,plotf2b,common.legend=T,legend="bottom",ncol=3,nrow=1)
fig

dev.off()


jpeg("C:\\Users\\gnune\\Dropbox\\Roots&food-webs\\Manuscript\\Figures\\Q1Regressions_appendix.jpeg", width = 14, height = 5, units = 'in', res = 800)

fig = ggarrange(plotm1,plotm2,plotm3,common.legend=T,ncol=3,nrow=1)
fig

dev.off()



