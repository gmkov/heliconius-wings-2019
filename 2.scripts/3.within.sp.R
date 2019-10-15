##### Altitude and life-history shape the evolution of Heliconius wings ####
#### Gabriela Montejo-Kovacevich 2019 ####
####### Heliconius wing shape/size PGLS altitude models ######
#### packages #####
rm(list=ls())
dev.off()
library(dplyr)
library(data.table)
library(ape)
library(nlme)
library(rptR)
library(rr2)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(lattice)
library(lme4)
library(lmSupport)
library(car)
library(relaimpo)
library(rgl)
library(forcats)
library(cowplot)
library(ggExtra)
library(AICcmodavg)
library(caper)
library(MuMIn)
library(broom)
library(broom.mixed)
library(pwr)



#### data  #####
setwd("../heliconius-wings-2019/")
master <- read.csv("1.data/master.analyses.csv")
mean.df <- read.csv("1.data/mean.df.analyses.csv")
nat.hist <- read.csv("1.data/species.nat.hist.csv")
pruned.tree <- read.tree("1.data/pruned.tree.phylo")

#### data handling, functions  #####
master$sp.short <- substr(master$species,0,3)
master.sex <- subset(master, sex!="")
rownames(mean.df) <- mean.df$genus.sp
era.mel <- subset(master, (species=="erato"|species=="melpomene")&wing.pattern!="hybrid"&wing.pattern!="hybrid")
era.mel.sex <- subset(era.mel, (species=="erato"|species=="melpomene")&sex!="")
era <- subset(era.mel, species=="erato")
mel <- subset(era.mel, species=="melpomene")
capitalize <- function(string) {
  substr(string, 1, 1) <- toupper(substr(string, 1, 1))
  string
}

############ 1. Plots era/mel #############
# Fig. 5 raw data shape vs. alt 
options(digits = 3, scipen = -2)
shape.raw.alt.mel.era <- ggplot(era.mel.sex, aes(x=altitude, y=aspect.ratio,colour=species, shape=sex, linetype=sex )) +
  geom_point(size=1, alpha=.3)  +
  #geom_text(label=mel.master$cam.id)+
  geom_smooth(method='gam',formula =  y ~ poly(x, 1), aes(fill=species))+
  scale_color_manual(name="Species", values = c("melpomene"="#E69F00","erato"="#56B4E9"), labels=c(expression(italic("H. erato")), expression(italic("H. melpomene"))))+
  scale_fill_manual(name="Species", values = c("melpomene"="#E69F00","erato"="#56B4E9"), labels=c(expression(italic("H. erato")), expression(italic("H. melpomene"))))+
  #scale_shape_manual(name="Species", values = c("melpomene"=16,"erato"=17), labels=c(expression(italic("H. erato")), expression(italic("H. melpomene")))) +
  xlab("Altitude (m)")+ 
  scale_linetype_manual(name="Sex", 
                          breaks=c( "female","male"), values=c("dashed","solid"),
                          labels = c("female","male"))+
  ylab("Wing Shape (aspect ratio)") + 
  theme_classic()+
  stat_cor(show.legend = FALSE, label.y=c(1.82,1.85,1.82,1.85), 
           aes(label =paste(..r.label.., cut(..p.., 
                                              breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                              labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), 
                                              sep = "~")))+
  facet_wrap(~species, labeller =labeller(species=c("erato"=c("H. erato"),"melpomene"="H. melpomene")))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),
        axis.line=element_blank(),
        panel.background = element_blank(), legend.text = element_text(size=12), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(), axis.text=element_text(size=12),
        strip.text = element_text(size=12, face="bold.italic"),plot.margin=grid::unit(c(3,10,3,3), "mm"),
        panel.grid.minor = element_blank(), axis.title=element_text(size=14,face="bold"),
        legend.text.align = 0, legend.position = "bottom"); shape.raw.alt.mel.era 

# Fig. S10 - raw data size vs alt
size.raw.alt.mel.era <- ggplot(era.mel.sex, aes(x=altitude, y=area.mm2,colour=species )) +
  geom_point(size=.5, alpha=.4)  +
  #geom_text(label=mel.master$cam.id)+
  #facet_wrap(~sp.subsp)+
  geom_smooth(method='gam',formula =  y ~ poly(x, 1), aes(fill=species))+
  scale_color_manual(name="Species", values = c("melpomene"="#E69F00","erato"="#56B4E9"), labels=c(expression(italic("H. erato")), expression(italic("H. melpomene"))))+
  scale_fill_manual(name="Species", values = c("melpomene"="#E69F00","erato"="#56B4E9"), labels=c(expression(italic("H. erato")), expression(italic("H. melpomene"))))+
  #scale_shape_manual(name="Species", values = c("melpomene"=16,"erato"=17), labels=c(expression(italic("H. erato")), expression(italic("H. melpomene")))) +
  xlab("Altitude (m)")+ 
  ylab(bquote(bold('Wing area ('*mm^2*')'))) + 
  theme_classic()+ stat_cor(show.legend = FALSE)+
  facet_wrap(~sex, labeller =labeller(sex=capitalize))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),
        panel.background = element_blank(), legend.text = element_text(size=12), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(), axis.text=element_text(size=12),
        strip.text = element_text(size=14, face="bold"),plot.margin=grid::unit(c(3,10,3,3), "mm"),
        panel.grid.minor = element_blank(), axis.title=element_text(size=18,face="bold"),
        legend.text.align = 0, legend.position = "bottom"); size.raw.alt.mel.era 

# Fig. S11 A - raw data shape vs. distance from the Equator
shape.raw.lat.mel.era <- ggplot(era.mel.sex, aes(x=abs(latitude), y=aspect.ratio,colour=species )) +
  geom_point(size=.5, alpha=.4)  +
  #geom_text(label=mel.master$cam.id)+
  #facet_wrap(~sp.subsp)+
  geom_smooth(method='gam',formula =  y ~ poly(x, 1), aes(fill=species))+
  scale_color_manual(name="Species", values = c("melpomene"="#E69F00","erato"="#56B4E9"), labels=c(expression(italic("H. erato")), expression(italic("H. melpomene"))))+
  scale_fill_manual(name="Species", values = c("melpomene"="#E69F00","erato"="#56B4E9"), labels=c(expression(italic("H. erato")), expression(italic("H. melpomene"))))+
  #scale_shape_manual(name="Species", values = c("melpomene"=16,"erato"=17), labels=c(expression(italic("H. erato")), expression(italic("H. melpomene")))) +
  xlab("Distance from Equator (degrees)")+ 
  ylab("Wing aspect ratio") + 
  theme_classic()+ stat_cor(show.legend = FALSE)+
  #facet_wrap(~sex, labeller =labeller(sex=capitalize))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),
        panel.background = element_blank(), legend.text = element_text(size=12), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(), axis.text=element_text(size=12),
        strip.text = element_text(size=14, face="bold"),
        panel.grid.minor = element_blank(), axis.title=element_text(size=18,face="bold"),
        legend.text.align = 0, legend.position = "bottom"); shape.raw.lat.mel.era 

# Fig. S11 B
size.raw.lat.mel.era <- ggplot(era.mel.sex, aes(x=abs(latitude), y=area.mm2,colour=species )) +
  geom_point(size=.5, alpha=.4)  +
  #geom_text(label=mel.master$cam.id)+
  #facet_wrap(~sp.subsp)+
  geom_smooth(method='gam',formula =  y ~ poly(x, 1), aes(fill=species))+
  scale_color_manual(name="Species", values = c("melpomene"="#E69F00","erato"="#56B4E9"), labels=c(expression(italic("H. erato")), expression(italic("H. melpomene"))))+
  scale_fill_manual(name="Species", values = c("melpomene"="#E69F00","erato"="#56B4E9"), labels=c(expression(italic("H. erato")), expression(italic("H. melpomene"))))+
  #scale_shape_manual(name="Species", values = c("melpomene"=16,"erato"=17), labels=c(expression(italic("H. erato")), expression(italic("H. melpomene")))) +
  xlab("Distance from Equator (degrees)")+ 
  ylab(bquote(bold('Wing area ('*mm^2*')'))) + 
  theme_classic()+ stat_cor(show.legend = FALSE)+
  #facet_wrap(~sex, labeller =labeller(sex=capitalize))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),
        panel.background = element_blank(), legend.text = element_text(size=12), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(), axis.text=element_text(size=12),
        strip.text = element_text(size=14, face="bold"),
        panel.grid.minor = element_blank(), axis.title=element_text(size=18,face="bold"),
        legend.text.align = 0, legend.position = "bottom"); size.raw.lat.mel.era 

# Fig. S12 - wing.pattern 
wing.pattern.info <- summarise (group_by(era.mel, species, wing.pattern),
                           n=n())
era.mel <- subset(era.mel, !(is.na(wing.pattern)))
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.08, 1), symbols = c("****", "***", "**", "*", "•", "ns"))

wing.pattern.era.mel.AR <- ggplot(data=era.mel, aes(x=wing.pattern,  y=aspect.ratio, fill=species))+
  xlab("") + ylab("Aspect Ratio")+
  scale_fill_manual(name="Species" , values=c("erato"="#56B4E9", "melpomene"="#E69F00"), breaks=c("erato","melpomene"))+
  geom_point(pch = 21, position = position_jitterdodge(), alpha=.6, size=.6)+
  geom_boxplot(outlier.shape = NA, alpha=.7) +
  theme_bw()+
  stat_compare_means(method = "t.test", label = "p.signif", symnum.args=symnum.args, size=5)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),
        panel.background = element_blank(), legend.text = element_text(size=12), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(), axis.text.x=element_text(size=14,face="bold", colour="black"),
        panel.grid.minor = element_blank(), axis.title=element_text(size=14,face="bold", colour="black"),
        legend.text.align = 0, legend.position = "none",strip.text = element_text(size=14,face="bold")); wing.pattern.era.mel.AR

wing.pattern.era.mel.size  <- ggplot(data=era.mel, aes(x=wing.pattern,  y=area.mm2, fill=species))+
  xlab("") + ylab(bquote(bold('Wing area ('*mm^2*')')))+
  scale_fill_manual(name="Species" , values=c("erato"="#56B4E9", "melpomene"="#E69F00"), breaks=c("erato","melpomene"))+
  geom_point(pch = 21, position = position_jitterdodge(), alpha=.6, size=.6)+
  geom_boxplot(outlier.shape = NA, alpha=.7) +
  theme_bw()+
  stat_compare_means(method = "t.test", label = "p.signif", symnum.args=symnum.args, size=5)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),
        panel.background = element_blank(), legend.text = element_text(size=12), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(), axis.text.x=element_text(size=14,face="bold", colour="black"),
        panel.grid.minor = element_blank(), axis.title=element_text(size=14,face="bold", colour="black"),
        legend.text.align = 0, legend.position = "bottom",strip.text = element_text(size=14,face="bold")); wing.pattern.era.mel.size 

plot_grid(wing.pattern.era.mel.AR , wing.pattern.era.mel.size,ncol = 1, labels=c("A","B"))


############ 2. Analyses ############
#### era shape ####
lm1 <- lm(aspect.ratio ~ area.mm2*altitude+abs(latitude)+longitude+sex, 
          data=era); summary(lm1)
lm0 <- lm(aspect.ratio ~ 1,  data=era); summary(lm0)
shape.final<-stepAIC(lm1, direction = "backward",scope=list(upper=lm1,lower=lm0)); summary(shape.final)
shape.final$anova

modelEffectSizes(shape.final); #plot(era.shape.final)
relimp.shape.final <- calc.relimp.lm(shape.final, type=c("lmg"), rela=TRUE ); print(relimp.shape.final)
relimp.shape.final.df <- data_frame(rel.r2=relimp.shape.final@lmg, variable= rownames(as.data.frame(relimp.shape.final@lmg)))

era.shape.plot <- ggplot(data=relimp.shape.final.df, aes(x=reorder(variable, -rel.r2), y=rel.r2)) +
  geom_bar(fill="grey",colour="grey",stat="identity")+
  theme_classic()+ ylab(bquote(bold('Relative '*R^2*'')))+
  xlab(bquote(bold('Fixed effect')))+
  scale_x_discrete(labels=c("area.mm2:altitude" = "area*altitude", "altitude:abs(latitude)"="altitude*dist.Eq.", 
                            "abs(latitude):longitude"="longitude*dist.Eq.","abs(latitude)"="dist.Eq","area.mm2"="area"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),
        axis.line=element_blank(),
        panel.background = element_blank(), legend.text = element_text(size=12), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(), axis.text=element_text(size=14),
        strip.text = element_text(size=12, face="bold"),plot.margin=grid::unit(c(3,10,3,3), "mm"),
        panel.grid.minor = element_blank(), axis.title=element_text(size=14,face="bold"),
        legend.text.align = 0, legend.position = "bottom"); era.shape.plot


#boot <- boot.relimp(shape.final, b=1000, type = c("lmg"), rank = TRUE, diff = TRUE, rela = TRUE); booteval.relimp(boot) # print result
#png("3.plots/SOM/rela.impo.era.mel/era.shape.rel.imp.boot1000.mod.png", res = 150,width = 2380, height = 980, units = "px")
#plot(booteval.relimp(boot,sort=TRUE),names.abbrev=17,main="Relative importances for wing shape H. erato", cex.title=1.5)
#dev.off()

#### mel shape #####
lm1 <- lm(aspect.ratio ~area.mm2*altitude+abs(latitude)+longitude+sex, 
          data=mel); summary(lm1)
lm0 <- lm(aspect.ratio ~ 1,  data=mel); summary(lm0)
shape.final<-stepAIC(lm1, direction = "both",scope=list(upper=lm1,lower=lm0)); summary(shape.final)
shape.final$anova
#plot(lm)
anova(shape.final)

modelEffectSizes(shape.final); #plot(era.shape.final)
relimp.shape.final <- calc.relimp.lm(shape.final, type=c("lmg"), rela=TRUE ); print(relimp.shape.final)
relimp.shape.final.df <- data_frame(rel.r2=relimp.shape.final@lmg, variable= rownames(as.data.frame(relimp.shape.final@lmg)))

mel.shape.plot <- ggplot(data=relimp.shape.final.df, aes(x=reorder(variable, -rel.r2), y=rel.r2)) +
  geom_bar(fill="grey",colour="grey",stat="identity")+
  theme_classic()+ ylab(bquote(bold('Relative '*R^2*'')))+
  xlab(bquote(bold('Fixed effect')))+
  scale_x_discrete(labels=c("area.mm2:altitude" = "area*altitude", "altitude:abs(latitude)"="altitude*dist.Eq.", 
                            "abs(latitude):longitude"="longitude*dist.Eq.","abs(latitude)"="dist.Eq","area.mm2"="area", 
                            "altitude:longitude"="altitude*longitude"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),
        axis.line=element_blank(),
        panel.background = element_blank(), legend.text = element_text(size=12), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(), axis.text=element_text(size=14),
        strip.text = element_text(size=12, face="bold"),plot.margin=grid::unit(c(3,10,3,3), "mm"),
        panel.grid.minor = element_blank(), axis.title=element_text(size=14,face="bold"),
        legend.text.align = 0, legend.position = "bottom"); mel.shape.plot

#### era size ####
lm1 <- lm(area.mm2 ~ aspect.ratio*altitude+abs(latitude)+longitude+sex, 
          data=era); summary(lm1)
lm0 <- lm(aspect.ratio ~ 1,  data=era); summary(lm0)
size.final <-stepAIC(lm1, direction = "backward",scope=list(upper=lm1,lower=lm0)); summary(size.final)
size.final$anova
anova(size.final)

modelEffectSizes(size.final); #plot(era.size.final)
relimp.size.final <- calc.relimp.lm(size.final, type=c("lmg"), rela=TRUE ); print(relimp.size.final)
relimp.size.final.df <- data_frame(rel.r2=relimp.size.final@lmg, variable= rownames(as.data.frame(relimp.size.final@lmg)))

era.size.plot <- ggplot(data=relimp.size.final.df, aes(x=reorder(variable, -rel.r2), y=rel.r2)) +
  geom_bar(fill="grey",colour="grey",stat="identity")+
  theme_classic()+ ylab(bquote(bold('Relative '*R^2*'')))+
  xlab(bquote(bold('Fixed effect')))+
  scale_x_discrete(labels=c("area.mm2:altitude" = "area*altitude", "altitude:abs(latitude)"="altitude*dist.Eq.", 
                            "abs(latitude):longitude"="longitude*dist.Eq.","abs(latitude)"="dist.Eq","area.mm2"="area", 
                            "altitude:longitude"="altitude*longitude", "aspect.ratio:altitude"="aspect.ratio\n*altitude"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),
        axis.line=element_blank(),
        panel.background = element_blank(), legend.text = element_text(size=12), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(), axis.text=element_text(size=14),
        strip.text = element_text(size=12, face="bold"),plot.margin=grid::unit(c(3,10,3,3), "mm"),
        panel.grid.minor = element_blank(), axis.title=element_text(size=14,face="bold"),
        legend.text.align = 0, legend.position = "bottom");era.size.plot 

#### mel size ####
lm1 <- lm(area.mm2 ~ aspect.ratio*altitude+abs(latitude)+longitude+sex, 
          data=mel); summary(lm1)
lm0 <- lm(aspect.ratio ~ 1,  data=mel); summary(lm0)
size.final <-stepAIC(lm1, direction = "backward",scope=list(upper=lm1,lower=lm0)); summary(size.final)
size.final$anova
anova(size.final)

modelEffectSizes(size.final); #plot(era.size.final)
relimp.size.final <- calc.relimp.lm(size.final, type=c("lmg"), rela=TRUE ); print(relimp.size.final)
relimp.size.final.df <- data_frame(rel.r2=relimp.size.final@lmg, variable= rownames(as.data.frame(relimp.size.final@lmg)))

mel.size.plot <- ggplot(data=relimp.size.final.df, aes(x=reorder(variable, -rel.r2), y=rel.r2)) +
  geom_bar(fill="grey",colour="grey",stat="identity")+
  theme_classic()+ ylab(bquote(bold('Relative '*R^2*'')))+
  xlab(bquote(bold('Fixed effect')))+
  scale_x_discrete(labels=c("area.mm2:altitude" = "area*altitude", "altitude:abs(latitude)"="altitude*dist.Eq.", 
                            "abs(latitude):longitude"="longitude*dist.Eq.","abs(latitude)"="dist.Eq","area.mm2"="area", 
                            "altitude:longitude"="altitude*longitude", "aspect.ratio:altitude"="aspect.ratio\n*altitude"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),
        axis.line=element_blank(),
        panel.background = element_blank(), legend.text = element_text(size=12), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(), axis.text=element_text(size=14),
        strip.text = element_text(size=12, face="bold"),plot.margin=grid::unit(c(3,10,3,3), "mm"),
        panel.grid.minor = element_blank(), axis.title=element_text(size=14,face="bold"),
        legend.text.align = 0, legend.position = "bottom"); mel.size.plot

#### Fig. S13 ####
plot_grid(NULL,era.shape.plot, mel.shape.plot, era.size.plot, mel.size.plot, ncol = 1, scale = .95, label_x = 0,
          hjust = 0, label_size=14, label_y=1.03, align = "v", rel_heights = c(.1,1,1,1,1),
          labels = c("","A) Aspect ratio H. erato", "B) Aspect ratio H. melpomene", "C) Area H. erato", "D) Area H. melpomene"))


