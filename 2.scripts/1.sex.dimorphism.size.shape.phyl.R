##### Altitude and life-history shape the evolution of Heliconius wings ####
#### Gabriela Montejo-Kovacevich 2019 ####
####### wing shape/size sexual dimorphism ######
#### packages #####
rm(list=ls())
dev.off()
library(phylosignal)
library(adephylo)
library(ade4)
library(dplyr)
library(data.table)
library(geiger)
library(ape)
library(nlme)
library(rptR)
library(rr2)
library(phytools)
library(phylobase)
library(phylolm)
library(ggplot2)
library(caper)
library(ggpubr)
library(tidyr)
library(lattice)
library(lme4)
library(leaps)
library(MASS)
library(lmSupport)
library(car)
library(relaimpo)
library(scatterplot3d)
library(plotly)
library(rgl)
library(forcats)
library(cowplot)
library(ggExtra)
library(AICcmodavg)
library(caper)
library(MuMIn)
library(broom)
library(broom.mixed)
library(lmerTest)
library(pwr)
library(multcompView)

####  data  #####
setwd("../heliconius-wings-2019/")
master <- read.csv("1.data/master.analyses.csv")
mean.df <- read.csv("1.data/mean.df.analyses.csv")
nat.hist <- read.csv("1.data/species.nat.hist.csv")
pruned.tree <- read.tree("1.data/pruned.tree.phylo")

####  data handling  #####
master$sp.short <- substr(master$species,0,3)
# subset 
master.sex <- subset(master, sex!="")
rownames(mean.df) <- mean.df$genus.sp
master$species <- factor(master$species, levels=mean.df$species)

####  functions  #####
# stepAICc function from Read, Quentin D. et al. (2018), Data from: Tropical bird species have less variable body sizes, Dryad, Dataset, https://doi.org/10.5061/dryad.gd290
stepAICc <- function (object, scope, scale = 0, direction = c("both", "backward", 
                                                              "forward"), trace = 1, keep = NULL, steps = 1000, use.start = FALSE, 
                      k = 2, ...) 
{
  mydeviance <- function(x, ...) {
    dev <- deviance(x)
    if (!is.null(dev)) 
      dev
    else MuMIn::AICc(x, k=0)
  }
  cut.string <- function(string) {
    if (length(string) > 1L) 
      string[-1L] <- paste("\n", string[-1L], sep = "")
    string
  }
  re.arrange <- function(keep) {
    namr <- names(k1 <- keep[[1L]])
    namc <- names(keep)
    nc <- length(keep)
    nr <- length(k1)
    array(unlist(keep, recursive = FALSE), c(nr, nc), list(namr, 
                                                           namc))
  }
  step.results <- function(models, fit, object, usingCp = FALSE) {
    change <- sapply(models, "[[", "change")
    rd <- sapply(models, "[[", "deviance")
    dd <- c(NA, abs(diff(rd)))
    rdf <- sapply(models, "[[", "df.resid")
    ddf <- c(NA, abs(diff(rdf)))
    AIC <- sapply(models, "[[", "AIC")
    heading <- c("Stepwise Model Path \nAnalysis of Deviance Table", 
                 "\nInitial Model:", deparse(formula(object)), "\nFinal Model:", 
                 deparse(formula(fit)), "\n")
    aod <- if (usingCp) 
      data.frame(Step = change, Df = ddf, Deviance = dd, 
                 `Resid. Df` = rdf, `Resid. Dev` = rd, Cp = AIC, 
                 check.names = FALSE)
    else data.frame(Step = change, Df = ddf, Deviance = dd, 
                    `Resid. Df` = rdf, `Resid. Dev` = rd, AIC = AIC, 
                    check.names = FALSE)
    attr(aod, "heading") <- heading
    class(aod) <- c("Anova", "data.frame")
    fit$anova <- aod
    fit
  }
  Terms <- terms(object)
  object$formula <- Terms
  if (inherits(object, "lme")) 
    object$call$fixed <- Terms
  else if (inherits(object, "gls")) 
    object$call$model <- Terms
  else object$call$formula <- Terms
  if (use.start) 
    warning("'use.start' cannot be used with R's version of 'glm'")
  md <- missing(direction)
  direction <- match.arg(direction)
  backward <- direction == "both" | direction == "backward"
  forward <- direction == "both" | direction == "forward"
  if (missing(scope)) {
    fdrop <- numeric()
    fadd <- attr(Terms, "factors")
    if (md) 
      forward <- FALSE
  }
  else {
    if (is.list(scope)) {
      fdrop <- if (!is.null(fdrop <- scope$lower)) 
        attr(terms(update.formula(object, fdrop)), "factors")
      else numeric()
      fadd <- if (!is.null(fadd <- scope$upper)) 
        attr(terms(update.formula(object, fadd)), "factors")
    }
    else {
      fadd <- if (!is.null(fadd <- scope)) 
        attr(terms(update.formula(object, scope)), "factors")
      fdrop <- numeric()
    }
  }
  models <- vector("list", steps)
  if (!is.null(keep)) 
    keep.list <- vector("list", steps)
  n <- nobs(object, use.fallback = TRUE)
  fit <- object
  bAIC <- extractAIC(fit, scale, k = k, ...)
  edf <- bAIC[1L]
  bAIC <- MuMIn::AICc(fit, k=k)
  if (is.na(bAIC)) 
    stop("AIC is not defined for this model, so 'stepAIC' cannot proceed")
  if (bAIC == -Inf) 
    stop("AIC is -infinity for this model, so 'stepAIC' cannot proceed")
  nm <- 1
  Terms <- terms(fit)
  if (trace) {
    cat("Start:  AIC=", format(round(bAIC, 2)), "\n", cut.string(deparse(formula(fit))), 
        "\n\n", sep = "")
    utils::flush.console()
  }
  models[[nm]] <- list(deviance = mydeviance(fit), df.resid = n - 
                         edf, change = "", AIC = bAIC)
  if (!is.null(keep)) 
    keep.list[[nm]] <- keep(fit, bAIC)
  usingCp <- FALSE
  while (steps > 0) {
    steps <- steps - 1
    AIC <- bAIC
    ffac <- attr(Terms, "factors")
    if (!is.null(sp <- attr(Terms, "specials")) && !is.null(st <- sp$strata)) 
      ffac <- ffac[-st, ]
    scope <- factor.scope(ffac, list(add = fadd, drop = fdrop))
    aod <- NULL
    change <- NULL
    if (backward && length(scope$drop)) {
      aod <- dropterm(fit, scope$drop, scale = scale, 
                      trace = max(0, trace - 1), k = k, ...)
      rn <- row.names(aod)
      row.names(aod) <- c(rn[1L], paste("-", rn[-1L], 
                                        sep = " "))
      if (any(aod$Df == 0, na.rm = TRUE)) {
        zdf <- aod$Df == 0 & !is.na(aod$Df)
        nc <- match(c("Cp", "AIC"), names(aod))
        nc <- nc[!is.na(nc)][1L]
        ch <- abs(aod[zdf, nc] - aod[1, nc]) > 0.01
        if (any(is.finite(ch) & ch)) {
          warning("0 df terms are changing AIC")
          zdf <- zdf[!ch]
        }
        if (length(zdf) > 0L) 
          change <- rev(rownames(aod)[zdf])[1L]
      }
    }
    if (is.null(change)) {
      if (forward && length(scope$add)) {
        aodf <- addterm(fit, scope$add, scale = scale, 
                        trace = max(0, trace - 1), k = k, ...)
        rn <- row.names(aodf)
        row.names(aodf) <- c(rn[1L], paste("+", rn[-1L], 
                                           sep = " "))
        aod <- if (is.null(aod)) 
          aodf
        else rbind(aod, aodf[-1, , drop = FALSE])
      }
      attr(aod, "heading") <- NULL
      if (is.null(aod) || ncol(aod) == 0) 
        break
      nzdf <- if (!is.null(aod$Df)) 
        aod$Df != 0 | is.na(aod$Df)
      aod <- aod[nzdf, ]
      if (is.null(aod) || ncol(aod) == 0) 
        break
      nc <- match(c("Cp", "AIC"), names(aod))
      nc <- nc[!is.na(nc)][1L]
      o <- order(aod[, nc])
      if (trace) {
        print(aod[o, ])
        utils::flush.console()
      }
      if (o[1L] == 1) 
        break
      change <- rownames(aod)[o[1L]]
    }
    usingCp <- match("Cp", names(aod), 0) > 0
    fit <- update(fit, paste("~ .", change), evaluate = FALSE)
    fit <- eval.parent(fit)
    nnew <- nobs(fit, use.fallback = TRUE)
    if (all(is.finite(c(n, nnew))) && nnew != n) 
      stop("number of rows in use has changed: remove missing values?")
    Terms <- terms(fit)
    bAIC <- extractAIC(fit, scale, k = k, ...)
    edf <- bAIC[1L]
    bAIC <- MuMIn::AICc(fit, k=k)
    if (trace) {
      cat("\nStep:  AIC=", format(round(bAIC, 2)), "\n", 
          cut.string(deparse(formula(fit))), "\n\n", sep = "")
      utils::flush.console()
    }
    if (bAIC >= AIC + 1e-07) 
      break
    nm <- nm + 1
    models[[nm]] <- list(deviance = mydeviance(fit), df.resid = n - 
                           edf, change = change, AIC = bAIC)
    if (!is.null(keep)) 
      keep.list[[nm]] <- keep(fit, bAIC)
  }
  if (!is.null(keep)) 
    fit$keep <- re.arrange(keep.list[seq(nm)])
  step.results(models = models[seq(nm)], fit, object, usingCp)
}

# function to create tukey post-hoc groups for plotting
generate_label_df <- function(TUKEY, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$treatment=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
  return(Tukey.labels)
}
#########################   1. SEXUAL DIMORPHISM, Fig. 2 ######
######## WING SIZE SD ####
# prep phyl data
rownames(mean.df) <- mean.df$genus.sp
p4d <- phylo4d(pruned.tree, mean.df)

# barplot prep
tip.col <- rep(1, nTips(p4d)); tip.col <- c("black", "black", "black", "black", "#CC79A7","#CC79A7", "#CC79A7", "#CC79A7", 
                                            "#CC79A7", "#CC79A7", "black", "black", "black")
# create colors
mat.col <- ifelse(tdata(p4d, type="tip", label.type="row.names") =="era.clade", "#56B4E9", "#E69F00")

# create tip labels
tip.lab<- list() ; for (i in 1:13){
  a <- strsplit(paste(tipLabels(p4d)),split='_')[[i]][2]
  print(tip.lab[[i]] <- a)
} ; tip.lab <- paste("H. ", tip.lab, sep = "")


#### fig 2 A and B - boxplot raw data #### 
# specify symbols for significance cutoffs
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.095, 1), symbols = c("****", "***", "**", "*", "•", "ns"))
# order
master.sex$sex <- factor(master.sex$sex, levels =c("male", "female"))
master.sex$sp.short<- factor(master.sex$sp.short, levels=mean.df$sp.short)
levels.soli <- c('tel','cly','era','tim','cyd','mel','num')                         
levels.greg <- c('ele','sar','xan','hie','dor','wal')

# plot fig. 2 A and B
soli <- ggplot(data=subset(master.sex, larva=="solitary"), aes(x=factor(sp.short, levels=levels.soli),  y=area.mm2, fill=sex))+
  xlab("") + ylab("")+
  scale_fill_manual(name="Sex" , values=c("male"="grey", "female"="white"), breaks=c("male","female"), labels=c("Male", "Female"))+
  scale_y_continuous(breaks=c(200,400,600, 800))+
  geom_point(pch = 21, position = position_jitterdodge(), alpha=.6, size=.6)+
  geom_boxplot(outlier.shape = NA, alpha=.7) +
  theme_bw()+
  #facet_wrap(~larva, nrow = 2, scales = "free_x")
  stat_compare_means(method = "t.test", label = "p.signif", symnum.args=symnum.args, size=5,label.y=820)+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),
          panel.background = element_blank(), legend.text = element_text(size=12), legend.title =  element_text(size=12),
          panel.grid.major = element_blank(), axis.text.x=element_text(size=14,face="bold", colour="black"),
          panel.grid.minor = element_blank(), axis.title=element_text(size=14,face="bold", colour="black"),
          legend.text.align = 0, legend.position = "right",strip.text = element_text(size=14,face="bold")); soli

 greg <- ggplot(data=subset(master.sex, larva=="gregarious"), aes(x=factor(sp.short, levels=levels.greg),  y=area.mm2, fill=sex))+
  xlab("") + ylab("")+
  scale_fill_manual(name="Sex" , values=c("male"="grey", "female"="white"), breaks=c("male","female"), labels=c("Male", "Female"))+
  geom_point(pch = 21, position = position_jitterdodge(), alpha=.6, size=.6)+
  scale_y_continuous(breaks=c(200,400,600, 800))+
  geom_boxplot(outlier.shape = NA, alpha=.7) +
  #facet_wrap(~larva, nrow = 2, scales = "free_x")
  stat_compare_means(method = "t.test", label = "p.signif", symnum.args=symnum.args, size=5,label.y=700)+
  theme_bw()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),
        panel.background = element_blank(), legend.text = element_text(size=12), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(), axis.text.x=element_text(size=14,face="bold",colour ="#CC79A7" ),
        panel.grid.minor = element_blank(), axis.title=element_text(size=14,face="bold"),
        legend.text.align = 0, legend.position = "right",strip.text = element_text(size=14,face="bold")); greg

# combine
plot_grid(soli, greg,ncol = 1,scale=.95)+ 
draw_label(bquote(bold('Wing area ('*mm^2*')')), x=0, y=  .5, vjust= 2,  angle= 90)

#### fig 2 C - barplot sexual size dimorphism  #### 
dev.off(); plot.new(); barplot.phylo4d(p4d, trait = c(27), bar.lwd = 15, scale = FALSE,trait.font = 2, show.data.axis = TRUE, trait.cex = 1.3,
                                       tree.xlim = c(-.5,11),tree.ratio = 0.2, tip.labels = tip.lab, tree.ladderize = FALSE, center = FALSE, data.xlim = c(-11,11),
                                       tip.cex = 1.8,  tip.col = tip.col, grid.vertical = FALSE, trait.bg.col = FALSE, tip.font = 3,
                                       #error.bar.sup = mat.e, error.bar.inf = mat.e, error.bar.col = "bisque4",
                                       show.box = TRUE,trait.labels = c("Sexual Size Dimorphism"), 
                                       bar.col = c("black", "black", "black", "black", "#CC79A7","#CC79A7", "#CC79A7", "#CC79A7", "#CC79A7", "#CC79A7", "black", "black", "black"))

fig2b <- recordPlot() ; fig2b
plot_grid(fig2b, nrow=1,labels = c(''), scale = .8)

#### t-tests size sexual dimorphism  ####
sp.list <-unique(master.sex$genus.sp)
t.test.df <- list()
t.test.master <- data_frame(species=sp.list, p.value=c(0), mean.female=c(0), 
                            mean.male=c(0), t.value=c(0), df=c(0)) 
# loop through species and keep t-test values
options(scipen = 999)
for (i in 1:length(sp.list)){
  temp <- master.sex[master.sex$genus.sp==sp.list[i],]
  t.test.df[[i]] <-t.test(area.mm2 ~ sex, data=temp)
  t.test.master[i,]$p.value <-t.test.df[[i]]$p.value
  t.test.master[i,]$mean.female <-t.test.df[[i]]$estimate[1]
  t.test.master[i,]$mean.male <-t.test.df[[i]]$estimate[2]
  t.test.master[i,]$t.value <-t.test.df[[i]]$statistic
  t.test.master[i,]$df <-t.test.df[[i]]$parameter[1]
}
t.test.master <- as.data.frame(t.test.master); row.names(t.test.master)<- t.test.master$species
t.test.master <- t.test.master[(rev(pruned.tree$tip.label)),]; t.test.master$species = factor(t.test.master$species)
# write.csv(t.test.master, "1.data/t.test.sisd.csv", row.names = FALSE)


######## WING SHAPE SD ####
#### fig S4 A and B - boxplot raw data #### 
names(master)
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.08, 1), symbols = c("****", "***", "**", "*", "•", "ns"))
master.sex$sex <- factor(master.sex$sex, levels =c("male", "female"))
master.sex$sp.short<- factor(master.sex$sp.short, levels=mean.df$sp.short)

era <- ggplot(data=subset(master.sex, clade=="era.clade"), aes(x=sp.short,  y=aspect.ratio, fill=sex))+
  xlab("") + ylab("Aspect Ratio")+
  scale_fill_manual(name="Sex" , values=c("male"="grey", "female"="white"), breaks=c("male","female"), labels=c("Male", "Female"))+
  geom_point(pch = 21, position = position_jitterdodge(), alpha=.6, size=.6)+
  geom_boxplot(outlier.shape = NA, alpha=.7) +
  theme_bw()+
  #facet_wrap(~larva, nrow = 2, scales = "free_x")
  stat_compare_means(method = "t.test", label = "p.signif", symnum.args=symnum.args, size=5)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),
        panel.background = element_blank(), legend.text = element_text(size=12), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(), axis.text.x=element_text(size=14,face="bold", colour="#56B4E9"),
        panel.grid.minor = element_blank(), axis.title=element_text(size=14,face="bold", colour="black"),
        legend.text.align = 0, legend.position = "right",strip.text = element_text(size=14,face="bold")); era

mel <- ggplot(data=subset(master.sex, clade=="mel.clade"), aes(x=sp.short,  y=aspect.ratio, fill=sex))+
  xlab("") + ylab("Aspect Ratio")+
  scale_fill_manual(name="Sex" , values=c("male"="grey", "female"="white"), breaks=c("male","female"), labels=c("Male", "Female"))+
  geom_point(pch = 21, position = position_jitterdodge(), alpha=.6, size=.6)+
  geom_boxplot(outlier.shape = NA, alpha=.7) +
  #facet_wrap(~larva, nrow = 2, scales = "free_x")
  stat_compare_means(method = "t.test", label = "p.signif", symnum.args=symnum.args, size=5)+
  theme_bw()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),
        panel.background = element_blank(), legend.text = element_text(size=12), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(), axis.text.x=element_text(size=14,face="bold",colour ="#E69F00" ),
        panel.grid.minor = element_blank(), axis.title=element_text(size=14,face="bold"),
        legend.text.align = 0, legend.position = "right",strip.text = element_text(size=14,face="bold")); mel

plot_grid(era, mel,ncol = 1,scale=1)

#### t-tests shsd ####
sp.list <-unique(master.sex$genus.sp)
t.test.df <- list()
t.test.master <- data_frame(species=sp.list, p.value=c(0), mean.female=c(0), 
                            mean.male=c(0), t.value=c(0), df=c(0)) 
options(scipen = 999)
for (i in 1:length(sp.list)){
  temp <- master.sex[master.sex$genus.sp==sp.list[i],]
  t.test.df[[i]] <-t.test(aspect.ratio ~ sex, data=temp)
  t.test.master[i,]$p.value <-t.test.df[[i]]$p.value
  t.test.master[i,]$mean.female <-t.test.df[[i]]$estimate[1]
  t.test.master[i,]$mean.male <-t.test.df[[i]]$estimate[2]
  t.test.master[i,]$t.value <-t.test.df[[i]]$statistic
  t.test.master[i,]$df <-t.test.df[[i]]$parameter[1]
}
t.test.master <- as.data.frame(t.test.master); row.names(t.test.master)<- t.test.master$species
t.test.master <- t.test.master[(rev(pruned.tree$tip.label)),]; t.test.master$species = factor(t.test.master$species)

write.csv(t.test.master, "1.data/t.test.shsd.csv", row.names = FALSE)



######## SD MODELS ######
#### SIZE #### 
# 2 ways 
# LM stepAIC
hist(mean.df$sisd.raw)
fit1 <- lm(sisd.raw ~  larva+shsd.raw+ size.mean+shape.mean+clade , mean.df); summary(fit1)
fit2 <- lm(sisd.raw ~  1, mean.df); summary(fit2)
sisd.lm <-stepAICc(fit1,direction="backward",scope=list(upper=fit1,lower=fit2)); summary(sisd.lm ); sisd.lm$anova 
summary(sisd.lm)

boot <- boot.relimp(sisd.lm, b=200, type = c("lmg"), rank = TRUE, diff = TRUE, rela = TRUE); booteval.relimp(boot) # print result
png("3.plots/SOM/SISD.rel.imp.boot1000.mod.png", res = 150,width = 980, height = 980, units = "px")
plot(booteval.relimp(boot,sort=TRUE),names.abbrev=17,main="Relative importances for wing SiSD", cex.title=1.5)
dev.off()

modelEffectSizes(sisd.lm); #plot(sisd.lm)
relimp.lm.sisd <- calc.relimp.lm(sisd.lm, type=c("lmg"), rela=TRUE ); print(relimp.lm.sisd) # of the variation explained by the model, 34% is latitude
relimp.lm.sisd.df <- data_frame(rel.r2=relimp.lm.sisd@lmg, variable= rownames(as.data.frame(relimp.lm.sisd@lmg)))
ggplot(data=relimp.lm.sisd.df, aes(x=variable, y=rel.r2)) +
  geom_bar(stat="identity")
anova(sisd.lm)

# PGLS (for SI)
bm.tree<-corBrownian(phy=pruned.tree)
vf<-varComb(varFixed(~size.n),varFixed(~1/sqrt(size.se)),
            varFixed(~1/sqrt(alt.se)))

# stepwise model selection AIC
fit1 <- gls(sisd.raw ~  larva+shsd.raw+ size.mean+shape.mean+clade, mean.df, weights = vf, correlation = bm.tree,method = "ML"); summary(fit1)
fit2 <- gls(sisd.raw ~ 1, mean.df, weights = vf, method = "ML", correlation = bm.tree); summary(fit2)
step.size <-stepAICc(fit1, direction = "backward",scope=list(upper=fit1,lower=fit2)); summary(step.size)
step.size$anova

#### SHAPE #### 
hist(mean.df$shsd.raw)
fit1 <- lm(shsd.raw ~  larva+shape.mean+sisd.raw+size.mean +clade, mean.df); summary(fit1)
fit2 <- lm(shsd.raw ~  1, mean.df); summary(fit2)
shsd.lm <-stepAICc(fit1,direction="backward",scope=list(upper=fit1,lower=fit2)); summary(step.shape ); shsd.lm$anova 
summary(shsd.lm)

modelEffectSizes(shsd.lm); #plot(sisd.lm)
relimp.lm.shsd <- calc.relimp.lm(shsd.lm, type=c("lmg"), rela=TRUE ); print(relimp.lm.shsd) # of the variation explained by the model, 34% is latitude
relimp.lm.shsd.df <- data_frame(rel.r2=relimp.lm.shsd@lmg, variable= rownames(as.data.frame(relimp.lm.shsd@lmg)))
ggplot(data=relimp.lm.shsd.df, aes(x=variable, y=rel.r2)) +
  geom_bar(stat="identity")
boot <- boot.relimp(shsd.lm, b = 100, type = c("lmg"), rank = TRUE, diff = TRUE, rela = TRUE); booteval.relimp(boot) # print result
png("3.plots/SOM/shsD.rel.imp.boot1000.mod.png", res = 150,width = 980, height = 980, units = "px")
plot(booteval.relimp(boot,sort=TRUE),names.abbrev=17,main="Relative importances for wing size era", cex.title=1.5)
dev.off()

#########################   2. WING SIZE PHYL SIGNAL ######
#### repeatability ######
# anova
shape.aov<-aov(master$area.mm2 ~ master$species, data = master); summary(shape.aov)

# lmm method, 0.48
rep.shape<- rpt(area.mm2 ~ 1 +(1|genus.sp), grname=c("genus.sp"), datatype="Gaussian", nboot = 100, data=master); print(rep.shape) ; plot(rep.shape, cex.main = 1) # shape repeatability R=0.74 SE = 0.069 P=0

#### phyl signal ######
rownames(mean.df) <- mean.df$genus.sp
p4d <- phylo4d(pruned.tree, mean.df)

# global estimate, of residuals too
names(p4d@data)
phylo.glob <- phyloSignal(p4d = p4d[,c(2,7,12,27,28,36,35,30,29,39,40)], method = "all"); phylo.glob # mostly not, low coeff (high p), median yes K-marginally Cmean

#### Fig S3 - abouheif plots ####
abouheif.test <- abouheif.moran(p4d[,c(2,7,27,28)],method="oriAbouheif")
abouheif.test
# You can see the observed value of Moran's I for each variable and its p value. 
# some traits suggest strong phylogenetic signal effects.
dev.off(); plot.new(); plot(abouheif.test)
figS1.phyl.plot <- recordPlot() ; figS1.phyl.plot 
plot_grid(figS1.phyl.plot); 

#### Fig 3 - barplot phyl signal male size mean, trait2 ######
p4d <- phylo4d(pruned.tree, mean.df)
mat.col <- ifelse(tdata(p4d, type="tip", label.type="row.names") =="era.clade", "#56B4E9", "#E69F00")

tip.lab<- list() ; for (i in 1:13){
  a <- strsplit(paste(tipLabels(p4d)),split='_')[[i]][2]
  print(tip.lab[[i]] <- a)} ; tip.lab <- paste("H. ", tip.lab, sep = "")

# obtain real error bars for size
se.df <- matrix(p4d@data$size.se)
se.df <- matrix(p4d@data$size.male.se)
# col name must be trait name, not se
rownames(se.df) <- p4d@data$genus.sp; colnames(se.df) <- c("size.mean")
rownames(se.df) <- p4d@data$genus.sp; colnames(se.df) <- c("size.male.mean")
# plot
names(p4d@data)
dev.off(); plot.new(); barplot.phylo4d(p4d, trait = c(35), bar.lwd = 15, scale = FALSE,trait.font = 2, show.data.axis = TRUE, trait.cex = 1.3,
                                       tree.xlim = c(-.5,11), tree.ratio = 0.2, tip.labels = tip.lab, tree.ladderize = FALSE, center = TRUE, data.xlim = c(-150,110),
                                       tip.cex = 1.8, error.bar.sup = se.df, error.bar.inf = se.df, tip.col = c(rep("#E69F00",8),rep("#56B4E9",5)), grid.vertical = FALSE, trait.bg.col = FALSE, tip.font = 3,
                                       show.box = TRUE,trait.labels = c("Wing Size"), 
                                       bar.col = c(rep("#E69F00",8),rep("#56B4E9",5)))

#focusTraits(); text(x=c(-1.8,1.5),y=c(1.5,1.5),labels = c(paste(sprintf('\u2190'),"\nsmaller") , paste(sprintf('\u2192'),"\nbigger")), cex = 1.8); focusStop()
fig3.size.clade <- recordPlot() ; fig3.size.clade


#### Fig. S6A - barplot phyl signal size res, trait40 ######
#  residuals dont have SE
# plot
names(p4d@data)
dev.off(); plot.new(); barplot.phylo4d(p4d, trait = c(40), bar.lwd = 15, scale = FALSE,trait.font = 2, show.data.axis = TRUE, trait.cex = 1.3,
                                       tree.xlim = c(-.5,11), tree.ratio = 0.2, tip.labels = tip.lab, tree.ladderize = FALSE, center = TRUE,
                                       tip.cex = 1.8, tip.col = c(rep("#E69F00",8),rep("#56B4E9",5)), grid.vertical = FALSE, trait.bg.col = FALSE, tip.font = 3,
                                       show.box = TRUE,trait.labels = c("Size residuals adj."), 
                                       bar.col = c(rep("#E69F00",8),rep("#56B4E9",5)))
size.res.clade <- recordPlot() ; size.res.clade

#### Fig. S6B and 6C - correologram size mean, size res #####
options(scipen=999, digits = 2)
dev.off(); plot.new(); size.corr<- phyloCorrelogram(p4d, trait ="size.mean"); plot.phylocorrelogram(size.corr, main = "Phylogenetic correlogram\n size mean"); size.corr <- recordPlot() # signal in size.res
dev.off(); plot.new(); size.res.corr<-  phyloCorrelogram(p4d, trait ="size.res"); plot.phylocorrelogram(size.res.corr, main = "Phylogenetic correlogram\n size residuals");  size.res.corr <- recordPlot()  # signal in size.res

dev.off();
dev.off();
size.corr
size.res.corr

#### Fig. S7 - morans I plots (size and shape) ####
# size
dev.off(); plot.new()
local.i <- lipaMoran(p4d, trait = "size.mean", prox.phylo = "nNodes", as.p4d = TRUE)
points.col <- lipaMoran(p4d, trait = "size.mean",prox.phylo = "nNodes")$p.value
points.col <- ifelse(points.col < 0.05, "red", "black")
dotplot.phylo4d(local.i, dot.col = points.col,  trait.labels = c("Size mean Moran's I"), tip.cex = 1.5,
                grid.vertical = FALSE, trait.bg.col = FALSE,show.box = TRUE, tip.labels = tip.lab)
points.col <- data.frame(points.col)
mean.df$size.mean.moransi <- ifelse((points.col$size.mean[match(mean.df$genus.sp, rownames(points.col))])=="red","sig", "ns")
size.mean.morans.plot <- recordPlot() ; size.mean.morans.plot

# shape
dev.off(); plot.new()
local.i <- lipaMoran(p4d, trait = "shape.mean", prox.phylo = "nNodes", as.p4d = TRUE)
points.col <- lipaMoran(p4d, trait = "shape.mean",prox.phylo = "nNodes")$p.value
points.col <- ifelse(points.col < 0.05, "red", "black")
dotplot.phylo4d(local.i, dot.col = points.col,  trait.labels = c("Shape mean Moran's I"), tip.cex = 1.5,
                grid.vertical = FALSE, trait.bg.col = FALSE,show.box = TRUE, tip.labels = tip.lab)
points.col <- data.frame(points.col)
mean.df$shape.mean.moransi <- ifelse((points.col$shape.mean[match(mean.df$genus.sp, rownames(points.col))])=="red","sig", "ns")
shape.mean.morans.plot <- recordPlot() ; shape.mean.morans.plot

dev.off()
size.shape <- plot_grid(size.mean.morans.plot,NULL,shape.mean.morans.plot, nrow = 1); size.shape

#########################   3. WING SHAPE PHYL SIGNAL ######
#### repeatability ######
# anova
shape.aov<-aov(aspect.ratio ~ species, data = master); summary(shape.aov)

# lmm method, 0.74, 0.09
rep.shape<- rpt(aspect.ratio ~ 1 +(1|genus.sp), grname=c("genus.sp"), datatype="Gaussian", nboot = 100, data=master); print(rep.shape) ; plot(rep.shape, cex.main = 1) # shape repeatability R=0.74 SE = 0.069 P=0

#### phyl signal ######
p4d <- phylo4d(pruned.tree, mean.df)

mat.e <- matrix(abs(rnorm(length(tipLabels(p4d)) * ncol(mean.df), 0, 0.5)), ncol = ncol(mean.df),
                dimnames = list(tipLabels(p4d), names(tdata(p4d))))

# global estimate, of residuals too
names(p4d@data)
p4d@data
phylo.glob <- phyloSignal(p4d = p4d[,c(2,7,12,27,28,34,35)], method = "all"); phylo.glob # mostly not, low coeff (high p), median yes K-marginally Cmean

#### barplot phyl signal shape mean, 7 ######
p4d <- phylo4d(pruned.tree, mean.df)

mat.e <- matrix(abs(rnorm(length(tipLabels(p4d)) * ncol(mean.df), 0, 0.5)), ncol = ncol(mean.df),
                dimnames = list(tipLabels(p4d), names(tdata(p4d))))

# 
tip.col <- rep(1, nTips(p4d)); tip.col <- c(rep("black",2),"#E69F00" , rep("black", 7), "#56B4E9", rep("black",2))
tip.col <- rep(1, nTips(p4d)); tip.col <- c("black", "black", "black", "black", "#CC79A7","#CC79A7", "#CC79A7", "#CC79A7", 
                                            "#CC79A7", "#CC79A7", "black", "black", "black")
mat.col <- ifelse(tdata(p4d, type="tip", label.type="row.names") =="era.clade", "#56B4E9", "#E69F00")

tip.lab<- list() ; for (i in 1:13){
  a <- strsplit(paste(rownames(mat.e)),split='_')[[i]][2]
  print(tip.lab[[i]] <- a)
} ; tip.lab <- paste("H. ", tip.lab, sep = "")
mat.e <- as.data.frame(mat.e)
names(mean.df) # all traits  trait = c(2,7,27,28,29,30)

# obtain real error bars for shape
se.df <- matrix(p4d@data$shape.se)
# col name must be trait name, not se
rownames(se.df) <- p4d@data$genus.sp; colnames(se.df) <- c("shape.mean")

# plot
names(p4d@data)
tip.col <- rep(1, nTips(p4d)); tip.col <- c(rep("black",2),"#E69F00" , rep("black", 7), "#56B4E9", rep("black",2))
dev.off(); plot.new(); barplot.phylo4d(p4d, trait = c(7), bar.lwd = 15, scale = TRUE,trait.font = 2, show.data.axis = TRUE, trait.cex = 1.3,
                                       tree.xlim = c(-.5,11),tree.ratio = 0.2, tip.labels = tip.lab, tree.ladderize = FALSE, 
                                       tip.cex = 1.8, error.bar.sup = se.df, error.bar.inf = se.df, tip.col = tip.col, grid.vertical = FALSE, trait.bg.col = FALSE, tip.font = 3,
                                       data.xlim = c(-2.7,2.3), show.box = TRUE, trait.labels = c("Mean Shape adj."), bar.col = c(rep("#E69F00",8),rep("#56B4E9",5)))
focusTraits(); text(x=c(-1.8,1.5),y=c(1.5,1.5),labels = c(paste(sprintf('\u2190'),"\nrounder") , paste(sprintf('\u2192'),"\nlonger")), cex = 1.8); focusStop()
shape.clade <- recordPlot() ; shape.clade


#### Fig. S5A -barplot phyl signal shape  res,  39 ######
#  residuals dont have SE
# plot
tip.col <- rep(1, nTips(p4d)); tip.col <- c(rep("black",2),"#E69F00" , rep("black", 7), "#56B4E9", rep("black",2))
dev.off(); plot.new(); barplot.phylo4d(p4d, trait = c(39), bar.lwd = 15, scale = TRUE,trait.font = 2, show.data.axis = TRUE, trait.cex = 1.3,
                                       tree.xlim = c(-.5,11),tree.ratio = 0.2, tip.labels = tip.lab, tree.ladderize = FALSE, 
                                       tip.cex = 1.8, tip.col = tip.col, grid.vertical = FALSE, trait.bg.col = FALSE, tip.font = 3,
                                       data.xlim = c(-2.7,2.3), show.box = TRUE, trait.labels = c("Shape residuals adj."), bar.col = c(rep("#E69F00",8),rep("#56B4E9",5)))
focusTraits(); text(x=c(-1.8,1.5),y=c(1.5,1.5),labels = c(paste(sprintf('\u2190'),"\nrounder") , paste(sprintf('\u2192'),"\nlonger")), cex = 1.8); focusStop()
shape.res.clade <- recordPlot() ; shape.res.clade


#### Fig. S5B and 5C - correologram shape mean, shape res #####
options(scipen=999, digits = 2)
dev.off(); plot.new(); shape.corr<- phyloCorrelogram(p4d, trait ="shape.mean"); plot.phylocorrelogram(shape.corr, main = "Phylogenetic correlogram\n shape mean"); shape.corr <- recordPlot() # signal in shape.res
dev.off(); plot.new(); shape.res.corr<-  phyloCorrelogram(p4d, trait ="shape.res"); plot.phylocorrelogram(shape.res.corr, main = "Phylogenetic correlogram\n shape residuals");  shape.res.corr <- recordPlot()  # signal in shape.res

#########################   4. repeatability within sexes ######
# anova
shape.aov<-aov(aspect.ratio ~ species, data = master); summary(shape.aov)
size.aov<-aov(area.mm2 ~ species, data = master); summary(size.aov)
alt.aov<-aov(altitude ~ species, data = master); summary(alt.aov)

# shape a bit more repeatable in males
# male 0.79, 0.05, 0
rep.shape<- rpt(aspect.ratio ~ 1 +(1|genus.sp), grname=c("genus.sp"), datatype="Gaussian", nboot = 100, data=subset(master, master$sex=="male")); print(rep.shape) ; plot(rep.shape, cex.main = 1) # shape repeatability R=0.74 SE = 0.069 P=0
#female , 0.74, 0.06, 0
rep.shape<- rpt(aspect.ratio ~ 1 +(1|genus.sp), grname=c("genus.sp"), datatype="Gaussian", nboot = 100, data=subset(master, master$sex=="female")); print(rep.shape) ; plot(rep.shape, cex.main = 1) # shape repeatability R=0.74 SE = 0.069 P=0

# size more repeatable in males
# male 0.51, 0.11, 0
rep.size<- rpt(area.mm2 ~ 1 +(1|genus.sp), grname=c("genus.sp"), datatype="Gaussian", nboot = 100, data=subset(master, master$sex=="male")); print(rep.size) ; plot(rep.size, cex.main = 1) # size repeatability R=0.74 SE = 0.069 P=0
# female, 0.43
rep.size<- rpt(area.mm2 ~ 1 +(1|genus.sp), grname=c("genus.sp"), datatype="Gaussian", nboot = 100, data=subset(master, master$sex=="female")); print(rep.size) ; plot(rep.size, cex.main = 1) # size repeatability R=0.74 SE = 0.069 P=0

# varance is lower in females, but is it just because of N?
rep.shape<- rpt(area.mm2 ~ 1 +(1|genus.sp), grname=c("genus.sp"), datatype="Gaussian", nboot = 100, data=subset(master, master$sex=="male"), ratio = FALSE); print(rep.shape) ; plot(rep.shape, cex.main = 1) # shape repeatability R=0.74 SE = 0.069 P=0
rep.shape<- rpt(area.mm2 ~ 1 +(1|genus.sp), grname=c("genus.sp"), datatype="Gaussian", nboot = 100, data=subset(master, master$sex=="female"), ratio = FALSE); print(rep.shape) ; plot(rep.shape, cex.main = 1) # shape repeatability R=0.74 SE = 0.069 P=0
rep.shape<- rpt(aspect.ratio ~ 1 +(1|genus.sp), grname=c("genus.sp"), datatype="Gaussian", nboot = 100, data=subset(master, master$sex=="male"), ratio = FALSE); print(rep.shape) ; plot(rep.shape, cex.main = 1) # shape repeatability R=0.74 SE = 0.069 P=0
rep.shape<- rpt(aspect.ratio ~ 1 +(1|genus.sp), grname=c("genus.sp","Residual"), datatype="Gaussian", nboot = 100, data=subset(master, master$sex=="female"), ratio = FALSE); print(rep.shape) ; plot(rep.shape, cex.main = 1) # shape repeatability R=0.74 SE = 0.069 P=0


#########################   5. compare species sizes and shapes ###################
#### Fig. S2 wing size and shape variation across species ######
size.aov<-aov(master$area.mm2 ~ master$species, data = master); summary(size.aov)
# tukey groupings for species
# size
TUKEY <- TukeyHSD(x=size.aov, 'master$species', conf.level=0.95);TUKEY
labels<- generate_label_df(TUKEY, "master$species")
names(labels)<-c('Letters', 'species')
yvalue<-summarise(group_by(master, species),
                  mean=mean(area.mm2))
final<-merge(labels,yvalue) 

wing.size.sp <- ggplot(data=master, aes(x=species,  y=area.mm2, fill=species))+
  xlab("") + ylab("Wing area")+
  scale_fill_manual(name="Species" , values=c("#56B4E9","#56B4E9","#56B4E9","#56B4E9","#56B4E9", "#E69F00", "#E69F00", "#E69F00", "#E69F00", "#E69F00", "#E69F00", "#E69F00", "#E69F00"))+
  geom_point(pch = 21, position = position_jitterdodge(), alpha=.6, size=.6)+
  geom_boxplot(outlier.shape = NA, alpha=.7) +
  theme_bw()+
  geom_text(data = final, aes(x = species, y = mean, label = Letters),vjust=-8,hjust=-.8) +
  stat_compare_means(method = "anova",size=5)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),
        panel.background = element_blank(), legend.text = element_text(size=12), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(), axis.text.x=element_text(size=10,face="bold", colour="black",angle = 25, hjust = 1),
        panel.grid.minor = element_blank(), axis.title=element_text(size=14,face="bold", colour="black"),
        legend.text.align = 0, legend.position = "none",strip.text = element_text(size=14,face="bold")); wing.size.sp

# shape
shape.aov<-aov(master$aspect.ratio ~ master$species, data = master); summary(shape.aov)
TUKEY <- TukeyHSD(x=shape.aov, 'master$species', conf.level=0.95);TUKEY
labels<- generate_label_df(TUKEY, "master$species")
names(labels)<-c('Letters', 'species')
yvalue<-summarise(group_by(master, species),
                  mean=mean(aspect.ratio))
final<-merge(labels,yvalue) 

wing.shape.sp <- ggplot(data=master, aes(x=species,  y=aspect.ratio, fill=species))+
  xlab("") + ylab("Wing aspect ratio")+
  scale_fill_manual(name="Species" , values=c("#56B4E9","#56B4E9","#56B4E9","#56B4E9","#56B4E9", "#E69F00", "#E69F00", "#E69F00", "#E69F00", "#E69F00", "#E69F00", "#E69F00", "#E69F00"))+
  geom_point(pch = 21, position = position_jitterdodge(), alpha=.6, size=.6)+
  geom_boxplot(outlier.shape = NA, alpha=.7) +
  theme_bw()+
  geom_text(data = final, aes(x = species, y = mean, label = Letters),vjust=-8,hjust=-.8) +
  stat_compare_means(method = "anova",size=5)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),
        panel.background = element_blank(), legend.text = element_text(size=12), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(), axis.text.x=element_text(size=12,face="bold", colour="black",angle = 25, hjust = 1),
        panel.grid.minor = element_blank(), axis.title=element_text(size=14,face="bold", colour="black"),
        legend.text.align = 0, legend.position = "none",strip.text = element_text(size=14,face="bold")); wing.shape.sp

plot_grid(wing.size.sp  , wing.shape.sp ,ncol = 1, labels=c("A","B"))

#### Fig. S8 - clade Wing area variation with altitude ###### 
ggplot(data=subset(master, altitude>00), aes(x=altitude, y=area.mm2, fill=clade, color=clade))+
  geom_point(alpha=.5) + stat_cor()+
  #facet_wrap(~species+clade) +
  ylab("Wing area")+
  scale_fill_manual(values=c("era.clade"="#56B4E9", "mel.clade"="#E69F00"))+
  scale_color_manual(values=c("era.clade"="#56B4E9", "mel.clade"="#E69F00"))+
  geom_smooth(method="lm")+
  stat_compare_means(method = "anova",size=5)+
  theme_bw()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),
        panel.background = element_blank(), legend.text = element_text(size=12), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(), axis.text.x=element_text(size=10,face="bold", colour="black"),
        panel.grid.minor = element_blank(), axis.title=element_text(size=14,face="bold", colour="black"),
        legend.text.align = 0, legend.position = "none",strip.text = element_text(size=14,face="bold"))

