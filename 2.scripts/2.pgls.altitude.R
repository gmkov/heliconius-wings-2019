##### Altitude and life-history shape the evolution of Heliconius wings ####
#### Gabriela Montejo-Kovacevich 2019 ####
####### Heliconius wing shape/size PGLS altitude models ######
#### packages #####
rm(list=ls())
dev.off()
library(phylosignal)
library(adephylo)
library(ade4)
library(dplyr)
library(data.table)
library(geiger)o
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

#### data  #####
setwd("../heliconius-wings-2019/")
master <- read.csv("1.data/master.analyses.csv")
mean.df <- read.csv("1.data/mean.df.analyses.csv")
nat.hist <- read.csv("1.data/species.nat.hist.csv")
pruned.tree <- read.tree("1.data/pruned.tree.phylo")

#### data handling  #####
master$sp.short <- substr(master$species,0,3)
master.sex <- subset(master, sex!="")
rownames(mean.df) <- mean.df$genus.sp

#### checks, functions #####
hist(mean.df$size.mean)
hist(mean.df$shape.mean)
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

######################### 1. WING SIZE PGLS #############
###### PGLS ######  
# make trees
bm.tree<-corBrownian(phy=pruned.tree)
plot(pruned.tree)
# prep error for heterogenous variance in fixed effects 
# and sample size
vf<-varComb(varFixed(~size.n),varFixed(~1/sqrt(size.se)),varFixed(~1/sqrt(shape.se)),
            varFixed(~1/sqrt(alt.se)),varFixed(~1/sqrt(lat.se)))

# stepwise model selection AIC
fit1 <- gls(size.mean ~ shape.mean *sex.ratio + alt.mean*lat.mean, mean.df, weights = vf, correlation = bm.tree,method = "ML"); summary(fit1)
fit2 <- gls(size.mean ~ 1, mean.df, weights = vf, method = "ML", correlation = bm.tree); summary(fit2)
step.size <-stepAICc(fit1, direction = "backward",scope=list(upper=fit1,lower=fit2)); summary(step.size)
step.size$anova

MuMIn::AICc(step.size)
caper::anova(step.size)

# check fit
plot(size.gls, resid(., type="n")~fitted(.), main="Normalized Residuals v Fitted Values",
     abline=c(0,0))
res <- resid(size.gls , type="n"); dev.off(); qqnorm(res); qqline(res)
res[which.max(res)] # clysonymus kind of an outlier


###### FIGS ######

# order genus.sp phylogenetically 
mean.df <- mean.df[(rev(pruned.tree$tip.label)),]; mean.df$genus.sp = factor(mean.df$genus.sp); tip.col <- rep(1, length(mean.df$genus.sp)); tip.col <- c(rep("black",2), "#56B4E9", rep("black", 7), "#E69F00", rep("black",2))
mean.df$row_num <- 1:nrow(mean.df); mean.df$legend_entry <- paste("  ", mean.df$row_num, mean.df$genus.sp)

# Fig. 4B - size alt raw fig
fig.size.raw <-   ggplot(data = mean.df, aes(x = alt.mean, y = size.mean, fill=clade)) + 
  geom_smooth(data=subset(mean.df, clade=="era.clade"),  aes(x = alt.mean, y = size.mean), method = "gam",se=TRUE, alpha = 0.3, size=.4,formula =  y ~ poly(x, 1), inherit.aes = FALSE, colour="#56B4E9", fill="#56B4E9")+ 
  stat_cor(inherit.aes = FALSE,data=subset(mean.df, clade=="era.clade"),  aes(x = alt.mean, y = size.mean), colour="#56B4E9", size=4)+
  scale_fill_manual(name="Clades" ,values=c("era.clade"="#56B4E9", "mel.clade"="#E69F00"), labels=c("erato", "melpomene"))+
  #geom_rect(mapping=aes(xmin=with(mean.df[1,], (alt.mean-alt.sd)-20), xmax=with(mean.df[1,], (alt.mean+alt.sd)+20), ymin=with(mean.df[2,], (size.mean-size.sd)-.005), ymax=with(mean.df[2,], (size.mean+size.sd)+.01)), fill="#FF000020", alpha=0.01, inherit.aes = FALSE) +
  geom_errorbar(aes(ymin = size.mean-size.sd,ymax = size.mean+size.sd, colour = clade))+ 
  geom_errorbarh(aes(xmin = alt.mean-alt.sd, xmax = alt.mean+alt.sd, colour = clade)) +
  geom_text(label=mean.df$sp.short, size=5,vjust =- 0.3, hjust=-.2)+
  geom_point(size = 4,colour="black",shape=21) +
  xlab("Altitude of species (m)") + ylab(bquote(bold('Wing area ('*mm^2*')')))+
  scale_color_manual(name="Clades" ,values=c("era.clade"="#56B4E9", "mel.clade"="#E69F00"), labels=c("erato", "melpomene"))+
  scale_fill_manual(name="Clades" ,values=c("era.clade"="#56B4E9", "mel.clade"="#E69F00"), labels=c("erato", "melpomene"))+
  theme_bw()+
  theme(panel.border = element_rect(colour = "black", fill="transparent", size=.8),
        panel.background = element_blank(), legend.text = element_text(size=10), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(), axis.text=element_text(size=8),
        panel.grid.minor = element_blank(), axis.title=element_text(size=18,face="bold"),
        legend.text.align = 0, legend.position = c("none")); fig.size.raw

# Fig. S11 - size lat raw fig
lat<- summarise(group_by(master, genus.sp),
          lat.raw.mean=mean(latitude))
mean.df$lat.raw.mean <- lat$lat.raw.mean[match(mean.df$genus.sp, lat$genus.sp)]

fig.size.lat.raw <-   ggplot(data = mean.df, aes(x = lat.mean, y = size.mean, fill=clade)) + 
  #stat_smooth(data=mean.df,  aes(x = lat.mean, y = size.mean), method = "gam",se=TRUE, alpha = 0.3, size=.4,formula =  y ~ poly(x, 1), inherit.aes = FALSE, colour="black", fill="grey")+ 
  scale_fill_manual(name="Clades" ,values=c("era.clade"="#56B4E9", "mel.clade"="#E69F00"), labels=c("erato", "melpomene"))+
  geom_errorbar(aes(ymin = size.mean-size.sd,ymax = size.mean+size.sd, colour = clade))+ 
  geom_errorbarh(aes(xmin = lat.mean-lat.sd, xmax = lat.mean+lat.sd, colour = clade)) +
  geom_text(label=mean.df$sp.short, size=5,vjust =- 0.3, hjust=-.2)+
  geom_point(size = 4,colour="black",shape=21) +
  stat_cor(inherit.aes = FALSE,data=mean.df,  aes(x = lat.mean, y = size.mean))+
  xlab("Distance from Equator") + ylab(bquote(bold('Wing area ('*mm^2*')')))+
  scale_color_manual(name="Clades" ,values=c("era.clade"="#56B4E9", "mel.clade"="#E69F00"), labels=c("erato", "melpomene"))+
  scale_fill_manual(name="Clades" ,values=c("era.clade"="#56B4E9", "mel.clade"="#E69F00"), labels=c("erato", "melpomene"))+
  theme_bw()+
  theme(panel.border = element_rect(colour = "black", fill="transparent", size=.8),
        panel.background = element_blank(), legend.text = element_text(size=10), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(), axis.text=element_text(size=8),
        panel.grid.minor = element_blank(), axis.title=element_text(size=18,face="bold"),
        legend.text.align = 0, legend.position = c("none")); fig.size.lat.raw

######################### 2. WING SHAPE PGLS #############
###### PGLS #####
# make trees
bm.tree<-corBrownian(phy=pruned.tree)

# adds row names to the frame
row.names(mean.df)<-mean.df$genus.sp

# prep error for heterogenous variance in fixed effects 
# and sample size
vf<-varComb(varFixed(~shape.n),
            varFixed(~1/sqrt(shape.se)), varFixed(~1/sqrt(alt.se)), 
            varFixed(~1/sqrt(size.se)), varFixed(~1/sqrt(lat.se)))

# stepwise model selection AIC
fit1 <- gls(shape.mean ~  size.mean*sex.ratio +alt.mean*lat.mean , mean.df, weights = vf, method = "ML", correlation = bm.tree); summary(fit1)
fit2 <- gls(shape.mean ~ 1, mean.df, weights = vf, method = "ML", correlation = bm.tree); summary(fit2)
step.shape <-stepAICc(fit1, direction="backward",scope=list(upper=fit1,lower=fit2)); summary(step.shape ); step.shape$anova 
coef(summary(step.shape))

# check fit
hist(residuals(shape.gls))
plot(shape.gls, resid(., type="n")~fitted(.), main="Normalized Residuals v Fitted Values",
     abline=c(0,0))
res <- resid(shape.gls , type="n"); dev.off(); qqnorm(res); qqline(res)
res[which.max(res)] 

###### FIGS #### 
# order genus.sp phylogenetically 
mean.df <- mean.df[(rev(pruned.tree$tip.label)),]; mean.df$genus.sp = factor(mean.df$genus.sp)
mean.df$row_num <- 1:nrow(mean.df); mean.df$legend_entry <- paste("  ", mean.df$row_num, mean.df$genus.sp)

# Fig. 4A shape raw vs alt
fig.shape.raw <-  ggplot(data = mean.df, aes(x = alt.mean, y = shape.mean, fill=clade)) + 
  stat_smooth(data=mean.df[3:13,],  aes(x = alt.mean, y = shape.mean), method = "gam",se=TRUE, alpha = 0.3, size=.4,formula =  y ~ poly(x, 1), inherit.aes = FALSE, colour="black", fill="grey")+ 
  scale_fill_manual(name="Clades" ,values=c("era.clade"="#56B4E9", "mel.clade"="#E69F00"), labels=c("erato", "melpomene"))+
  geom_rect(mapping=aes(xmin=with(mean.df[1,], (alt.mean-alt.sd)-20), xmax=with(mean.df[1,], (alt.mean+alt.sd)+20), ymin=with(mean.df[2,], (shape.mean-shape.sd)-.005), ymax=with(mean.df[2,], (shape.mean+shape.sd)+.01)), fill="#FF000020", alpha=0.01, inherit.aes = FALSE) +
  geom_errorbar(aes(ymin = shape.mean-shape.sd,ymax = shape.mean+shape.sd, colour = clade))+ 
  geom_errorbarh(aes(xmin = alt.mean-alt.sd, xmax = alt.mean+alt.sd, colour = clade)) +
  geom_text(label=mean.df$sp.short, size=5,vjust =- 0.3, hjust=-.2)+
  stat_cor(inherit.aes = FALSE,data=subset(mean.df, genus.sp!="Heliconius_telesiphe"&genus.sp!="Heliconius_clysonymus"),aes(x = alt.mean, y = shape.mean), size=4)+
  geom_point(size = 4,colour="black",shape=21) +
  xlab("Altitude of species (m)") + ylab("Wing aspect ratio")+
  scale_color_manual(name="Clades" ,values=c("era.clade"="#56B4E9", "mel.clade"="#E69F00"), labels=c("erato", "melpomene"))+
  scale_fill_manual(name="Clades" ,values=c("era.clade"="#56B4E9", "mel.clade"="#E69F00"), labels=c("erato", "melpomene"))+
  theme_bw()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),
        panel.background = element_blank(), legend.text = element_text(size=18), legend.title =  element_text(size=18),
        panel.grid.major = element_blank(), axis.text=element_text(size=12),
        panel.grid.minor = element_blank(), axis.title=element_text(size=18,face="bold"),
        legend.text.align = 0, legend.position = "none"); fig.shape.raw 
