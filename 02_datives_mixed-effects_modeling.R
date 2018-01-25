#=======================================================================
#   Script: 02a_datives_mixed-effects_modeling.R
#   Date last changed: 9 January 2018
#
#   Rversion: 3.3.2
#   Author: MR
#   This script provides the Rcode for the mixed-effects model in the
#   dissertation that was run on the full dataset.
#=======================================================================

#----SETUP: LIBRARIES-----------------------------------------
library(reshape2);library(plyr);library(dplyr);library(magrittr);library(lme4);library(ggplot2);library(gridExtra); library(JGmisc); library(beepr); library(effects); library(Hmisc); library(car)

setwd("/Users/melanierothlisberger/Dropbox/Universitaet/04_Statistics/statistical analyses/20170209_THESIS_full_model")
source("/Users/melanierothlisberger/Dropbox/Universitaet/03_Programming/R/helper.R")

options(scipen = 999)

#-----SETUP: DATA------------------------
data <- read.delim("OA/datives.txt")


#---SETUP: PLOTTING----------------
theme_mr = theme_set(theme_light())
theme_mr = theme_update(axis.text = element_text(size = rel(1.0), color="black"),
                        axis.ticks = element_line(colour = "grey90", size = 0.25),
                        axis.title = element_text(size=rel(0.9)),
                        panel.border = element_rect(color = "black"),
                        strip.background=element_rect(fill="grey95", color="black"), 
                        strip.text.x=element_text(color="black"))


# -- Set the reference levels and the contrasts -------------------------------
# predictions are for the prepositional dative, i.e. all levels are set for the prototypical ditransitive dative
data$Resp <- relevel(data$Resp, ref = "do")
data$Nativity <- relevel(data$Nativity, ref = "L1")
data$Mode <- relevel(data$Mode, ref = "spoken")
data$Corpus <- relevel(data$Corpus, ref = "ice")
data$VerbSemantics <- relevel(data$VerbSemantics, ref = "a")
data$RecGivenness <- relevel(data$RecGivenness, ref = "given")
data$ThemeGivenness <- relevel(data$ThemeGivenness, ref = "new")
data$RecDefiniteness <- relevel(data$RecDefiniteness, ref = "def")
data$ThemeDefiniteness <- relevel(data$ThemeDefiniteness, ref = "indef")
data$RecBinComplexity <- relevel(data$RecBinComplexity, ref = "simple")
data$ThemeBinComplexity <- relevel(data$ThemeBinComplexity, ref = "complex")
data$RecBinAnimacy <- relevel(data$RecBinAnimacy, ref = "animate")
data$ThemeBinAnimacy <- relevel(data$ThemeBinAnimacy, ref = "inanimate")
data$RecPerson <- relevel(data$RecPerson, ref = "local")
data$ThemeConcreteness <- relevel(data$ThemeConcreteness, ref="non-concrete")
data$Persistence <- relevel(data$Persistence, ref = "none")
data$RecPron <- relevel(data$RecPron, ref="pron")
data$ThemePron <- relevel(data$ThemePron, ref="non-pron")
data$PrimeType <- relevel(data$PrimeType, ref="none")
data$PrimeTypePruned <- relevel(data$PrimeTypePruned, ref="none")

# relevel Variety
data$Variety <- factor(data$Variety, levels=c("CAN", "GB", "HK", "IND", "IRE", "JA", "NZ", "PHI", "SIN"))
data$Variety.Sum <- data$Variety
contrasts(data$Variety.Sum) = contr.sum(9)



# -- Filter infrequent ranef levels -------------------------------

head(sort(table(data$RecHeadPlain), decreasing=T))
hist(as.vector(table(data$RecHeadPlain)),
     breaks=length(levels(data$RecHeadPlain)), 
     col="gray",
     xlim=c(0,25))

quantile(as.vector(table(data$RecHeadPlain)),
         probs = seq(0,1,0.05))
quantile(as.vector(table(data$ThemeHeadPlain)),
         probs = seq(0,1,0.05))
# go for the 90% threshold
data$RecHeadFilter <- filter.infrequent(data$RecHeadLemma, threshold = 4)
data$ThemeHeadFilter <- filter.infrequent(data$ThemeHeadLemma, threshold = 8)


# -- Random effect structure ------------------------------------

# RANEFS: Subcorpus, FileID, TextID, SpeakerID, GenreFine, GenreCoarse, Verb, VerbForm, RecHeadPlain, ThemeHeadPlain, VerbThemeLemma, VerbSense 

# (Variety | GenreFine)
# (1 | FileID/SpeakerID)
# (RecBinAnimacy|Verb)
# (Verb|VerbSense)
# ThemeHeadLemma
# RecHeadLemma
# (Verb | ThemeHeadFilter)

# -- Fixed effect structure -------------------------------------
# RecBinAnimacy
# ThemeBinAnimacy
# RecBinComplexity
# ThemeBinComplexity
# RecGivenness
# ThemeGivenness
# RecDefiniteness
# ThemeDefiniteness
# RecThematicity
# ThemeThematicity
# RecHeadFreq
# ThemeHeadFreq
# TypeTokenRatio
# RecPerson: exclude as coding is not adequate
# ThemeConcreteness: exclude due to overlap with verb semantics
# PrimeTypePruned
# RecPron
# ThemePron
# logWeightRatio
# ThemePron:ThemeDefiniteness?
# VerbSemantics
# Corpus
# Use all language external factors as interaction terms to see how internal factors vary across those

# -- Model selection ----------------------------------------------

# -- NOT CONV:BEYOND OPTIMAL MODEL------------------------------------
t1 <- proc.time()
r00 <- glmer(Resp~
               # Random factors
               # ( 0 + Variety | Verb)
               (1 | Verb/VerbSense)
             + (1 | GenreFine/FileID/SpeakerID)
             + (1 | GenreCoarse)
             + (1 | ThemeHeadFilter) 
             + (1 | RecHeadFilter)
             
             # Fixed factors with all interactions
             + Variety.Sum * (
               Corpus   
               + RecBinAnimacy
               #+ ThemeBinAnimacy
               + RecBinComplexity
               + ThemeBinComplexity
               + RecGivenness
               + ThemeGivenness
               + RecDefiniteness
               + ThemeDefiniteness
               + z.RecThematicity
               + z.ThemeThematicity
               + z.RecHeadFreq
               + z.ThemeHeadFreq
               + z.TypeTokenRatio
               + PrimeType
               + RecPron
               + ThemePron
               + z.logWeightRatio
               + VerbSemantics
               + Mode)
             + Mode:RecBinAnimacy
             + Mode:RecBinComplexity
             + Mode:ThemeBinComplexity
             + Mode:RecGivenness
             + Mode:ThemeGivenness
             + Mode:RecDefiniteness
             + Mode:ThemeDefiniteness
             + Mode:z.RecThematicity
             + Mode:z.ThemeThematicity
             + Mode:z.RecHeadFreq
             + Mode:z.ThemeHeadFreq
             + Mode:z.TypeTokenRatio
             + Mode:PrimeType
             + Mode:RecPron
             + Mode:ThemePron
             + Mode:z.logWeightRatio
             + Mode:VerbSemantics, 
             data=data, family = binomial,
             glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 2e5), 
                          calc.derivs = FALSE))

beep(3)
t2 <- proc.time() - t1
# -- NOT CONV: intercept-only model ----------------------------------------------------
t1 <- proc.time()
r01 <- glmer(Resp~
               # Random factors
               # ( 0 + Variety | Verb)
               (1 | Verb/VerbSense)
             + (1 | GenreFine/FileID/SpeakerID)
             + (1 | GenreCoarse)
             + (1 | ThemeHeadFilter) 
             + (1 | RecHeadFilter)
             + 1,
             data=data, family = binomial,
             glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6)))

t2 <- proc.time() - t1


# -- Model with same call structure as ICE (see Röthlisberger et al. 2017) -------------------
# converges successively --> add random factors next
t1 <- proc.time()
r00 <- glmer(Resp~
               # Random factors
               # ( 0 + Variety | Verb)
               (1 | Verb/VerbSense)
             + (1 | GenreCoarse/GenreFine/FileID)
             + (1 | ThemeHeadFilter) 
             #+ (1 | RecHeadFilter)
             
             # Fixed factors with all interactions
            # + Corpus   
             + RecBinAnimacy
             #+ ThemeBinAnimacy
             + RecBinComplexity
             + ThemeBinComplexity
             + RecGivenness
             + ThemeGivenness
             + RecDefiniteness
             + ThemeDefiniteness
            # + z.RecThematicity
            # + z.ThemeThematicity
            # + z.RecHeadFreq
            # + z.ThemeHeadFreq
            # + z.TypeTokenRatio
            # + PrimeType
             + ThemePron
            # + VerbSemantics
             + Mode
             + Variety.Sum * (
                 z.logWeightRatio
               + RecPron
               ),
             data=data, family = binomial,
             glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6)))

beep(3)
t2 <- proc.time() - t1

# add ranef RecHead
t1 <- proc.time()
r01 <- glmer(Resp~
               # Random factors
               # ( 0 + Variety | Verb)
               (1 | Verb/VerbSense)
             + (1 | GenreCoarse/GenreFine/FileID)
             + (1 | ThemeHeadFilter) 
             + (1 | RecHeadFilter)
             
             # Fixed factors with all interactions
             # + Corpus   
             + RecBinAnimacy
             #+ ThemeBinAnimacy
             + RecBinComplexity
             + ThemeBinComplexity
             + RecGivenness
             + ThemeGivenness
             + RecDefiniteness
             + ThemeDefiniteness
             # + z.RecThematicity
             # + z.ThemeThematicity
             # + z.RecHeadFreq
             # + z.ThemeHeadFreq
             # + z.TypeTokenRatio
             # + PrimeType
             + ThemePron
             # + VerbSemantics
             + Mode
             + Variety.Sum * (
               z.logWeightRatio
               + RecPron
             ),
             data=data, family = binomial,
             glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6)))

beep(3)
t2 <- proc.time() - t1

anova(r00, r01) # RecHead adds significantly, keep in model (even with Bonferroni-Correction, 0.05/2=0.025)

# add ranef SpeakerID --> leave out SpeakerID as it takes too long to converge
t1 <- proc.time()
r02 <- glmer(Resp~
               # Random factors
               # ( 0 + Variety | Verb)
               (1 | Verb/VerbSense)
             + (1 | GenreCoarse/GenreFine/FileID)
             + (1 | SpeakerID)
             + (1 | ThemeHeadFilter) 
             + (1 | RecHeadFilter)
             
             # Fixed factors with all interactions
             # + Corpus   
             + RecBinAnimacy
             #+ ThemeBinAnimacy
             + RecBinComplexity
             + ThemeBinComplexity
             + RecGivenness
             + ThemeGivenness
             + RecDefiniteness
             + ThemeDefiniteness
             # + z.RecThematicity
             # + z.ThemeThematicity
             # + z.RecHeadFreq
             # + z.ThemeHeadFreq
             # + z.TypeTokenRatio
             # + PrimeType
             + ThemePron
             # + VerbSemantics
             + Mode
             + Variety.Sum * (
               z.logWeightRatio
               + RecPron
             ),
             data=data, family = binomial,
             glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 2e5)))

beep(3)
t2 <- proc.time() - t1

# -- f01: full model with reduced random structure ----------------------
# add all fixef and all possible interactions
t1 <- proc.time()
f00 <- glmer(Resp~
               # Random factors
               (1 | Verb/VerbSense)
             + (1 | GenreCoarse/GenreFine/FileID)
             + (1 | ThemeHeadFilter) 
             + (1 | RecHeadFilter)
             
             # Fixed factors 
             + ThemeBinAnimacy
             + RecBinComplexity
             + RecGivenness
             + ThemeGivenness
              + z.RecThematicity
              + z.ThemeThematicity
              + z.RecHeadFreq
              + z.ThemeHeadFreq
              + z.TypeTokenRatio
             + ThemePron
             + VerbSemantics
             + Mode
             
             # Interaction effects with mode
             + Mode:PrimeTypePruned
             + Mode:Variety.Sum
             + Mode:RecGivenness
             + Mode:ThemeDefiniteness
             
             # Interaction effects with Variety
             + Variety.Sum * (
               z.logWeightRatio
               + RecPron
               + RecDefiniteness
               + ThemeDefiniteness
               + Corpus
               + PrimeTypePruned
               + ThemeBinComplexity
               + RecBinAnimacy
             ),
             data=data, family = binomial,
             glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 2e5)))

beep(3)
t2 <- proc.time() - t1
summary(f00)
Anova(f00)
summary(f00)$varcor # GenreCoarse could be left out but I retain here since it is important as interaction with other Corpus-related ranefs

round(summary(f00)$coefficients,3)
# discard interactions whose estimates (absolute values) are less than their standard errors, i.e.: Mode and Variety, Mode and ThemeDefiniteness, Variety and ThemeDefiniteness, Variety and ThemeComplexity

# verify the inclusion of interactions
f01 <- glmer(Resp~
               # Random factors
               (1 | Verb/VerbSense)
             + (1 | GenreCoarse/GenreFine/FileID)
             + (1 | ThemeHeadFilter) 
             + (1 | RecHeadFilter)
             
             # Fixed factors 
             + ThemeBinAnimacy
             + RecBinComplexity
             + RecGivenness
             + ThemeGivenness
             + z.RecThematicity
             + z.ThemeThematicity
             + z.RecHeadFreq
             + z.ThemeHeadFreq
             + z.TypeTokenRatio
             + ThemePron
             + VerbSemantics
             + Mode
             + ThemeDefiniteness
             + ThemeBinComplexity
             
             # Interaction effects with mode
             + Mode:PrimeTypePruned
             #+ Mode:Variety.Sum
             + Mode:RecGivenness
             #+ Mode:ThemeDefiniteness
             
             # Interaction effects with Variety
             + Variety.Sum * (
               z.logWeightRatio
               + RecPron
               + RecDefiniteness
               + Corpus
               + PrimeTypePruned
               + RecBinAnimacy
             ),
             data=data, family = binomial,
             glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6)))

beep(3)
anova(f00, f01) # no sign. difference, f01 has lower AIC
Anova(f01)
summary(f01)

somers.mer(f01)
ggLogit.plot(f01, data) # looks like the model is under-predicting the lower proportion and overpredicting the upper proportion. That's not good.

# -- .. evaluation of the f01 model ----------------

collin.fnc.mer(f01)$cnumber  # WTF! 89.5
kappa.mer(f01) # 40.1
as.data.frame(sort(vif.mer(f01))) # use deviation coding for Semantics? And what is going on with ThemeConcreteness??
overdisp.mer(f01) # good

dat.fits <- fitted(f01)
sum(ifelse(dat.fits > .5, "pd", "do") == data$Resp)/nrow(data) # 93.3 %


# bin and plot residuals
# add binary 1-0 response
data$NumResp <- as.numeric(data$Resp)-1
dat.resid.bin <- as.data.frame(binned.resids(dat.fits, data$NumResp - dat.fits, nclass = 30)$binned)
ggplot(dat.resid.bin, aes(xbar, ybar)) + 
  geom_ribbon(aes(ymin = ybar - 2*se, ymax = ybar + 2*se), alpha = 0.1) + 
  geom_point() + labs(x = "estimated Pr(pd-dative)", y = "mean residual") +
  geom_smooth(method="loess",se=F) # plot a loess smoother
# that looks quite good!

# outliers:
plot(f01, id = .01, idLabels = ~.obs)
# might not be too bad, but there are some over resid=5

# plot hcluster to see which factors go together
#subset for to visualize clusters of variance: Hierarchical cluster analysis on variables using Spearman correlation using varclus() in the Hmisc package
attach(data)
plot(varclus(~ RecPron  + RecGivenness + ThemeGivenness + RecDefiniteness + ThemeDefiniteness + RecBinComplexity + VerbSemantics + ThemeBinComplexity  + RecBinAnimacy + ThemePron + z.logWeightRatio + Variety.Sum + Mode + PrimeTypePruned + z.RecHeadFreq + z.ThemeHeadFreq + z.RecThematicity + z.ThemeThematicity + z.TypeTokenRatio + Corpus +  Variety.Sum:RecPron + Variety.Sum:z.logWeightRatio + Variety.Sum:Corpus + Variety.Sum:RecBinAnimacy + Variety.Sum:PrimeTypePruned + Variety.Sum:RecDefiniteness - 1), abbrev=TRUE) # it uses spearman

coll <- colldiag.mer(f01, scale = TRUE, center = FALSE, add.intercept = TRUE)


# NEXT: take out more interactions: Mode:PrimeType, Mode:RecGivenness
# THEN NEXT fixed effects: z.RecThematicity; z.ThemeHeadFreq;

# -- f02 ------------------------
f02 <- glmer(Resp~
               # Random factors
               (1 | Verb/VerbSense)
             + (1 | GenreCoarse/GenreFine/FileID)
             + (1 | ThemeHeadFilter) 
             + (1 | RecHeadFilter)
             
             # Fixed factors 
             + ThemeBinAnimacy
             + RecBinComplexity
             + RecGivenness
             + ThemeGivenness
             + z.RecThematicity
             + z.ThemeThematicity
             + z.RecHeadFreq
             + z.ThemeHeadFreq
             + z.TypeTokenRatio
             + ThemePron
             + VerbSemantics
             + Mode
             + ThemeDefiniteness
             + ThemeBinComplexity
             
             # Interaction effects with mode
             #+ Mode:PrimeTypePruned
             #+ Mode:Variety.Sum
             #+ Mode:RecGivenness
             #+ Mode:ThemeDefiniteness
             
             # Interaction effects with Variety
             + Variety.Sum * (
               z.logWeightRatio
               + RecPron
               + RecDefiniteness
               + Corpus
               + PrimeTypePruned
               + RecBinAnimacy
             ),
             data=data, family = binomial,
             glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6)))

beep(3)
anova(f01, f02) # n.s. different, justified, f02 has lower AIC
Anova(f02)
summary(f02)

# -- f03 ------------------------
# NEXT: consider leaving interactions out with PrimeType, RecDef and RecAni
# Then NEXT: VerbSemantics strongly correlates with ThemeConcreteness!! consider leaving ThemeConcreteness out, other fixed effects: z.RecThematicity and z.ThemeHeadFreq

f03 <- glmer(Resp~
               # Random factors
               (1 | Verb/VerbSense)
             + (1 | GenreCoarse/GenreFine/FileID)
             + (1 | ThemeHeadFilter) 
             + (1 | RecHeadFilter)
             
             # Fixed factors 
             + ThemeBinAnimacy
             + RecBinComplexity
             + RecGivenness
             + ThemeGivenness
             + z.RecThematicity
             + z.ThemeThematicity
             + z.RecHeadFreq
             + z.ThemeHeadFreq
             + z.TypeTokenRatio
             + ThemePron
             + VerbSemantics
             + Mode
             + ThemeDefiniteness
             + ThemeBinComplexity
             + RecDefiniteness
             + PrimeTypePruned
             + RecBinAnimacy
             
             # Interaction effects with mode
             #+ Mode:PrimeTypePruned
             #+ Mode:Variety.Sum
             #+ Mode:RecGivenness
             #+ Mode:ThemeDefiniteness
             
             # Interaction effects with Variety
             + Variety.Sum * (
               z.logWeightRatio
               + RecPron
               + Corpus),
             data=data, family = binomial,
             glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6)))

anova(f02,f03) # n.s. different, interactions can be left out
anova(f01,f03) # n.s.

# -- .. evaluation of the f03 model ----------------

collin.fnc.mer(f03)$cnumber  # WTF! 74.74!
kappa.mer(f03) # 
as.data.frame(sort(vif.mer(f03))) # use deviation coding for Semantics? And what is going on with ThemeConcreteness??
overdisp.mer(f03) # good

dat.fits <- fitted(f03)
sum(ifelse(dat.fits > .5, "pd", "do") == data$Resp)/nrow(data) # 93.2 %


# bin and plot residuals
# add binary 1-0 response
data$NumResp <- as.numeric(data$Resp)-1
dat.resid.bin <- as.data.frame(binned.resids(dat.fits, data$NumResp - dat.fits, nclass = 30)$binned)
ggplot(dat.resid.bin, aes(xbar, ybar)) + 
  geom_ribbon(aes(ymin = ybar - 2*se, ymax = ybar + 2*se), alpha = 0.1) + 
  geom_point() + labs(x = "estimated Pr(pd-dative)", y = "mean residual") +
  geom_smooth(method="loess",se=F) # plot a loess smoother
# that looks quite good!

# outliers:
plot(f03, id = .01, idLabels = ~.obs)
# might not be too bad, but there are some over resid=5
ggLogit.plot(f03, data) # good, no S-curve
somers.mer(f03)
MuMIn::r.squaredGLMM(f03) # 

# -- f04 ------------------------
# f04: remove fixed effects: z.ThemeThem, RecHeadFreq, ThemeHeadFreq
summary(f03)
Anova(f03)
f04 <- glmer(Resp~
               # Random factors
               (1 | Verb/VerbSense)
             + (1 | GenreCoarse/GenreFine/FileID)
             + (1 | ThemeHeadFilter) 
             + (1 | RecHeadFilter)
             
             # Fixed factors 
             + ThemeBinAnimacy
             + RecBinComplexity
             + RecGivenness
             + ThemeGivenness
             + z.RecThematicity
             #+ z.ThemeThematicity
             #+ z.RecHeadFreq
             #+ z.ThemeHeadFreq
             + z.TypeTokenRatio
             + ThemePron
             + VerbSemantics
             + Mode
             + ThemeDefiniteness
             + ThemeBinComplexity
             + RecDefiniteness
             + PrimeTypePruned
             + RecBinAnimacy
             
             
             # Interaction effects with Variety
             + Variety.Sum * (
               z.logWeightRatio
               + RecPron
               + Corpus),
             data=data, family = binomial,
             glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6)))
beep(3)
anova(f03, f04) # ns. different
summary(f04)
Anova(f04)

# -- f05 ------------------------
# consider leaving out Mode completely
f05 <- glmer(Resp~
               # Random factors
               (1 | Verb/VerbSense)
             + (1 | GenreCoarse/GenreFine/FileID)
             + (1 | ThemeHeadFilter) 
             + (1 | RecHeadFilter)
             
             # Fixed factors 
             + ThemeBinAnimacy
             + RecBinComplexity
             + RecGivenness
             + ThemeGivenness
             + z.RecThematicity
             #+ z.ThemeThematicity
             #+ z.RecHeadFreq
             #+ z.ThemeHeadFreq
             + z.TypeTokenRatio
             + ThemePron
             + VerbSemantics
             #+ Mode
             + ThemeDefiniteness
             + ThemeBinComplexity
             + RecDefiniteness
             + PrimeTypePruned
             + RecBinAnimacy
             
             
             # Interaction effects with Variety
             + Variety.Sum * (
               z.logWeightRatio
               + RecPron
               + Corpus),
             data=data, family = binomial,
             glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6)))
beep(3)
anova(f04, f05) # ns. different
summary(f05)
Anova(f05)

# -- f06 ------------------------
# exclude RecThematicity and ThemePron
f06 <- glmer(Resp~
               # Random factors
               (1 | Verb/VerbSense)
             + (1 | GenreCoarse/GenreFine/FileID)
             + (1 | ThemeHeadFilter) 
             + (1 | RecHeadFilter)
             
             # Fixed factors 
             + ThemeBinAnimacy
             + RecBinComplexity
             + RecGivenness
             + ThemeGivenness
            # + z.RecThematicity
             #+ z.ThemeThematicity
             #+ z.RecHeadFreq
             #+ z.ThemeHeadFreq
             + z.TypeTokenRatio
             #+ ThemePron
             + VerbSemantics
             #+ Mode
             + ThemeDefiniteness
             + ThemeBinComplexity
             + RecDefiniteness
             + PrimeTypePruned
             + RecBinAnimacy
             
             
             # Interaction effects with Variety
             + Variety.Sum * (
               z.logWeightRatio
               + RecPron
               + Corpus),
             data=data, family = binomial,
             glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6)))
beep(3)
anova(f06, f05) # ns. different
summary(f06)
Anova(f06)

# -- f07 ------------------------
# leave out VerbSemantics
f07 <- glmer(Resp~
               # Random factors
               (1 | Verb/VerbSense)
             + (1 | GenreCoarse/GenreFine/FileID)
             + (1 | ThemeHeadFilter) 
             + (1 | RecHeadFilter)
             
             # Fixed factors 
             + ThemeBinAnimacy
             + RecBinComplexity
             + RecGivenness
             + ThemeGivenness
             # + z.RecThematicity
             #+ z.ThemeThematicity
             #+ z.RecHeadFreq
             #+ z.ThemeHeadFreq
             + z.TypeTokenRatio
             #+ ThemeConcreteness
             #+ ThemePron
             #+ VerbSemantics
             #+ Mode
             + ThemeDefiniteness
             + ThemeBinComplexity
             + RecDefiniteness
             + PrimeTypePruned
             + RecBinAnimacy
             
             
             # Interaction effects with Variety
             + Variety.Sum * (
               z.logWeightRatio
               + RecPron
               + Corpus),
             data=data, family = binomial,
             glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6)))
beep(3)
anova(f06, f07) # ns. different
summary(f07)
Anova(f07)


# -- .. evaluation of the f07 model ----------------

collin.fnc.mer(f07)$cnumber  # 13.5
kappa.mer(f07) # 13.8
as.data.frame(sort(vif.mer(f07))) 
overdisp.mer(f07) # good

dat.fits <- fitted(f07)
(d.acc <- sum(ifelse(dat.fits > .5, "pd", "do") == data$Resp)/nrow(data)) # 93.2 %


# bin and plot residuals
# add binary 1-0 response
data$NumResp <- as.numeric(data$Resp)-1
dat.resid.bin <- as.data.frame(binned.resids(dat.fits, data$NumResp - dat.fits, nclass = 30)$binned)
ggplot(dat.resid.bin, aes(xbar, ybar)) + 
  geom_ribbon(aes(ymin = ybar - 2*se, ymax = ybar + 2*se), alpha = 0.1) + 
  geom_point() + labs(x = "estimated Pr(pd-dative)", y = "mean residual") +
  geom_smooth(method="loess",se=F) # plot a loess smoother
# that looks quite good!

# outliers:
plot(f07, id = .01, idLabels = ~.obs)
# might not be too bad, but there are some over resid=5
ggLogit.plot(f07, data) # good, slight S-curve
somers.mer(f07)
(r.sq.f07 <- MuMIn::r.squaredGLMM(f07)) # R2m = 0.3833460, R2c = 0.8364792

# one-tailed binomial test to compare accuracy to baseline
# Assuming you have a vector of predictions, "preds", and a baseline probability, "base.p", You can do it like so - testing the null hypothesis that the baseline is greater than the model's accuracy:
# Baseline accuracy is 73.84% (for do)
preds <- ifelse(dat.fits > .5, "pd", "do")
hits <- sum(preds == data$Resp)
binom.test(x = hits, n = nrow(data), p = 0.7384, alternative = "greater")


data$fitted.probs <- fitted(f07, type='response')
data$pred.resp <- ifelse(data$fitted.probs > .5, 'pd', 'do')
tab.fit <- table(obs=data$Resp, pred=data$pred.resp)
tab.fit
prop.table(tab.fit, margin=2)

sum(diag(tab.fit))/nrow(data)


# QQ plots of ranefs

rf.th <- ranef(f07)[[2]]
names(rf.th) <- "Intercept"
ggQQ.plot(rf.th, Intercept)

th <- rf.th[order(rf.th$Intercept), , drop=F]
# tables of most extreme adjustments
head(th, 10)
tail(th, 10)

rf.verb <- ranef(f07)[[5]]
names(rf.verb) <- "Intercept"
ggQQ.plot(rf.verb, Intercept)


# -- .. check for SinE ----------------
data$Variety.Sum2 <- data$Variety
data$Variety.Sum2 <- factor(data$Variety.Sum2, levels=c("CAN","GB","HK","IND","IRE", "JA","NZ","SIN", "PHI"))
contrasts(data$Variety.Sum2) <- contr.sum(9)

f07.SIN <- glmer(Resp~
               # Random factors
               (1 | Verb/VerbSense)
             + (1 | GenreCoarse/GenreFine/FileID)
             + (1 | ThemeHeadFilter) 
             + (1 | RecHeadFilter)
             
             # Fixed factors 
             + ThemeBinAnimacy
             + RecBinComplexity
             + RecGivenness
             + ThemeGivenness
             # + z.RecThematicity
             #+ z.ThemeThematicity
             #+ z.RecHeadFreq
             #+ z.ThemeHeadFreq
             + z.TypeTokenRatio
             #+ ThemeConcreteness
             #+ ThemePron
             #+ VerbSemantics
             #+ Mode
             + ThemeDefiniteness
             + ThemeBinComplexity
             + RecDefiniteness
             + PrimeTypePruned
             + RecBinAnimacy
             
             
             # Interaction effects with Variety
             + Variety.Sum2 * (
               z.logWeightRatio
               + RecPron
               + Corpus),
             data=data, family = binomial,
             glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e8)))
beep(3)

summary(f07.SIN)


# -- Validation / Bootstrapping-------------------------------

split(data, factor(runif(nrow(data))>0.25)) -> xx # split dataset approximately 75 vs. 25%
table(xx$'TRUE'$Resp)  # training set for model
table(xx$'FALSE'$Resp) # test set for model

pcorrect <- array(dim=100) # allocate space for results
for (i in 1:100) {
  cat('run', i + 1, '...') # count out the runs
  
  # step 1: split the data into a training and test set
  split(data, factor(runif(nrow(data))>0.25)) -> xx
  
  # step 2: estimate model parameters with the training set
  xx$'TRUE'$Resp <- relevel(xx$'TRUE'$Resp, ref = "do")
  xx$'FALSE'$Resp <- relevel(xx$'FALSE'$Resp, ref = "do")
  mod <- update(f07, data = xx$'TRUE', na.rm=T)
  
  # step 3: get model predictions for the test set, allow new levels of predictors
  predict(mod, xx$'FALSE', allow.new.levels=T) -> pred
  
  # bookkeeping 1: tabulate the responses
  table((pred > 0.0)==(xx$'FALSE'$Resp=="pd")) -> tab
  
  # bookkeeping 2: save the probability of a correct answer
  pcorrect[i] <- tab[2]/(tab[1]+tab[2])
}
mean(pcorrect)
hist(pcorrect)

table((predict(f07)>0.0)==(data$Resp=="pd")) # that works

# ==> this took ages and was run on the Amazon Server

# -- Plotting of data structure f07 --------------------------
# random intercept for verb
ranverbs <- as.data.frame(ranef(f07)$Verb)
ranverbs$verb <- rownames(ranverbs)
names(ranverbs) <- c("adjustment", "verb")
ranverbs <- ranverbs[with(ranverbs, order(-adjustment)),]
ranverbs$abs <- abs(ranverbs$adjustment)
ranverbs$variant <- ifelse(ranverbs$adjustment<0, "ditransitive", "prepositional")
ranverbs$jitterd <- jitter(as.numeric(factor(ranverbs$variant)))

require("ggrepel")
set.seed(42)

(plot.rf.verb <- ggplot(ranverbs, aes(abs,adjustment)) +
  geom_point(aes(shape=variant), size=5) +
  geom_text_repel(data = subset(ranverbs, abs(adjustment) > 2), aes(label=verb), size=6, box.padding = unit(1, "lines"), force=2) +
  geom_hline(yintercept=0, linetype="longdash") +
  annotate("text", x = 0.4, y = -2.5, label = "ADJUSTMENT FOR\n DITRANSITIVE", fontface="bold", size=6) +
  annotate("text", x=0.4, y=2.5, label = "ADJUSTMENT FOR\n PREPOSITIONAL", fontface="bold", size=6) +
  scale_shape_discrete(name="Variant adjusted for") +
  theme(legend.position = "none") +
  labs(x="absolute adjustment to intercept", y = "adjustment to intercept") +
  theme(axis.title = element_text(size=16)))

pdf("Rplot_glmm_adjustment_verb.pdf", width = 12, height = 10)
plot.rf.verb
dev.off()

# Random intercept for theme
ranth <- as.data.frame(ranef(f07)$ThemeHeadFilter)
ranth$theme <- rownames(ranth)
names(ranth) <- c("adjustment", "theme")
ranth <- ranth[with(ranth, order(-adjustment)),]
ranth$abs <- abs(ranth$adjustment)
ranth$variant <- ifelse(ranth$adjustment<0, "ditransitive", "prepositional")
ranth$jitterd <- jitter(as.numeric(factor(ranth$variant)))

require("ggrepel")
set.seed(42)

(plot.rf.theme <- ggplot(ranth, aes(abs,adjustment)) +
  geom_point(aes(shape=variant), size=5, alpha=0.7) +
  geom_text_repel(data = subset(ranth, abs(adjustment) > 1.71), aes(label=theme), size=6, box.padding = unit(1, "lines"), force=3) +
  geom_hline(yintercept=0, linetype="longdash") +
  annotate("text", x = 0.6, y = -2.5, label = "ADJUSTMENT FOR\n DITRANSITIVE", fontface="bold", size=6) +
  annotate("text", x=0.6, y=2.6, label = "ADJUSTMENT FOR\n PREPOSITIONAL", fontface="bold", size=6) +
  theme(legend.position = "none") +
  labs(x="absolute adjustment to intercept", y = "adjustment to intercept") +
  theme(axis.title = element_text(size=16)))

pdf("Rplot_glmm_adjustment_theme.pdf", width = 9, height = 7)
plot.rf.theme
dev.off()

# mosaicplots  
require(gridExtra)
require(vcd)
data$Variant <- data$Resp
levels(data$Variant) <- c("ditransitive", "prepositional")
data$Variant <- factor(data$Variant, levels=c("prepositional", "ditransitive"))

# tests:
mosaicplot(with(data, table(Variant, RecPron)), main = "RecPron : Variant", color = c("gray", "black"), cex.axis = 1.5, xlab="", ylab="")

mosaic(with(data, table(Variant, RecPron)), main = "", shade = TRUE, legend=FALSE, gp=gpar(fill=matrix(c("dimgray", "dimgray", "gainsboro" ,"gainsboro"),2)))

# create mosaicplots
mp.recpron <- grid.grab(mosaic(with(data, table(Variant, RecPron)), main = "", shade = TRUE, legend=FALSE, gp=gpar(fill=matrix(c("dimgray", "dimgray", "gainsboro" ,"gainsboro"),2))))
mp.corpus <- grid.grab(mosaic(with(data, table(Variant, Corpus)), main = "", shade = TRUE, legend=FALSE, gp=gpar(fill=matrix(c("dimgray", "dimgray", "gainsboro" ,"gainsboro"),2))))

pdf("Rplot_mosaicplots_recpronXcorpus", width=7, height = 5)
grid.arrange(mp.recpron, mp.corpus, ncol=2)
dev.off()

data$Variant <- factor(data$Variant, levels=c("ditransitive", "prepositional"))
cdplot(Variant ~ z.logWeightRatio, data = data, xlab="relative length\n (log-scaled, centred and standardised)")

pdf("Rplot_mosaicplots_length.pdf", width=7, height = 4)
cdplot(Variant ~ z.logWeightRatio, data = data, xlab="relative length\n (log-scaled, centred and standardised)")
dev.off()




# barplots with percentages
ggplot(data, aes(Variant, fill=RecPron)) +
  geom_bar(aes(fill=RecPron), position="fill", width=0.7, color="black") +
  labs(title = "", y="proportion of tokens", x="dative variant") + 
  theme(legend.position="bottom") + 
  theme(axis.text.x = element_text(size=16)) +
  theme(legend.text=element_text(size=20)) +
  #geom_text(data = freq_sorted, aes(x=Variety,y=Freq, label = Tokens), vjust=2, color=ifelse(freq_sorted$Resp=="pd", "black", "white"), size=8) + 
  scale_fill_manual(values = c("white", "black"))




# -- .. Effects plots ---------------------------
eff <- Effect(c("Variety.Sum", "RecPron"), f07) %>% as.data.frame() %>% droplevels()
plot.recpron <- plot(eff, type="link", zvar="Variety", multiline = TRUE, ylab = "log odds of prepositional dative", xlab = "Variety")

eff$Variety.Sum <- factor(eff$Variety.Sum, levels=c("GB", "CAN", "IRE", "NZ", "JA", "SIN", "IND", "HK", "PHI"))


(plot.recpron2 <- ggplot(eff, aes(Variety.Sum, fit, group=RecPron, shape=RecPron)) +
  geom_hline(yintercept = 0.5, linetype="dashed") +
  ylim(0,1) +
  geom_line(aes(linetype=RecPron)) +
  geom_point(size=3) +
  #scale_linetype_discrete(name="recipient pronominality", breaks = c("non-pron", "pron"), labels = c("non-pron", "pron")) +
  scale_linetype_manual(values = c("dashed", "solid"), name= 'recipient pronominality', labels = c('non-pron', 'pron'), guide=FALSE) +
  scale_shape_manual(values = c(17, 16), name="Recipient pronominality", labels = c("non-pron", "pron")) +
  labs(y="probability of prepositional dative", x="") +
  theme(legend.position = "bottom", legend.text = element_text(size=13), legend.title = element_text(size=13, face="bold")) +
  guides(shape=guide_legend(order=1, title.position = "top", title.hjust = 0.5), linetype=FALSE) +
  #theme(legend.background = element_rect(linetype="solid")) +
    theme(axis.text.x = element_text(size=16)) +
    theme(axis.title = element_text(size=14)))

pdf("Rplot_effects_f07_VarXRecPron.pdf", width = 7, height = 5)
plot.recpron2
dev.off()


eff <- Effect(c("Variety.Sum", "Corpus"), f07) %>% as.data.frame() %>% droplevels()
eff$Variety.Sum <- factor(eff$Variety.Sum, levels=c("GB", "CAN", "IRE", "NZ", "JA", "SIN", "IND", "HK", "PHI"))

(plot.corpus <- ggplot(eff, aes(Variety.Sum, fit, group=Corpus, shape=Corpus)) +
    geom_hline(yintercept = 0.5, linetype="dashed") +
    ylim(0,1) +
    geom_line(aes(linetype=Corpus)) +
    geom_point(size=3) +
    #scale_linetype_discrete(name="recipient pronominality", breaks = c("non-pron", "pron"), labels = c("non-pron", "pron")) +
    scale_linetype_discrete(name= 'Corpus', labels = c('glowbe', 'ice'), guide=FALSE) +
    scale_shape_discrete(name="Corpus", labels = c('glowbe', 'ice')) +
    labs(y="probability of prepositional dative", x="") +
    theme(legend.position = "bottom", legend.text = element_text(size=13), legend.title = element_text(size=13, face="bold")) +
    guides(shape=guide_legend(order=1, title.position = "top", title.hjust = 0.5), linetype=FALSE) +
    #theme(legend.background = element_rect(linetype="solid")) +
    theme(axis.text.x = element_text(size=16)) +
    theme(axis.title = element_text(size=14)))


pdf("Rplot_effects_f07_VarXCorpus.pdf", width = 7, height = 5)
plot.corpus
dev.off()


levs <- quantile(data$z.logWeightRatio, seq(0,1,.1))
d <- data.frame(Effect(c("Variety.Sum", "z.logWeightRatio"), f07,
                       xlevels = list(z.logWeightRatio = levs))
) %>% droplevels

(plot.weight <- ggplot(d, aes(z.logWeightRatio, fit)) + facet_wrap(~ Variety.Sum) +
  geom_hline(yintercept = 0.5, linetype="dotted") +
  geom_rect(data=subset(d, Variety.Sum=="IND"), aes(fill=Variety.Sum), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, alpha=0.05) +
  geom_rect(data=subset(d, Variety.Sum=="JA"), aes(fill=Variety.Sum), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, alpha=0.05) +
  geom_rect(data=subset(d, Variety.Sum=="IRE"), aes(fill=Variety.Sum), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, alpha=0.05) +
  scale_fill_manual(values=c(rep('gainsboro',3)), guide=FALSE) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "dimgray", alpha = .4) +
  ylim(-.3,1) +
  geom_line() +
  theme(strip.text = element_text(size=13)) +
  geom_point(data = data, aes(x = z.logWeightRatio, y = min(d$lower) - .25),
             position = position_jitter(width = .75, height = 0.05),
             size = .4, col = 'dimgray') +
  labs(x = "ratio of recipient/theme length in letters\n(log-scaled, centered and scaled by 2 SDs)", 
       y = "probability of prepositional dative")+
    theme(axis.title = element_text(size=16)))

pdf("Rplot_effects_f07_VarXWeight.pdf", width = 10, height = 7)
plot.weight
dev.off()


# -- LOOK at the DATA -----------------------

# random effects in f07

#29
data[data$ThemeHeadPlain == "try", "RecHeadPlain"]
data[data$ThemeHeadPlain == "try", "WholeConstructionPlain"]
data[data$ThemeHeadPlain == "try", "FileID"]
data[data$ThemeHeadPlain == "try", "SentencePlain"]

data[data$ThemeHeadPlain == "go", "RecHeadPlain"]
data[data$ThemeHeadPlain =="go", "WholeConstructionPlain"]
data[data$ThemeHeadPlain =="go", "FileID"]
data[data$ThemeHeadPlain =="go", "SentencePlain"]

# 13
data[data$ThemeHeadPlain == "choice", "RecHeadPlain"]
data[data$ThemeHeadPlain =="choice", "WholeConstructionPlain"]
data[data$ThemeHeadPlain =="choice", "FileID"]
data[data$ThemeHeadPlain =="choice", "SentencePlain"]
data[data$ThemeHeadPlain =="choice", "VerbSemantics"]

#14
data[data$ThemeHeadPlain == "fee", "RecHeadPlain"]
data[data$ThemeHeadPlain =="fee", "WholeConstructionPlain"]
data[data$ThemeHeadPlain =="fee", "FileID"]
data[data$ThemeHeadPlain =="fee", "SentencePlain"]
data[data$ThemeHeadPlain =="fee", "VerbSemantics"]


ranef(f07)
data[data$VerbSense=="gives.a", "FileID"]
data[data$VerbSense=="gives.t", "FileID"]

#3
data[data$ThemeHeadPlain =="birth", "RecHeadPlain"]
data[data$ThemeHeadPlain =="birth", "Resp"]
data[data$ThemeHeadPlain =="birth", "WholeConstructionPlain"]
data[data$ThemeHeadPlain =="birth", "SentencePlain"]
data[data$ThemeHeadPlain =="birth", "FileID"]

# 11
data[data$ThemeHeadPlain =="it", "WholeConstructionPlain"]
data[data$ThemeHeadPlain =="it", "Resp"]
data[data$ThemeHeadPlain =="it", "SentencePlain"]
data[data$ThemeHeadPlain =="it", "FileID"]

# 37
data[data$ThemeHeadPlain =="them", "WholeConstructionPlain"]
data[data$ThemeHeadPlain =="them", "FileID"]
data[data$ThemeHeadPlain =="them", "SentencePlain"]

# 158
data[data$ThemeHeadPlain =="attention", "WholeConstructionPlain"]
data[data$ThemeHeadPlain =="attention", "FileID"]
data[data$ThemeHeadPlain =="attention", "SentencePlain"]


# examples for main effects

# 997
data[data$RecPron=="pron" & data$ThemePron=="non-pron" & data$Corpus=="ice", "WholeConstructionPlain"]
data[data$RecPron=="pron" & data$ThemePron=="non-pron" & data$Corpus=="ice", "SentencePlain"]
data[data$RecPron=="pron" & data$ThemePron=="non-pron" & data$Corpus=="ice", "FileID"]

# 124
data[data$RecPron=="non-pron" & data$ThemePron=="pron" & data$Corpus=="ice", "SentencePlain"]
data[data$RecPron=="non-pron" & data$ThemePron=="pron" & data$Corpus=="ice", "FileID"]

# 990 
data[data$RecBinAnimacy=="animate" & data$ThemeBinAnimacy=="inanimate" & data$RecPron == "non-pron" & data$Corpus=="ice", "SentencePlain"]

# 13
data[data$RecBinAnimacy=="inanimate" & data$ThemeBinAnimacy=="animate" & data$RecPron == "non-pron" & data$Corpus=="ice", "SentencePlain"]

# 999
data[data$RecDefiniteness=="def" & data$ThemeDefiniteness=="indef" & data$RecPron == "non-pron" & data$Corpus=="ice", "FileID"]

# 218
data[data$RecDefiniteness=="indef" & data$ThemeDefiniteness=="def" & data$ThemePron == "non-pron" & data$Corpus=="ice", "SentencePlain"]






save.image("datives_GLMM.RData")