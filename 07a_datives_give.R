#=======================================================================
#   Script: 07a_datives_give.R
#   Date last changed: 11 January 2018
#
#   Rversion: 3.3.2
#   Author: MR
#   This script provides the code used for the comparison of a model
#   reduced to the verb 'give' with the full model.
#=======================================================================

#----SETUP: LIBRARIES-----------------------------------------
library(reshape2);library(plyr);library(dplyr);library(magrittr);library(lme4);library(ggplot2);library(gridExtra); library(JGmisc); library(beepr); library(effects); library(Hmisc); library(car)

setwd("/Users/melanierothlisberger/Dropbox/Universitaet/04_Statistics/statistical analyses/20170209_THESIS_full_model")
source("/Users/melanierothlisberger/Dropbox/Universitaet/03_Programming/R/helper.R")

options(scipen = 999)

#---SETUP: PLOTTING----------------
theme_mr = theme_set(theme_light())
theme_mr = theme_update(axis.text = element_text(size = rel(1.0), color="black"),
                        axis.ticks = element_line(colour = "grey90", size = 0.25),
                        axis.title = element_text(size=rel(0.9)),
                        panel.border = element_rect(color = "black"),
                        strip.background=element_rect(fill="grey95", color="black"), 
                        strip.text.x=element_text(color="black"))

#-----SETUP: DATA------------------------
data <- read.delim("OA/datives.txt")

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
data$Persistence <- relevel(data$Persistence, ref = "none")
data$RecPron <- relevel(data$RecPron, ref="pron")
data$ThemePron <- relevel(data$ThemePron, ref="non-pron")
data$PrimeType <- relevel(data$PrimeType, ref="none")
data$PrimeTypePruned <- relevel(data$PrimeTypePruned, ref="none")

# relevel Variety
data$Variety <- factor(data$Variety, levels=c("CAN", "GB", "HK", "IND", "IRE", "JA", "NZ", "PHI", "SIN"))
data$Variety.Sum <- data$Variety
contrasts(data$Variety.Sum) = contr.sum(9)


#----------------GIVE ONLY-------------
#-----------.. replication of Bresnan and Hay (2008)----------
# Run analysis with restriction to GIVE verbs only, add Semantics of verb to fixed effect structure and follow Bresnan and Hay 2008:
# factors = log length of rec and theme, pronominality of rec and theme, givenness of rec and theme (theme in interaction with verb semantics), animacy of recipient (and in interaction with variety), semantic class, variety
# coding: pronouns incl. demonstratives => Pron3, 
# coding: length = log length of constituent in words

# 1. Restrict dataset to give only
datives3 <- droplevels(subset(data, Verb=="give")) # 7030 instances

# 2. Centralize and standardize based on the new dataset, add reference levels and contrast coding
datives3$VarietySum <- datives3$Variety
contrasts(datives3$VarietySum) = contr.sum(9)

datives3$z.logRecLetterLth <- z.(datives3$logRecLetterLth, factor=2)
datives3$z.logThemeLetterLth <- z.(datives3$logThemeLetterLth, factor=2)
datives3$z.logWeightRatio <- z.(datives3$logWeightRatio, factor=2)

datives3$Resp <- relevel(datives3$Resp, ref = "do")
datives3$Nativity <- relevel(datives3$Nativity, ref = "L1")
datives3$Mode <- relevel(datives3$Mode, ref = "spoken")
datives3$Corpus <- relevel(datives3$Corpus, ref = "ice")
datives3$VerbSemantics <- relevel(datives3$VerbSemantics, ref = "a")
datives3$RecGivenness <- relevel(datives3$RecGivenness, ref = "given")
datives3$ThemeGivenness <- relevel(datives3$ThemeGivenness, ref = "new")
datives3$RecDefiniteness <- relevel(datives3$RecDefiniteness, ref = "def")
datives3$ThemeDefiniteness <- relevel(datives3$ThemeDefiniteness, ref = "indef")
datives3$RecBinComplexity <- relevel(datives3$RecBinComplexity, ref = "simple")
datives3$ThemeBinComplexity <- relevel(datives3$ThemeBinComplexity, ref = "complex")
datives3$RecBinAnimacy <- relevel(datives3$RecBinAnimacy, ref = "animate")
datives3$ThemeBinAnimacy <- relevel(datives3$ThemeBinAnimacy, ref = "inanimate")
datives3$RecPerson <- relevel(datives3$RecPerson, ref = "local")
datives3$ThemeConcreteness <- relevel(datives3$ThemeConcreteness, ref="non-concrete")
datives3$Persistence <- relevel(datives3$Persistence, ref = "none")
datives3$RecPron <- relevel(datives3$RecPron, ref="pron")
datives3$ThemePron <- relevel(datives3$ThemePron, ref="non-pron")
datives3$PrimeType <- relevel(datives3$PrimeType, ref="none")
datives3$PrimeTypePruned <- relevel(datives3$PrimeTypePruned, ref="none")



# 3. Run model including random effect of speaker

give.a <- glmer(Resp ~ (1|SpeakerID) +
                  RecGivenness +
                  RecPron +
                  ThemePron +
                  z.logRecLetterLth +
                  z.logThemeLetterLth +
                  (ThemeGivenness * VerbSemantics) +
                  (RecBinAnimacy * VarietySum),
                data = datives3, family = binomial,
                control = glmerControl(optimizer="bobyqa",
                                       optCtrl = list(maxfun=1e6)))
summary(give.a)
Anova(give.a)
summary(give.a)$coefficients 
# - Diagnostics -
somers.mer(give.a) 
ggLogit.plot(give.a, datives3) # S-shaped curve -> not good
collin.fnc.mer(give.a)$cnumber # 7.9
kappa.mer(give.a) # 7.5
as.data.frame(sort(vif.mer(give.a))) # 10 for rec length
overdisp.mer(give.a) # but no indication of overdispersion

# Run model to check for SIN (VarSum)
datives3$Var <- NULL
datives3$Var <- datives3$Variety
datives3$Var <- relevel(datives3$Var, ref="SIN")

datives3$VarSum <- datives3$Var
contrasts(datives3$VarSum) = contr.sum(9)

give.b <- glmer(Resp ~ (1|SpeakerID) +
                  RecGivenness +
                  RecPron +
                  ThemePron +
                  z.logRecLetterLth +
                  z.logThemeLetterLth +
                  (ThemeGivenness * VerbSemantics) +
                  (RecBinAnimacy * VarSum),
                data = datives3, family = binomial,
                control = glmerControl(optimizer="bobyqa",
                                       optCtrl = list(maxfun=1e6)))

summary(give.b)


#----..plot B&H GIVE model-------
cols <- c(do=kublue,pd=brewerblue)

eff <- Effect(c("VarietySum", "RecBinAnimacy"), give.a)
plot(eff, type="link", zvar="Variety", multiline = TRUE, ylab = "log odds of prepositional dative", xlab = "Variety")

# as ggplot
df.eff <- Effect(c("RecBinAnimacy", "VarietySum"), give.a) %>% as.data.frame

df.eff$VarietySum <- factor(df.eff$VarietySum, levels=c("CAN", "GB", "IRE", "NZ","JA", "SIN", "HK", "IND", "PHI"))

# print 4x5 inches # use shape instead of colour = 
(p1 <- ggplot(df.eff, aes(VarietySum, fit, group=RecBinAnimacy, shape=RecBinAnimacy)) +
  #geom_errorbar(aes(ymin = lower, ymax = upper, color = RecBinAnimacy)) +
  geom_line(aes(linetype=RecBinAnimacy)) +
  geom_point(size=3) +
  #geom_hline(yintercept = 0.5, linetype="dashed") +
  scale_linetype_discrete(name="recipient animacy", breaks = c("animate", "inanimate"), labels = c("animate", "inanimate")) +
  scale_shape_discrete(name="recipient animacy", breaks = c("animate", "inanimate"), labels = c("animate", "inanimate")) +
  labs(y="probability of prepositional dative", x="give") +
  theme(legend.position = "bottom") +
  theme(axis.text.x = element_text(size=16), axis.title = element_text(size=16), axis.title.x = element_text(face = "italic")))


#-------.. plot FINAL model-------
load("../datives_GLMM.RData") # for the f07 model
eff1 <- Effect(c("Variety.Sum", "RecBinAnimacy"), f07)
plot(eff1, type="link", zvar="Variety", multiline = TRUE, ylab = "log odds of prepositional dative", xlab = "Variety")

df.eff1 <- Effect(c("RecBinAnimacy", "Variety.Sum"), f07) %>% as.data.frame
df.eff1$Variety.Sum <- factor(df.eff1$Variety.Sum, levels=c("CAN", "GB", "IRE", "NZ","JA", "SIN", "HK", "IND", "PHI"))

(p2 <- ggplot(df.eff1, aes(Variety.Sum, fit, group=RecBinAnimacy, shape=RecBinAnimacy)) +
  geom_line(aes(linetype=RecBinAnimacy)) +
  geom_point(size=3) +
  scale_linetype_discrete(name="recipient animacy", breaks = c("animate", "inanimate"), labels = c("animate", "inanimate"), guide=FALSE) +
  scale_shape_discrete(name="recipient animacy", breaks = c("animate", "inanimate"), labels = c("animate", "inanimate")) +
  labs(y="probability of prepositional dative", x="all verbs") +
    guides(shape=guide_legend(ncol=2,title.position = "top", title.hjust = 0.5))+
    theme(legend.position = "bottom",legend.title = element_text(size=18,face = "bold"), legend.text = element_text(size=18)) +
  theme(axis.text.x = element_text(size=16), axis.title = element_text(size=16)))

# with ggplot:
require(gridExtra)
mylegend <- g_legend(p2)
bothplots <- grid.arrange(arrangeGrob(p1 + theme(legend.position = "none"),
                         p2 + theme(legend.position = "none") + ylab(""),
                         nrow=1), mylegend, nrow=2, heights=c(10,1))

pdf("Rplot_give_effects_animacy.pdf", width = 10, height = 6)
grid.arrange(arrangeGrob(p1 + theme(legend.position = "none"),
                         p2 + theme(legend.position = "none") + ylab(""),
                         nrow=1), mylegend, nrow=2, heights=c(10,1))
dev.off()




#-- Comparison to other studies---------------------------------------
# Szmrecsanyi et al. (forthcoming)

datives3$VerbSemanticsSum <- datives3$VerbSemantics
contrasts(datives3$VerbSemanticsSum) <- contr.sum(3)

give.c <- glmer(Resp ~ (1|FileID) +
                  RecGivenness +
                  ThemePron +
                  ThemeGivenness + 
                  VarietySum * 
                  (RecPron +
                  z.logWeightRatio +
                  VerbSemanticsSum +
                  RecBinAnimacy),
                data = datives3, family = binomial,
                control = glmerControl(optimizer="bobyqa",
                                       optCtrl = list(maxfun=1e6)))

summary(give.c)

canbre.df <- subset(datives3, datives3$Variety %in% c("CAN", "GB")) %>% droplevels
canbre.df <- subset(canbre.df, canbre.df$Register=="SpokInf") %>% droplevels
brenz.df <- subset(datives3, datives3$Variety %in% c("NZ", "GB")) %>% droplevels
brenz.df <- subset(brenz.df, brenz.df$Register=="SpokInf") %>% droplevels
cannz.df <- subset(datives3, datives3$Variety %in% c("CAN", "NZ")) %>% droplevels
cannz.df <- subset(cannz.df, cannz.df$Register=="SpokInf") %>% droplevels

cb01 <- glm(Resp ~ z.logWeightRatio + Variety + Variety:z.logWeightRatio, family=binomial, data=canbre.df)
cb02 <- glm(Resp ~ RecPron + Variety + Variety:RecPron, family=binomial, data=canbre.df)
summary(cb01)
summary(cb02)

bz01 <- glm(Resp ~ VerbSemantics + Variety + Variety:VerbSemantics, family=binomial, data= brenz.df)
summary(bz01)
            
cz01 <- glm(Resp ~ z.logWeightRatio + Variety + Variety:z.logWeightRatio, family=binomial, data= cannz.df)              
summary(cz01)

# n.s. differences in spoken informal speech

#--..other effects plot------------------------

# RecPron
f <- Effect(c("VarietySum", "RecPron"), give.a)
plot(f, type="link", zvar="Variety", multiline = TRUE, ylab = "log odds of prepositional dative", xlab = "Variety")



save.image("datives_give.RData")