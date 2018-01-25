#=======================================================================
#   Script: 05_datives_register_effects.R
#   Date last changed: 11 January 2018
#
#   Rversion: 3.3.2
#   Author: MR
#   This script provides the analysis of the register-specificity of
#   the English dative alternation by computing the influence of register
#   on variation and comparing probabilistic grammars across different
#   registers.
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

# -- PREPARE THE DATA --------------------------

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

# add filter to recipient and theme
data$VerbFilter <- filter.infrequent(data$Verb, threshold = 5)
data$ThemeHeadFilter <- filter.infrequent(data$ThemeHeadLemma, threshold = 5)
data$RecHeadFilter <- filter.infrequent(data$RecHeadLemma, threshold = 5)


names(data)
vars <- levels(data$Variety)


# -- REGIONAL VARIATION IN THE EFFECT OF REGISTER ----------

#--model on genre---------------
# Quick-dirty test whether logWeightRatio is interacting with Register
#formula(f <- Resp ~ z.logWeightRatio + RecPron + ThemeBinComplexity + ThemePron + z.ThemeHeadFreq)
# add contrasts to Register
data$Register.Sum <- data$GenreCoarse
contrasts(data$Register.Sum) <- contr.sum(5)

mg01 <- glmer(Resp~
                # Random factors
                (1 | VerbFilter)
              + (1 | ThemeHeadFilter) 
              + (1 | RecHeadFilter) 
              + (1 | SpeakerID)  
              + z.ThemeHeadFreq
              + ThemeBinComplexity
              + RecPron
              + ThemePron
              + z.logWeightRatio
              + Register.Sum
              + Variety.Sum
              + Variety.Sum:Register.Sum,
              data=data, family = binomial,
              glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6)))
beep(3)
summary(mg01)

fits <- fitted(mg01)
sum(ifelse(fits > .5, "pd", "do") == data$Resp)/nrow(data) # 0.9330347

somers.mer(mg01) # C=0.9790382
overdisp.mer(mg01)
collin.fnc.mer(mg01)$cnumber # 12.6
contrasts(data$Register.Sum)
contrasts(data$Variety.Sum)

library(effects)
eff.df <- as.data.frame(Effect(c("Variety.Sum", "Register.Sum"), mg01))
eff.df <- eff.df %>% droplevels

g.blue <- "#08306b" # printed
g.lightblue <- "#6baed6" # non-printed
g.green <- "#00441b" # monologue
g.lightgreen <- "#74c476" # dialogue
g.red <- '#CC0000' # online

eff.df$Register.Sum <- factor(eff.df$Register.Sum, levels=c("printed", "non-printed", "monologue", "dialogue", "online"))

eff.df$Variety.Sum <- factor(eff.df$Variety.Sum, levels=c("CAN", "GB", "IRE", "NZ", "JA", "SIN", "HK", "IND", "PHI"))


(plot2 <- ggplot(eff.df, aes(Variety.Sum, fit, group=Register.Sum)) +
    #geom_errorbar(aes(ymin = lower, ymax = upper, color = RecBinAnimacy)) +
    geom_hline(yintercept = 0.5, linetype="dashed") +
    geom_line(aes(col=Register.Sum), alpha=0.85) +
    geom_point(aes(col=Register.Sum), size=6, alpha=0.90) +
    geom_point(shape=1, size=6, color="black") +
    theme(axis.title = element_text(size=14)) +
    scale_color_manual(name="Register", values=c(g.blue, g.lightblue, g.green, g.lightgreen, g.red), guide=guide_legend(nrow=2, title.position="top", title.hjust = .5)) +
    theme(axis.title.x = element_text(size=12), axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=16), plot.title = element_text(size=16, face="bold", hjust = 0.5)) +
    labs(y="Probability of prepositional dative", x="") +
    theme(legend.justification="left", legend.position=c(0.05,0.89), legend.background = element_rect(colour="black"), legend.title=element_text(size=20), legend.text = element_text(size=16)))


# shape instead of colors
(plot3 <- ggplot(eff.df, aes(Variety.Sum, fit, group=Register.Sum)) +
    #geom_errorbar(aes(ymin = lower, ymax = upper, color = RecBinAnimacy)) +
    geom_hline(yintercept = 0.5, linetype="dashed") +
    geom_line(aes(linetype=Register.Sum), show.legend = FALSE) +
    geom_point(aes(shape=Register.Sum), size=6, fill="dimgray", alpha=0.8) +
    scale_shape_manual(name="Register", values=c(21,22,23,24,25), guide=guide_legend(nrow=2, title.position="top", title.hjust = .5)) +
    theme(axis.title.x = element_text(size=12), axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=16), plot.title = element_text(size=16, face="bold", hjust = 0.5)) +
    labs(y="Probability of prepositional dative", x="") +
    theme(legend.justification="left", legend.position=c(0.05,0.89), legend.background = element_rect(colour="black"), legend.title=element_text(size=20), legend.text = element_text(size=16)))



pdf("Rplot_genre_effects.pdf", width = 12, height = 6)
plot3
dev.off()

pdf("Rplot_genre_effects_col.pdf", width = 12, height = 6)
plot2
dev.off()

# Change coding to make unseen level visible (SinE, printed)
data$Variety.Sum2 <- factor(data$Variety, levels=c("SIN","GB","HK","IND","IRE", "JA","NZ","PHI", "CAN"))
contrasts(data$Variety.Sum2) <- contr.sum(9)

data$Register.Sum2 <- factor(data$GenreCoarse, levels=c("printed","dialogue", "monologue","non-printed","online"))
contrasts(data$Register.Sum2) <- contr.sum(5)

mg02 <- glmer(Resp~
                # Random factors
                (1 | VerbFilter)
              + (1 | ThemeHeadFilter) 
              + (1 | RecHeadFilter) 
              + (1 | SpeakerID)  
              + z.ThemeHeadFreq
              + ThemeBinComplexity
              + RecPron
              + ThemePron
              + z.logWeightRatio
              + Register.Sum2
              + Variety.Sum2
              + Variety.Sum2:Register.Sum2,
              data=data, family = binomial,
              glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6)))
beep(3)
summary(mg02)

contrasts(data$Variety.Sum2)


# -- random forest for relative ranking by variety ---------------------
library(party)
# formula
f.crf <- Resp ~
  GenreCoarse +
  VerbSemantics +
  z.logWeightRatio +
  RecGivenness +
  ThemeGivenness +
  RecDefiniteness +
  ThemeDefiniteness +
  z.RecHeadFreq +
  z.ThemeHeadFreq +
  z.RecThematicity +
  z.ThemeThematicity +
  PrimeTypePruned +
  RecPron +
  ThemePron +
  RecBinAnimacy +
  ThemeBinAnimacy +
  z.TypeTokenRatio

forest.controls = cforest_unbiased(ntree=3000, mtry=5)

t1 <- proc.time()
for (i in 1:9){
  cat("Working on", vars[i], "...\n")
  d <- droplevels(subset(data, data$Variety == vars[i]))
  
  # rescale predictors
  d$z.logWeightRatio = z.(d[,"logWeightRatio"], factor=2)
  d$z.RecHeadFreq = z.(d[,"RecHeadFreq"], factor=2)
  d$z.ThemeHeadFreq = z.(d[,"ThemeHeadFreq"], factor=2)
  d$z.RecThematicity = z.(d[,"RecThematicity"], factor=2)
  d$z.ThemeThematicity = z.(d[,"ThemeThematicity"], factor=2)
  d$z.TypeTokenRatio = z.(d[,"TypeTokenRatio"], factor=2)
  
  # run model
  rf <- cforest(f.crf, data = d, controls = forest.controls)
  assign(paste(vars[i], "rf", sep = '.'), rf)
  
  # add predictions
  d$preds <- predict(rf)
  name <- paste(vars[i], "df", sep= '.')
  assign(name,d)
  
  # get varimp
  varimp <- party::varimpAUC(rf)
  assign(paste(vars[i], "varimp", sep = '.'), varimp)
  
  # create dataframe with varimp
  parameters <- c(as.vector(varimp))
  if (i == 1) {
    pars <- data.frame(row.names = names(varimp))
    pars <- cbind(pars, parameters)
  }
  else {pars <- cbind(pars, parameters)}
  rm("d")
  rm("varimp")
  rm("rf")
}

t2 <- proc.time() - t1

names(pars) <- vars

pars_sum <- transform(pars, sum=rowSums(pars))
pars.sorted <- pars[order(pars_sum[,"sum"], decreasing = F),]

# rename factors
rownames(pars.sorted) <- c("ThemeAnimacy", "TTR", "PrimeType", "RecDefiniteness", "ThemeThematicity", "RecThematicity", "Register", "ThemeGivenness", "VerbSemantics","RecAnimacy", "ThemeDefiniteness", "ThemePron", "RecGivenness","ThemeHeadFreq" , "RecHeadFreq" ,"RecPron", "WeightRatio")

pars_long <- melt(t(pars.sorted))
names(pars_long) <- c("Variety", "Factor", "Varimp")



(plot.rf.varimp.VAR <- ggplot(pars_long, aes(Factor,Varimp)) +
    geom_bar(fill="dimgray", stat="identity", color="black", width=.6) +
    coord_flip() +
    theme(axis.text = element_text(size=10, colour = "black")) +
    theme(title = element_text(size=12)) +
    labs(x="", y="Conditional variable importance") +
    facet_wrap(~Variety) +
    # add the Shih line
    geom_hline(yintercept=abs(min(pars_long$Varimp)),linetype="dashed") )

pdf("Rplot_cox_crf_varimp_byvar.pdf", width = 9, height=10)
plot.rf.varimp.VAR
dev.off()



#--INTRA-SYSTEMIC VARIATION------------------------
f2 <- formula(f2<- Resp ~ (1 | VerbFilter) + (1 | ThemeHeadFilter) + (1 | RecHeadFilter) + Register.Sum*(z.ThemeHeadFreq + ThemeBinComplexity + RecPron + ThemePron + z.logWeightRatio))

byvars_models <- list()
for (i in 1:length(vars)) {
  cat("Modeling", vars[i], "...")
  v <- subset(data, Variety==vars[i]) %>% droplevels
  
  cat("Standardising and setting contrasts...")
  v$z.logWeightRatio = z.(v[,"logWeightRatio"], factor=2)
  v$z.ThemeHeadFreq = z.(v[,"ThemeHeadFreq"], factor=2)
  v$ThemeHeadFilter <- filter.infrequent(v$ThemeHeadLemma, threshold=8)
  v$RecHeadFilter <- filter.infrequent(v$RecHeadLemma, threshold=8)
  v$VerbFilter <- filter.infrequent(v$Verb, threshold=8)
  contrasts(v$Register.Sum) <- contr.sum(5)
  
  cat("fitting model", "\n")
  m <- glmer(f2, data=v, family = binomial,
             glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6)))
  byvars_models[[i]] <- m
}

names(byvars_models) <- vars

# check the models
summary(byvars_models$CAN)
ggLogit.plot(byvars_models$CAN, data=subset(data, Variety=="CAN"))
overdisp.mer(byvars_models$CAN)
Anova(byvars_models$CAN) # recpron & theme cox

summary(byvars_models$GB)
Anova(byvars_models$GB) # themepron
Anova(byvars_models$HK) # thcox in dialogue
Anova(byvars_models$IRE) # not converged
summary(byvars_models$IND)
Anova(byvars_models$IND) # recpron in online, less PD with non-pron; recpron in dia with more PD; thcox in online
summary(byvars_models$JA) #
Anova(byvars_models$JA) # WR in non-printed
Anova(byvars_models$NZ) # recpron in online
summary(byvars_models$NZ)
summary(byvars_models$PHI)
Anova(byvars_models$PHI) # less PD if non-pron recpron in online
Anova(byvars_models$SIN) # none
summary(byvars_models$SIN)


overdisp.mer(byvars_models$CAN)
overdisp.mer(byvars_models$GB)
overdisp.mer(byvars_models$HK)
overdisp.mer(byvars_models$IRE)
overdisp.mer(byvars_models$IND)
overdisp.mer(byvars_models$JA) # !
overdisp.mer(byvars_models$NZ)
overdisp.mer(byvars_models$PHI)
overdisp.mer(byvars_models$SIN)
collin.fnc.mer(byvars_models$CAN)$cnumber # 
collin.fnc.mer(byvars_models$GB)$cnumber # 
collin.fnc.mer(byvars_models$HK)$cnumber # 
collin.fnc.mer(byvars_models$IND)$cnumber # 
collin.fnc.mer(byvars_models$IRE)$cnumber # 
collin.fnc.mer(byvars_models$JA)$cnumber # 
collin.fnc.mer(byvars_models$NZ)$cnumber # 
collin.fnc.mer(byvars_models$PHI)$cnumber # 
collin.fnc.mer(byvars_models$SIN)$cnumber # 

# ==> all quite bad

# run same analysis with different register and variety levels
# Rerun model to check for printed

f3 <- formula(f3<- Resp ~ (1 | VerbFilter) + (1 | ThemeHeadFilter) + (1 | RecHeadFilter) + Register.Sum2*(z.ThemeHeadFreq + ThemeBinComplexity + RecPron + ThemePron + z.logWeightRatio))

byvars_models2 <- list()
for (i in 1:length(vars)) {
  cat("Modeling", vars[i], "...")
  d <- subset(data, Variety==vars[i]) %>% droplevels
  
  cat("Standardising and setting contrasts...")
  d$z.logWeightRatio = z.(d[,"logWeightRatio"], factor=2)
  d$z.ThemeHeadFreq = z.(d[,"ThemeHeadFreq"], factor=2)
  d$ThemeHeadFilter <- filter.infrequent(d$ThemeHeadLemma, threshold=8)
  d$RecHeadFilter <- filter.infrequent(d$RecHeadLemma, threshold=8)
  d$VerbFilter <- filter.infrequent(d$Verb, threshold=8)
  contrasts(d$Register.Sum2) <- contr.sum(5)
  
  cat("fitting model", "\n")
  m <- glmer(f3, data=d, family = binomial,
             glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6)))
  byvars_models2[[i]] <- m
}

names(byvars_models2) <- vars

# check the models
summary(byvars_models2$CAN)
ggLogit.plot(byvars_models2$CAN, data=subset(data, Variety=="CAN"))
overdisp.mer(byvars_models2$CAN)
Anova(byvars_models2$CAN) # recpron in general

summary(byvars_models2$GB)
Anova(byvars_models2$GB) # recpron in printed, themepron in printed & dialogue
Anova(byvars_models2$HK) # thcox
Anova(byvars_models2$IRE) # not converged
Anova(byvars_models2$IND) # recpron in dialogue
Anova(byvars_models2$JA) # less PD in printed theme pron; less PD WR online
Anova(byvars_models2$NZ) # more PD in printed recpron
Anova(byvars_models2$PHI) # none
Anova(byvars_models2$SIN) # none


overdisp.mer(byvars_models2$CAN)
overdisp.mer(byvars_models2$GB)
overdisp.mer(byvars_models2$HK)
overdisp.mer(byvars_models2$IRE)
overdisp.mer(byvars_models2$IND)
overdisp.mer(byvars_models2$JA) # !
overdisp.mer(byvars_models2$NZ)
overdisp.mer(byvars_models2$PHI)
overdisp.mer(byvars_models2$SIN)
collin.fnc.mer(byvars_models2$CAN)$cnumber # 15.6
collin.fnc.mer(byvars_models2$GB)$cnumber # 16.4
collin.fnc.mer(byvars_models2$HK)$cnumber # 17.5
collin.fnc.mer(byvars_models2$IND)$cnumber # 16.9
collin.fnc.mer(byvars_models2$IRE)$cnumber # 18.2
collin.fnc.mer(byvars_models2$JA)$cnumber # 15.5
collin.fnc.mer(byvars_models2$NZ)$cnumber # 16.3
collin.fnc.mer(byvars_models2$PHI)$cnumber # 17.6
collin.fnc.mer(byvars_models2$SIN)$cnumber # 14.0



#--..effects plots of by-variety models for recpron-------------------

# GB
effs.gb <- Effect(c("RecPron","Register.Sum"), mod = byvars_models$GB) %>% as.data.frame

# IND
effs.ind <- Effect(c("RecPron","Register.Sum"), mod = byvars_models$IND) %>% as.data.frame

# NZ
effs.nz <- Effect(c("RecPron","Register.Sum"), mod = byvars_models$NZ) %>% as.data.frame

# PHI
effs.phi <- Effect(c("RecPron","Register.Sum"), mod = byvars_models$PHI) %>% as.data.frame


# create dataframe with all data incl. variety and then facet_grid(Variety ~ Register.Sum)
effs.all <- rbind(effs.gb, effs.ind, effs.nz, effs.phi)
effs.all$Variety <- c(rep(c("BrE", "IndE", "NZE", "PhiE"), each=10))

cols <- c("BrE" = "gainsboro", "CanE"="gainsboro", "IndE"="gainsboro", "JamE"="gainsboro", "PhiE"="gainsboro")


(plot4 <- ggplot(effs.all, aes(RecPron, fit)) +
  geom_line(aes(group = Register.Sum), linetype = 1, col="black") +
  geom_errorbar(aes(ymin = lower, ymax = upper), col = "black", width=.1, size=0.35) +
  geom_hline(yintercept = 0.5, linetype="dashed") +
  geom_point(col="black", size=2)+
  facet_grid(Variety ~ Register.Sum) +
  theme(strip.text.y = element_text(colour="black")) +
  theme(strip.text = element_text(size=12)) +
  geom_rect(data=subset(effs.all, Variety=="BrE" & Register.Sum =="printed"), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, alpha=0.05) +
  geom_rect(data=subset(effs.all, Variety=="IndE" & Register.Sum %in% c("online", "dialogue")), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, alpha=0.05) +
  geom_rect(data=subset(effs.all, Variety=="NZE" & Register.Sum %in% c("online", "printed")), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, alpha=0.05) +
  geom_rect(data=subset(effs.all, Variety=="PhiE" & Register.Sum=="online"), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, alpha=0.05) +
  scale_fill_manual(values=cols) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=15)) +
  theme(axis.text.x = element_text(size=16)) +
  labs(x = "\n recipient pronominality", y = "probability of prepositional"))


pdf("Rplot_effects_genre_var.pdf", width = 10, height = 8)
plot4
dev.off()

#--..effects plots of by-variety models for themepron-------------------
# GB
effs.gb <- Effect(c("ThemePron","Register.Sum"), mod = byvars_models$GB) %>% as.data.frame

# JA
effs.ja <- Effect(c("ThemePron","Register.Sum"), mod = byvars_models$JA) %>% as.data.frame


# create dataframe with all data incl. variety and then facet_grid(Variety ~ Register.Sum)
effs.all <- rbind(effs.gb, effs.ja)
effs.all$Variety <- c(rep(c("BrE", "JamE"), each=10))

cols <- c("BrE" = "gainsboro", "JamE"="gainsboro")


(plot5 <- ggplot(effs.all, aes(ThemePron, fit)) +
    geom_line(aes(group = Register.Sum), linetype = 1, col="black") +
    geom_errorbar(aes(ymin = lower, ymax = upper), col = "black", width=.1, size=0.35) +
    geom_hline(yintercept = 0.5, linetype="dashed") +
    geom_point(col="black", size=2)+
    facet_grid(Variety ~ Register.Sum) +
    theme(strip.text.y = element_text(colour="black")) +
    theme(strip.text = element_text(size=12)) +
    geom_rect(data=subset(effs.all, Variety=="BrE" & Register.Sum %in% c("printed", "dialogue")), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, alpha=0.05) +
    geom_rect(data=subset(effs.all, Variety=="JamE" & Register.Sum == "printed"), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, alpha=0.05) +
    scale_fill_manual(values=cols) +
    theme(legend.position = "none") +
    theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=16)) +
    labs(x = "\n recipient pronominality", y = "probability of prepositional"))


pdf("Rplot_effects_genre_var_themepron.pdf", width = 10, height = 8)
plot5
dev.off()

# reason for reverse effect of theme pron:

data[data$Variety=="GB" & data$Resp=="pd" & data$ThemePron == "non-pron" & data$GenreCoarse=="printed" & data$RecPron=="pron", "FileID"]

data[data$Variety=="JA" & data$Resp=="do" & data$ThemePron == "pron" & data$GenreCoarse=="printed", "FileID"]




save.image("datives_register.RData")