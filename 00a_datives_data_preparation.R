#=======================================================================
#   Script: 00_datives_data_preparation.R
#   Date last changed: 5 January 2017
#
#   Rversion: 3.3.2
#   Author: MR
#   This script loads the original data file and conflates levels
#   where necessary. In the end, a final dataset ("datives.txt") is
#   created which serves as input for all subsequent analyses.
#=======================================================================

#----SETUP: LIBRARIES-----------------------------------------
library(reshape2);library(plyr);library(dplyr);library(magrittr);library(lme4);library(ggplot2);library(gridExtra); library(JGmisc); library(beepr); library(effects); library(Hmisc)

setwd("/Users/melanierothlisberger/Dropbox/Universitaet/04_Statistics/statistical analyses/20170209_THESIS_full_model")
source("/Users/melanierothlisberger/Dropbox/Universitaet/03_Programming/R/helper.R")

cols <- c(do="dimgray", pd="gainsboro")

#-----SETUP: DATA------------------------
datives <- read.delim("data/20170504_datives.txt")

#--Fix mistake in the data -----------------------
datives[datives$VerbSense=="gives.a", "FileID"]
datives[datives$FileID=="IND:B:3396632", "VerbSense"] <- "give.a"
datives[datives$FileID=="IND:B:3396632", "Verb"] <- "give"

datives[datives$VerbSense=="gives.t", "FileID"]
datives[datives$FileID=="SG:B:3527930", "VerbSense"] <- "give.t"
datives[datives$FileID=="SG:B:3527930", "Verb"] <- "give"

datives[datives$VerbSense =="wish .c", "FileID"]
datives[datives$FileID=="GB:B:3068009", "VerbSense"] <- "wish.c"


# -- X-TABULATIONS ---------------------------------------------
list = list()
predictors <- names(datives)
vector <- character()
x <- 0

for (i in 1:ncol(datives)) {
  pred <- predictors[i]
  if (class(datives[[pred]])==c("numeric") || length(levels(datives[[pred]])) > 50) {
    next
  } # skip numeric predictors and those with a lot of factors
  else {
    tabs <- xtabs(~Resp+datives[[pred]], drop.unused.levels = TRUE, data=datives)
    x <- x+1
    list[[x]] <- tabs
    vector[x] <- pred # create new vector with predictors used
    
  }
}
names(list) <- vector


# -- NPExprType --------------------------------------------------------

with(datives, ftable(Variety, Resp, RecNPExprType))

(p1 = ggplot(datives, aes(RecNPExprType)) + 
  geom_bar(aes(fill=Resp), position="fill", width = 0.8)+
  scale_fill_manual(values = cols)+
  geom_hline(aes(yintercept = 0.5), linetype=2))
p1 + facet_wrap(~ Variety)

with(datives, ftable(Variety, Resp, ThemeNPExprType))

(p1 = ggplot(datives, aes(ThemeNPExprType)) + 
  geom_bar(aes(fill=Resp), position="fill", width = 0.8)+
  scale_fill_manual(values = cols)+
  geom_hline(aes(yintercept = 0.5), linetype=2))
p1 + facet_wrap(~ Variety)


# binary split pron (pp, iprn) vs. non-pron (rest)
datives$RecPron <- "non-pron"
datives[datives$RecNPExprType %in% c("pp", "iprn"), "RecPron"] <- "pron"
datives$ThemePron <- "non-pron"
datives[datives$ThemeNPExprType %in% c("pp", "iprn"), "ThemePron"] <- "pron"
datives$RecPron <- factor(datives$RecPron)
datives$ThemePron <- factor(datives$ThemePron)

# -- Nativity -----------------------------
tab.nat <- list$Nativity
prop.table(tab.nat, margin=2)
chisq.test(tab.nat) # sign.!

(p1 = ggplot(datives, aes(Nativity)) + 
  geom_bar(aes(fill=Resp), position="fill", width = 0.8)+
  scale_fill_manual(values = cols)+
  geom_hline(aes(yintercept = 0.5), linetype=2))

# -- Mode -----------------------
tab.mode <- list$Mode
prop.table(tab.mode, margin=2)
chisq.test(tab.mode)

with(datives, ftable(Variety, Resp, Mode))

(p1 = ggplot(datives, aes(Mode)) + 
  geom_bar(aes(fill=Resp), position="fill", width = 0.8)+
  scale_fill_manual(values = cols)+
  geom_hline(aes(yintercept = 0.5), linetype=2))
p1 + facet_wrap(~ Variety)

with(datives, ftable(RecPron, Resp, Mode))

(p1 = ggplot(datives, aes(Mode)) + 
  geom_bar(aes(fill=Resp), position="fill", width = 0.8)+
  scale_fill_manual(values = cols)+
  geom_hline(aes(yintercept = 0.5), linetype=2))
p1 + facet_wrap(~ RecPron)

(p1 = ggplot(datives, aes(Mode)) + 
  geom_bar(aes(fill=Resp), position="fill", width = 0.8)+
  scale_fill_manual(values = cols)+
  geom_hline(aes(yintercept = 0.5), linetype=2))
p1 + facet_wrap(~ ThemePron)

spoken <- droplevels(subset(datives, datives$Mode=="spoken"))
written <- droplevels(subset(datives, datives$Mode=="written"))

(p1 = ggplot(spoken, aes(RecPron)) + 
  geom_bar(aes(fill=Resp), position="fill", width = 0.8)+
  scale_fill_manual(values = cols)+
  geom_hline(aes(yintercept = 0.5), linetype=2))
p1 + facet_wrap(~ Variety)

(p1 = ggplot(written, aes(RecPron)) + 
  geom_bar(aes(fill=Resp), position="fill", width = 0.8)+
  scale_fill_manual(values = cols)+
  geom_hline(aes(yintercept = 0.5), linetype=2))
p1 + facet_wrap(~ Variety)

(p1 = ggplot(spoken, aes(ThemePron)) + 
  geom_bar(aes(fill=Resp), position="fill", width = 0.8)+
  scale_fill_manual(values = cols)+
  geom_hline(aes(yintercept = 0.5), linetype=2))
p1 + facet_wrap(~ Variety)

(p1 = ggplot(written, aes(ThemePron)) + 
  geom_bar(aes(fill=Resp), position="fill", width = 0.8)+
  scale_fill_manual(values = cols)+
  geom_hline(aes(yintercept = 0.5), linetype=2))
p1 + facet_wrap(~ Variety)

(p1 = ggplot(written, aes(ThemePron)) + 
  geom_bar(aes(fill=Resp), position="fill", width = 0.8)+
  scale_fill_manual(values = cols)+
  geom_hline(aes(yintercept = 0.5), linetype=2))
p1 + facet_wrap(~ RecPron)

(p1 = ggplot(datives, aes(ThemePron)) + 
  geom_bar(aes(fill=Resp), position="fill", width = 0.8)+
  scale_fill_manual(values = cols)+
  geom_hline(aes(yintercept = 0.5), linetype=2))
p1 + facet_wrap(~ RecPron+Variety)


# -- VerbSemantics -------------------------
tab.sem <- list$VerbSemantics
prop.table(tab.sem, margin=2)
with(datives, ftable(Variety, Resp, VerbSemantics))

(p1 = ggplot(datives, aes(VerbSemantics)) + 
  geom_bar(aes(fill=Resp), position="fill", width = 0.8)+
  scale_fill_manual(values = cols)+
  geom_hline(aes(yintercept = 0.5), linetype=2))
p1 + facet_wrap(~ Variety)

# Create heatmap of verbs per variety including only verbs that occur more than 1

tab.verbs <- data.frame(verbs=names(summary(datives$Verb)), count_total = summary(datives$Verb), count_pd = summary(datives[datives$Resp=="pd",]$Verb))
tab.verbs <- droplevels(subset(tab.verbs, tab.verbs$count_total>1))



# -- Animacy -------------------------------------------------
tab.rani <- list$RecAnimacy
prop.table(tab.rani, margin=2)
plot(datives$RecAnimacy)
(p1 = ggplot(datives, aes(RecAnimacy)) + 
  geom_bar(aes(fill=Resp), position="fill", width = 0.8)+
  scale_fill_manual(values = cols)+
  geom_hline(aes(yintercept = 0.5), linetype=2))
p1 + facet_wrap(~ Variety)


tab.tani <- list$ThemeAnimacy
prop.table(tab.tani, margin=2)
plot(datives$ThemeAnimacy)
(p1 = ggplot(datives, aes(ThemeAnimacy)) + 
  geom_bar(aes(fill=Resp), position="fill", width = 0.8)+
  scale_fill_manual(values = cols)+
  geom_hline(aes(yintercept = 0.5), linetype=2))
p1 + facet_wrap(~ Variety)

# Conflating levels of animacy: animate versus inanimate
datives$RecBinAnimacy <- "inanimate"
datives[datives$RecAnimacy %in% c("a1", "a2"), "RecBinAnimacy"] <- "animate"
datives$ThemeBinAnimacy <- "inanimate"
datives[datives$ThemeAnimacy %in% c("a1", "a2"), "ThemeBinAnimacy"] <- "animate"
datives$RecBinAnimacy <- factor(datives$RecBinAnimacy)
datives$ThemeBinAnimacy <- factor(datives$ThemeBinAnimacy)


# -- Complexity -----------------------------------
tab.reccox <- list$RecBinComplexity
prop.table(tab.reccox, margin=2)
plot(datives$RecBinComplexity)
(p1 = ggplot(datives, aes(RecBinComplexity)) + 
  geom_bar(aes(fill=Resp), position="fill", width = 0.8)+
  scale_fill_manual(values = cols)+
  geom_hline(aes(yintercept = 0.5), linetype=2))
p1 + facet_wrap(~ Variety)

tab.thecox <- list$ThemeBinComplexity
prop.table(tab.thecox, margin=2)
plot(datives$ThemeBinComplexity)
(p1 = ggplot(datives, aes(ThemeBinComplexity)) + 
  geom_bar(aes(fill=Resp), position="fill", width = 0.8)+
  scale_fill_manual(values = cols)+
  geom_hline(aes(yintercept = 0.5), linetype=2))
p1 + facet_wrap(~ Variety)

# -- Givenness ---------------------------
tab.recgiv <- list$RecGivenness
prop.table(tab.recgiv, margin=2)
plot(datives$RecGivenness)
(p1 = ggplot(datives, aes(RecGivenness)) + 
  geom_bar(aes(fill=Resp), position="fill", width = 0.8)+
  scale_fill_manual(values = cols)+
  geom_hline(aes(yintercept = 0.5), linetype=2))
p1 + facet_wrap(~ Variety)


# both given and new theme prefer do
tab.thegiv <- list$ThemeGivenness
prop.table(tab.thegiv, margin=2)
plot(datives$ThemeGivenness)
(p1 = ggplot(datives, aes(ThemeGivenness)) + 
  geom_bar(aes(fill=Resp), position="fill", width = 0.8)+
  scale_fill_manual(values = cols)+
  geom_hline(aes(yintercept = 0.5), linetype=2))
p1 + facet_wrap(~ Variety)

# -- Definiteness -------------------------
tab.recdef <- list$RecDefiniteness
prop.table(tab.recdef, margin=2)
plot(datives$RecDefiniteness)
(p1 = ggplot(datives, aes(RecDefiniteness)) + 
  geom_bar(aes(fill=Resp), position="fill", width = 0.8)+
  scale_fill_manual(values = cols)+
  geom_hline(aes(yintercept = 0.5), linetype=2))
p1 + facet_wrap(~ Variety)

# both def and indef theme prefer do
tab.thedef <- list$ThemeDefiniteness
prop.table(tab.thedef, margin=2)
plot(datives$ThemeDefiniteness)
(p1 = ggplot(datives, aes(ThemeDefiniteness)) + 
  geom_bar(aes(fill=Resp), position="fill", width = 0.8)+
  scale_fill_manual(values = cols)+
  geom_hline(aes(yintercept = 0.5), linetype=2))
p1 + facet_wrap(~ Variety)



# -- Persistence ------------------------

tab.pers <- list$Persistence # more persistence in do
prop.table(tab.pers, margin=2)
plot(datives$Persistence)
(p1 = ggplot(datives, aes(Persistence)) + 
  geom_bar(aes(fill=Resp), position="fill", width = 0.8)+
  scale_fill_manual(values = cols)+
  geom_hline(aes(yintercept = 0.5), linetype=2))
p1 + facet_wrap(~ Variety)


# -- Corpus/Age --------------------------
tab.corpus <- list$Corpus
prop.table(tab.corpus, margin=2)
plot(datives$Corpus)
(p1 = ggplot(datives, aes(Corpus)) + 
  geom_bar(aes(fill=Resp), position="fill", width = 0.8)+
  scale_fill_manual(values = cols)+
  geom_hline(aes(yintercept = 0.5), linetype=2))
p1 + facet_wrap(~ Variety)



# -- Length ------------------------------------------
quantile(datives$RecWordLth, probs=seq(0,1,.02)) # 98% cutoff
boxplot(datives$RecWordLth, col = "gray")
boxplot(datives$ThemeWordLth, col = "gray")
hist(datives$ThemeWordLth, col="gray", 
     breaks = length(levels(as.factor(datives$ThemeWordLth))))
hist(datives$RecWordLth, col="gray", 
     breaks=length(levels(as.factor(datives$RecWordLth))))
hist(log(datives$ThemeWordLth), col="gray", 
     breaks=length(levels(as.factor(datives$ThemeWordLth))))
hist(log(datives$RecWordLth), col="gray", 
     breaks=length(levels(as.factor(datives$RecWordLth))))

hist(datives$ThemeLetterLth, col="gray", 
     breaks = length(levels(as.factor(datives$ThemeLetterLth))))
hist(datives$RecLetterLth, col="gray", 
     breaks=length(levels(as.factor(datives$RecLetterLth))))
hist(log(datives$ThemeLetterLth), col="gray", 
     breaks=length(levels(as.factor(datives$ThemeLetterLth))))
hist(log(datives$RecLetterLth), col="gray", 
     breaks=length(levels(as.factor(datives$RecLetterLth))))
hist(log(datives$ThemeLetterLth), col="gray")
hist(log(datives$RecLetterLth), col="gray")
hist(datives[datives$Resp=="do", "RecLetterLth"], col="gray",
     breaks=length(levels(as.factor(datives$RecLetterLth))))
hist(datives[datives$Resp=="pd", "RecLetterLth"], col="gray",
     breaks=length(levels(as.factor(datives$RecLetterLth))))


# add new predictors
datives$logRecLetterLth <- log(datives$RecLetterLth)
datives$logThemeLetterLth <- log(datives$ThemeLetterLth)
datives$WeightRatio <- datives$RecLetterLth/datives$ThemeLetterLth
datives$logWeightRatio <- log(datives$WeightRatio)

# explore length
ggplot(datives, aes(x=logRecLetterLth, fill=Resp)) + 
    geom_density(alpha=.2)
ggplot(datives, aes(x=logRecLetterLth, fill=Resp)) + 
  geom_histogram(binwidth=.2, alpha=0.5, position="identity")

datives[datives$logRecLetterLth<1.2 & datives$Resp=="do", "RecPron"]
#==> mostly pronouns

ggplot(datives, aes(x=logThemeLetterLth, fill=Resp)) + 
  geom_density(alpha=.2)
ggplot(datives, aes(x=logThemeLetterLth, fill=Resp)) + 
  geom_histogram(binwidth=.2, alpha=0.5, position="identity")
datives[datives$logThemeLetterLth<1, "ThemePron"]
#==> mostly cases of "it"

ggplot(datives, aes(x=logWeightRatio, fill=Resp)) + 
  geom_density(alpha=.2) +
  facet_wrap(~Variety)
ggplot(datives, aes(x=logWeightRatio, fill=Resp)) + 
  geom_histogram(binwidth=.2, alpha=0.5, position="identity")
# much better normal distribution ==> use weight ratio!!

ggplot(datives, aes(x=Resp, y=logWeightRatio, fill=Resp)) + 
  geom_boxplot() +
  facet_wrap(~Variety)
# WeightRatio: slightly higher for PD than for DO
ggplot(datives, aes(x=Resp, y=logWeightRatio, fill=Resp)) + 
  geom_boxplot() +
  facet_wrap(~Mode)




# -- restrict dataset due to length ----------------------------
summary(datives[datives$Resp=="do", "RecWordLth"]) # longest rec = 18 words long
summary(datives[datives$Resp=="pd", "ThemeWordLth"]) # longest th = 23 words

datives2 <- droplevels(subset(datives, datives$RecWordLth <= 18)) # loose 36
datives2 <- droplevels(subset(datives2, datives2$ThemeWordLth <=23)) # loose 34

# save data
data <- datives2

# -- Others

list$RecPerson
list$ThemeConcreteness # both concrete and non-concrete themes prefer do
summary(datives$TypeTokenRatio)



# -- filter primetype according to distance from previous---------
# change NA in PrimeType to "none"
data$PrimeType <- as.character(data$PrimeType)
data$PrimeType[is.na(data$PrimeType)] <- "none" 
data$PrimeType <- as.factor(data$PrimeType)

data$PrimeTypePruned <- "none"
data$PrimeTypePruned <- as.factor(data$PrimeTypePruned)
data$PrimeTypePruned <- factor(data$PrimeTypePruned, levels=c(levels(data$PrimeTypePruned), "do", "pd"))
data$NumDistanceToPrevious <- ifelse(data$DistanceToPrevious=="none", "NA", as.vector(data$DistanceToPrevious))
data$NumDistanceToPrevious <- as.numeric(data$NumDistanceToPrevious)

data[data$NumDistanceToPrevious < 11 & !is.na(data$NumDistanceToPrevious) & data$PrimeType=="do", "PrimeTypePruned"] <- "do"
data[data$NumDistanceToPrevious < 11 & !is.na(data$NumDistanceToPrevious) & data$PrimeType=="pd", "PrimeTypePruned"] <- "pd"
data$PrimeTypePruned <- factor(data$PrimeTypePruned)

data$PersistencePruned <- ifelse(as.character(data$Resp)==as.character(data$PrimeTypePruned), "yes", ifelse(data$PrimeTypePruned=="none", "none", "no"))
data$PersistencePruned <- factor(data$PersistencePruned)




# -- standardize and scale predictors -----------------------
# NOT DONE FOR CATEGORICAL PREDICTORS
#data2 <- cbind(data, as.numeric2(data[,c(XXX)], as.char = F))
#data2 <- cbind(data2, c.(data2[, 9:11]))
#colnames(data2)[9:11] = c("NumRecAnimacy", "NumThemeDef", "NumRecPron")


data <- cbind(data, z.logRecLetterLth = z.(data[, "logRecLetterLth"], factor = 2))
data <- cbind(data, z.logThemeLetterLth = z.(data[, "logThemeLetterLth"], factor = 2))
data <- cbind(data, z.logWeightRatio = z.(data[,"logWeightRatio"], factor=2))
data <- cbind(data, z.RecHeadFreq = z.(data[,"RecHeadFreq"], factor=2))
data <- cbind(data, z.ThemeHeadFreq = z.(data[,"ThemeHeadFreq"], factor=2))
data <- cbind(data, z.RecThematicity = z.(data[,"RecThematicity"], factor=2))
data <- cbind(data, z.ThemeThematicity = z.(data[,"ThemeThematicity"], factor=2))
data <- cbind(data, z.TypeTokenRatio = z.(data[,"TypeTokenRatio"], factor=2))

# assign mean to these TTR NAs: replace NA in data with mean
summary(data$TypeTokenRatio)
data$TypeTokenRatio[is.na(data$TypeTokenRatio)] <- as.numeric(0.8160)
data$z.TypeTokenRatio <- z.(data[,"TypeTokenRatio"], factor=2)

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



# -- save the data file-----------
write.table(data, "datives.txt", quote=F, sep="\t", row.names = F)


save.image("datives_prepared.RData")
