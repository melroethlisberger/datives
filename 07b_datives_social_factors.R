#=======================================================================
#   Script: 07b_datives_social_factors.R
#   Date last changed: 12 January 2018
#
#   Rversion: 3.3.2
#   Author: MR
#   This script provides the code for the analysis of the metadata in 
#   the discussion section of the dissertation.
#=======================================================================

#----SETUP: LIBRARIES-----------------------------------------
library(reshape2);library(plyr);library(dplyr);library(magrittr);library(lme4);library(ggplot2);library(gridExtra); library(JasonG); library(beepr); library(effects); library(Hmisc)

setwd("/Users/melanierothlisberger/Dropbox/Universitaet/04_Statistics/statistical analyses/20170209_THESIS_full_model")
source("/Users/melanierothlisberger/Dropbox/Universitaet/03_Programming/R/helper.R")


#---SETUP: PLOTTING----------------
theme_mr = theme_set(theme_light())
theme_mr = theme_update(axis.text = element_text(size = rel(1.0), color="black"),
                        axis.ticks = element_line(colour = "grey90", size = 0.25),
                        axis.title = element_text(size=rel(0.9)),
                        panel.border = element_rect(color = "black"),
                        strip.background=element_rect(fill="grey95", color="black"), 
                        strip.text.x=element_text(color="black"))


#-----SETUP: DATA ------------------------
#datives <- read.delim("data/20170818_datives.txt") # for creation of new dataset
# no further level conflating necessary for following datafile:
datives <- read.delim("OA/datives_social_spoken.txt") # spoken texts excl. SIN

# binary split pron (pp, iprn) vs. non-pron (rest)
datives$RecPron <- "non-pron"
datives[datives$RecNPExprType %in% c("pp", "iprn"), "RecPron"] <- "pron"
datives$ThemePron <- "non-pron"
datives[datives$ThemeNPExprType %in% c("pp", "iprn"), "ThemePron"] <- "pron"
datives$RecPron <- factor(datives$RecPron)
datives$ThemePron <- factor(datives$ThemePron)

# Conflating levels of animacy: animate versus inanimate
datives$RecBinAnimacy <- "inanimate"
datives[datives$RecAnimacy %in% c("a1", "a2"), "RecBinAnimacy"] <- "animate"
datives$ThemeBinAnimacy <- "inanimate"
datives[datives$ThemeAnimacy %in% c("a1", "a2"), "ThemeBinAnimacy"] <- "animate"
datives$RecBinAnimacy <- factor(datives$RecBinAnimacy)
datives$ThemeBinAnimacy <- factor(datives$ThemeBinAnimacy)

# add new predictors
datives$logRecLetterLth <- log(datives$RecLetterLth)
datives$logThemeLetterLth <- log(datives$ThemeLetterLth)
datives$WeightRatio <- datives$RecLetterLth/datives$ThemeLetterLth
datives$logWeightRatio <- log(datives$WeightRatio)

# -- restrict dataset due to length ----------------------------
summary(datives[datives$Resp=="do", "RecWordLth"]) # longest rec = 18 words long
summary(datives[datives$Resp=="pd", "ThemeWordLth"]) # longest th = 23 words

datives2 <- droplevels(subset(datives, datives$RecWordLth <= 18)) # loose 36
datives2 <- droplevels(subset(datives2, datives2$ThemeWordLth <=23)) # loose 34

data <- datives2


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


# -- standardize and scale predictors -----------------------

data <- cbind(data, z.logRecLetterLth = z.(data[, "logRecLetterLth"], factor = 2))
data <- cbind(data, z.logThemeLetterLth = z.(data[, "logThemeLetterLth"], factor = 2))
data <- cbind(data, z.logWeightRatio = z.(data[, "logWeightRatio"], factor=2))
data <- cbind(data, z.RecHeadFreq = z.(data[, 'RecHeadFreq'], factor=2))
data <- cbind(data, z.ThemeHeadFreq = z.(data[, 'ThemeHeadFreq'], factor=2))
data <- cbind(data, z.RecThematicity = z.(data[, 'RecThematicity'], factor=2))
data <- cbind(data, z.ThemeThematicity = z.(data[, 'ThemeThematicity'], factor=2))
data <- cbind(data, z.TypeTokenRatio = z.(data[, 'TypeTokenRatio'], factor=2))



# -- TREES: exploratory ------------------------------------
library(partykit)
set.seed(123456)
control = ctree_control(mincriterion = 0.95)

# tree with all external factors
tree1 <- ctree(Resp ~ 
                 Variety +  
                 Nativity +
                 Corpus +
                 Mode +
                 MD_Gender +
                 MD_Age +
                 MD_TextTime +
                 MD_Education +
                 MD_Occupation +
                 MD_FirstLanguage
                 , data=data, control=control)


plot(tree1, gp = gpar(fontsize = 8)) 

# leave out time related factors to see whether age turns up
tree2 <- ctree(Resp ~ 
                 Variety +  
                 Nativity +
                 #Corpus +
                 Mode +
                 MD_Gender +
                 MD_Age +
                 #MD_TextTime +
                 MD_Education +
                 MD_Occupation +
                 MD_FirstLanguage
               , data=data, control=control)


plot(tree2, gp = gpar(fontsize = 8)) 



# -- X-TABULATIONS ---------------------------------------------
list = list()
predictors <- names(data)
vector <- character()
x <- 0

for (i in 1:ncol(data)) {
  pred <- predictors[i]
  if (class(data[[pred]])==c("numeric") || length(levels(data[[pred]])) > 50) {
    next
  } # skip numeric predictors and those with a lot of factors
  else {
    tabs <- xtabs(~Resp+data[[pred]], drop.unused.levels = TRUE, data=data)
    x <- x+1
    list[[x]] <- tabs
    vector[x] <- pred # create new vector with predictors used
    
  }
}
names(list) <- vector

cols <- c(do="dimgray",pd="gainsboro")
# -- Age --------------------------------------------------------

with(data, ftable(Variety, Resp, MD_Age))

(p1 = ggplot(data, aes(MD_Age)) + 
  geom_bar(aes(fill=Resp), position="fill", width = 0.8)+
  scale_fill_manual(values = cols)+
  geom_hline(aes(yintercept = 0.5), linetype=2))
p1 + facet_wrap(~ Variety)

# New conflated coding for age:
levels(data$MD_Age)
str(data)
barplot(prop.table(table(data$MD_Age)))

# age combinations due to data sparseness: 60+, 0-19, 20-29, 30-39, 40-49, 50-59. If numbers are on the border, I took the bin with more years in, e.g. 19-24 -> 20-30

data$Age6 <- "teenies"
data[data$MD_Age %in% c("(?)26-30", "18-25", "19-24", "19-25", "20", "20-24", "21", "21-25", "22", "23", "24", "25", "25-29", "25-30", "26", "26-30", "26-33", "27", "28", "29"), "Age6"] <- "twenties"
data[data$MD_Age %in% c("26-45", "26-45?", "30-34", "31", "31-35", "31-40", "32", "34", "34-41", "35-39", "36", "36-40", "37"), "Age6"] <- "thirties"
data[data$MD_Age %in% c("40", "40-44", "40+", "41", "41-45", "41-46", "41-50", "42", "42-49", "43", "44", "45", "45-49", "45-50", "46", "46-50", "47", "48", "49"), "Age6"] <- "fourties"
data[data$MD_Age %in% c("50", "50-54", "50+", "51", "51-55", "46-65", "51-60", "53", "54" ,"55", "55-59", "56", "56-60", "57", "58"), "Age6"] <- "fifties"
data[data$MD_Age %in% c("60", "60-64", "60plus", "61", "61+", "63", "64", "65", "65-69", "65+", "66+", "68", "70-74", "71", "73", "75-79", "80", "85-89", "92"), "Age6"] <- "sixtiesPlus"

data$Age6 <- factor(data$Age6)
data$Age6 <- factor(data$Age6, levels=c("teenies", "twenties", "thirties", "fourties", "fifties", "sixtiesPlus"))

# -- Mode -----------------------
# restrict data to spoken data only
spoken <- droplevels(subset(data, data$Mode=="spoken"))

# -- New X-Tab with spoken only --------------------
list = list()
predictors <- names(spoken)
vector <- character()
x <- 0

for (i in 1:ncol(spoken)) {
  pred <- predictors[i]
  if (class(spoken[[pred]])==c("numeric") || length(levels(spoken[[pred]])) > 50) {
    next
  } # skip numeric predictors and those with a lot of factors
  else {
    tabs <- xtabs(~Resp+spoken[[pred]], drop.unused.levels = TRUE, data=spoken)
    x <- x+1
    list[[x]] <- tabs
    vector[x] <- pred # create new vector with predictors used
    
  }
}
names(list) <- vector


# -- Gender ----------------------
prop.table(list$MD_Gender, margin=2) # more pds with males

# -- Corpus/Age --------------------------
prop.table(list$Corpus, margin=2) # only ICE data

# -- TextTime ----------------------------
tab.tt <- as.data.frame(prop.table(list$MD_TextTime, margin=2))
names(tab.tt) <- c("Resp", "Time", "Freq")
do.tt <- droplevels(subset(tab.tt, tab.tt$Resp=="do"))


(p1 = ggplot(do.tt, aes(x=Time, y=Freq, group=1)) + 
    geom_point() +
    geom_line() +
    geom_smooth(method="lm") +
    geom_hline(aes(yintercept = 0.5), linetype=2))


# increase of the ditransitive dative over time roughly
list$MD_TextTime # consider conflating some levels..

# -- Age6 ---------
list$Age6
prop.table(list$Age6, margin=2) # indication of increase of ditransitive over time as well


# -- FirstLanguage -----------
levels(spoken$MD_FirstLanguage)

# add new level with English agains the rest
spoken$FirstLanguageBin <- "Other"
spoken[spoken$MD_FirstLanguage %in% c("English"), "FirstLanguageBin"] <- "English"
spoken$FirstLanguageBin <- as.factor(spoken$FirstLanguageBin)


# -- Occupation -------------
levels(spoken$MD_Occupation) # too many levels

# -- Education ------------
levels(spoken$MD_Education) # too many levels, could be conflated into University, Secondary school, ...?

# -- final tree with conflated levels ------------
# leave out SIN since no info given there:
spoken2 <- droplevels(subset(spoken, spoken$Variety != "SIN"))

tree.final <- ctree(Resp ~ 
                      Age6 +
                      FirstLanguageBin + 
                      MD_Gender +
                      Variety +
                      Nativity +
                      GenreFine
                      ,data=spoken2)

plot(tree.final, gp = gpar(fontsize = 6))

tree.give <- ctree(Resp ~ 
                      Age6 +
                      FirstLanguageBin +
                      MD_Gender +
                      Variety +
                      Nativity +
                      GenreFine
                    ,data=spoken2[spoken2$Verb=="give", ])

plot(tree.give, gp = gpar(fontsize = 8))


write.table(spoken2, "datives_social_spoken.txt", quote=F, sep="\t", row.names = F)
save.image("datives_meta_descriptive.RData")
