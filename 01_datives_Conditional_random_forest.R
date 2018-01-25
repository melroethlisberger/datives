#=======================================================================
#   Script: 01_datives_Conditional_random_forest.R
#   Date last changed: 5 January 2018
#
#   Rversion: 3.3.2
#   Author: MR
#
#   This script provides the Rcode for the random forest that was run
#   on the full dataset for the dissertation.
#=======================================================================

#----SETUP: LIBRARIES-----------------------------------------
library(reshape2);library(plyr);library(dplyr);library(magrittr);library(lme4);library(ggplot2);library(gridExtra); library(JGmisc); library(beepr); library(effects); library(Hmisc)

setwd("/Users/melanierothlisberger/Dropbox/Universitaet/04_Statistics/statistical analyses/20170209_THESIS_full_model")
source("/Users/melanierothlisberger/Dropbox/Universitaet/03_Programming/R/helper.R")

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

# -- Constraints to be included -------------------------------------
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
# PrimeTypePruned
# RecPron
# ThemePron
# logWeightRatio
# VerbSemantics
# Corpus
# Variety


# -- CTREES -------------------------------------------------------------
library(partykit)
set.seed(123456)
control = ctree_control(mincriterion = 0.998)
# bonferroni corrected alpha level, 0.05/22 (predictors)
# tree with all factors
tree1 <- ctree(Resp ~ 
                 Variety +  
                 Corpus +
                 Mode +
                 VerbSemantics +
                 z.logWeightRatio +
                 RecGivenness +
                 ThemeGivenness +
                 RecDefiniteness +
                 ThemeDefiniteness +
                 RecBinComplexity +
                 ThemeBinComplexity +
                 z.RecHeadFreq +
                 z.ThemeHeadFreq +
                 z.RecThematicity +
                 z.ThemeThematicity +
                 ThemeConcreteness +
                 PrimeTypePruned +
                 RecPron +
                 ThemePron +
                 RecBinAnimacy +
                 ThemeBinAnimacy +
                 z.TypeTokenRatio
               , data=data, control=control)
plot(tree1, gp = gpar(fontsize = 8)) 

# different seed
set.seed(784378)
tree2 <- ctree(Resp ~ 
                 Variety +  
                 Corpus +
                 Mode +
                 VerbSemantics +
                 z.logWeightRatio +
                 RecGivenness +
                 ThemeGivenness +
                 RecDefiniteness +
                 ThemeDefiniteness +
                 RecBinComplexity +
                 ThemeBinComplexity +
                 z.RecHeadFreq +
                 z.ThemeHeadFreq +
                 z.RecThematicity +
                 z.ThemeThematicity +
                 ThemeConcreteness +
                 PrimeTypePruned +
                 RecPron +
                 ThemePron +
                 RecBinAnimacy +
                 ThemeBinAnimacy +
                 z.TypeTokenRatio
               , data=data, control=control)
plot(tree2, gp = gpar(fontsize = 8)) 

control = ctree_control(mincriterion = 0.99)
tree3 <- ctree(Resp ~ 
                 Variety +  
                 Corpus +
                 Mode +
                 VerbSemantics +
                 z.logWeightRatio +
                 RecGivenness +
                 ThemeGivenness +
                 RecDefiniteness +
                 ThemeDefiniteness +
                 RecBinComplexity +
                 ThemeBinComplexity +
                 z.RecHeadFreq +
                 z.ThemeHeadFreq +
                 z.RecThematicity +
                 z.ThemeThematicity +
                 ThemeConcreteness +
                 PrimeTypePruned +
                 RecPron +
                 ThemePron +
                 RecBinAnimacy +
                 ThemeBinAnimacy +
                 z.TypeTokenRatio
               , data=data, control=control)
plot(tree3, gp = gpar(fontsize = 8)) 

# -- RANDOM FORESTS ------------------------------------------------------
detach("package:partykit", unload=TRUE)
library(party)

forest.controls = cforest_control(ntree=2000, mtry=3)
set.seed(12355)

rf1 <- cforest(Resp ~
                 Variety + 
                 Corpus + 
                 Mode +
                 VerbSemantics +
                 z.logRecLetterLth +
                 z.logThemeLetterLth +
                 RecGivenness +
                 ThemeGivenness +
                 RecDefiniteness +
                 ThemeDefiniteness +
                 RecBinComplexity +
                 ThemeBinComplexity +
                 z.RecHeadFreq +
                 z.ThemeHeadFreq +
                 z.RecThematicity +
                 z.ThemeThematicity +
                 PrimeTypePruned +
                 RecPron +
                 ThemePron +
                 RecBinAnimacy +
                 ThemeBinAnimacy +
                 z.TypeTokenRatio,
               data = data, controls = forest.controls)
varimp.rf1 = varimp(rf1, conditional=FALSE)
dotplot(sort(varimp.rf1))

# same forest with different random seed
set.seed(98638)
rf2 <- cforest(Resp ~
                 Variety + 
                 Corpus + 
                 Mode +
                 VerbSemantics +
                 z.logWeightRatio +
                 RecGivenness +
                 ThemeGivenness +
                 RecDefiniteness +
                 ThemeDefiniteness +
                 RecBinComplexity +
                 ThemeBinComplexity +
                 z.RecHeadFreq +
                 z.ThemeHeadFreq +
                 z.RecThematicity +
                 z.ThemeThematicity +
                 PrimeTypePruned +
                 RecPron +
                 ThemePron +
                 RecBinAnimacy +
                 ThemeBinAnimacy +
                 z.TypeTokenRatio,
               data = data, controls = forest.controls)
varimp.rf2 = varimp(rf2, conditional=FALSE)
dotplot(sort(varimp.rf2))
varimpAUC.rf2 <- party::varimpAUC(rf2)
dotplot(sort(varimp.rf2))

# -- train with caret package -----------------------------------------------
library(caret)
seeds <- vector(mode="list", length=nrow(data) +1)
seeds <- lapply(seeds, function(x) 1:20)

fitControl <- trainControl(
  method="cv",
  number=10,
  returnResamp = "all",
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  seeds = seeds)

fitControl2 <- trainControl(
  method="repeatedcv",
  number=10,
  repeats =3,
  search = "random",
  returnResamp = "all",
  classProbs = TRUE,
  summaryFunction = twoClassSummary)



set.seed(849)

fit.rf1 <- train(x=data[,c("Variety", "Corpus", "Mode", "VerbSemantics", "z.logWeightRatio", "RecGivenness", "ThemeGivenness", "RecDefiniteness", "ThemeDefiniteness", "RecBinComplexity", "ThemeBinComplexity", "z.RecHeadFreq", "z.ThemeHeadFreq", "z.RecThematicity", "z.ThemeThematicity", "PrimeTypePruned", "RecPron", "ThemePron", "RecBinAnimacy", "ThemeBinAnimacy", "z.TypeTokenRatio")], y=data[,"Resp"], method="cforest", trControl = fitControl2, controls=cforest_unbiased(ntree=500))

# -- FINAL RANDOM FOREST --------------------------------------------
#mtry = 5
forest.controls = cforest_unbiased(ntree=2000, mtry=5)
set.seed(12355)

rf3 <- cforest(Resp ~
                 Variety + 
                 Corpus + 
                 Mode +
                 VerbSemantics +
                 z.logWeightRatio +
                 RecGivenness +
                 ThemeGivenness +
                 RecDefiniteness +
                 ThemeDefiniteness +
                 RecBinComplexity +
                 ThemeBinComplexity +
                 z.RecHeadFreq +
                 z.ThemeHeadFreq +
                 z.RecThematicity +
                 z.ThemeThematicity +
                 PrimeTypePruned +
                 RecPron +
                 ThemePron +
                 RecBinAnimacy +
                 ThemeBinAnimacy +
                 z.TypeTokenRatio,
               data = data, controls = forest.controls)
# varimp.rf3 = varimp(rf3, conditional=TRUE) doesn't work
varimpAUC.rf3 <- party::varimpAUC(rf3)


#-- Evaluation of final forest ---------------------

probs <- treeresponse(rf3)
# Add predictions for the prepositional dative variant to dataframe
data$probs.rf3 <- sapply(probs, FUN = function(x) return(x[1])) %>%
  as.vector
somers2(1-data$probs.rf3, as.numeric(data$Resp)-1)
data$preds.rf3 <- predict(rf3)
(acc.rf3 <- sum(data$Resp==data$preds.rf3)/nrow(data)) # 89.2%

#-----..one-tailed binomial test to compare accuracy to baseline------
# Assuming you have a vector of predictions, "preds", and a baseline probability, "base.p", You can do it like so - testing the null hypothesis that the baseline is greater than the model's accuracy:
# baseline = 8987 do, 4184 pd = 13171 == 73.84% for do
hits <- sum(data$preds == data$Resp)
binom.test(x = hits, n = nrow(data), p = 0.7384, alternative = "greater")

# -- verify crf with different seed --------------------------------

forest.controls = cforest_unbiased(ntree=2000, mtry=5)
set.seed(657849)

rf3b <- cforest(Resp ~
                  Variety + 
                  Corpus + 
                  Mode +
                  VerbSemantics +
                  z.logWeightRatio +
                  RecGivenness +
                  ThemeGivenness +
                  RecDefiniteness +
                  ThemeDefiniteness +
                  RecBinComplexity +
                  ThemeBinComplexity +
                  z.RecHeadFreq +
                  z.ThemeHeadFreq +
                  z.RecThematicity +
                  z.ThemeThematicity +
                  PrimeTypePruned +
                  RecPron +
                  ThemePron +
                  RecBinAnimacy +
                  ThemeBinAnimacy +
                  z.TypeTokenRatio,
                data = data, controls = forest.controls)
# varimp.rf3 = varimp(rf3, conditional=TRUE) doesn't work
varimpAUC.rf3b <- party::varimpAUC(rf3b)
dot1 <- dotplot(sort(varimpAUC.rf3))
dot2 <- dotplot(sort(varimpAUC.rf3b))
require(gridExtra)
grid.arrange(dot1, dot2, ncol=2)
# LOOK THE SAME!!


# -- plot variable importance ----------------------------------------
varimpAUC.rf3 <- sort(varimpAUC.rf3)
varimp.all <- data.frame(vals = varimpAUC.rf3, 
                         preds = names(varimpAUC.rf3),
                         levels = names(sort(varimpAUC.rf3)),
                         var = rep("Variable Importance", 21))

# change names of predictors:
varimp.all$preds <- c("ThemeAnimacy", "Corpus", "TTR", "RecThematicity", "ThemeThematicity", "RecDefiniteness", "PrimeType", "Mode", "Variety", "ThemeGivenness", "RecComplexity", "RecAnimacy", "VerbSemantics", "RecGivenness", "ThemeDefiniteness", "RecHeadFreq", "ThemeHeadFreq", "ThemePron", "ThemeComplexity", "RecPron", "WeightRatio")
varimp.all$levels <- c("ThemeAnimacy", "Corpus", "TTR", "RecThematicity", "ThemeThematicity", "RecDefiniteness", "PrimeType", "Mode", "Variety", "ThemeGivenness", "RecComplexity", "RecAnimacy", "VerbSemantics", "RecGivenness", "ThemeDefiniteness", "RecHeadFreq", "ThemeHeadFreq", "ThemePron", "ThemeComplexity", "RecPron", "WeightRatio")
varimp.all$preds <- factor(varimp.all$preds, as.character(varimp.all$preds))



(plot.varimp <- ggplot(varimp.all, aes(preds,vals)) +
  geom_bar(fill="dimgray", stat="identity", color="black", width=.6) +
  coord_flip() +
  theme(axis.text = element_text(size=16, colour = "black")) +
  theme(title = element_text(size=20)) +
  labs(x="", y="Conditional variable importance") +
  # add the Shih line
  geom_hline(yintercept=abs(min(varimp.all$vals)),linetype="dashed"))

pdf("Rplot_crf_varimp_all.pdf", width = 12, height = 8)
plot.varimp
dev.off()


# -- BY-VARIETY RANDOM FORESTS -----------------------------------

f1.crf <- Resp ~
  Variety +  
  Mode +
  Corpus + 
  VerbSemantics +
  logWeightRatio +
  RecGivenness +
  ThemeGivenness +
  RecDefiniteness +
  ThemeDefiniteness +
  RecBinComplexity +
  ThemeBinComplexity +
  RecHeadFreq +
  ThemeHeadFreq +
  RecThematicity +
  ThemeThematicity +
  PrimeTypePruned +
  RecPron +
  ThemePron +
  RecBinAnimacy +
  ThemeBinAnimacy +
  TypeTokenRatio


data2 <- data[ ,all.vars(f1.crf)]
data2 <- cbind(data2, as.numeric2(data2[,c(3:5,7:12,17:21)], as.char = F))
data2 <- cbind(data2, c.(data2[, 23:36]))
colnames(data2)[23:36] = c("NumMode", "NumCorpus", "NumVerbSemantics", "NumRecGivenness", "NumThemeGivenness", "NumRecDefiniteness", "NumThemeDefiniteness", "NumRecComplexity", "NumThemeComplexity", "NumPrimeType", "NumRecPron", "NumThemePron", "NumRecAnimacy", "NumThemeAnimacy")

vars <- levels(data$Variety)

# formula
f2.crf <- Resp ~
  c.Mode +
  c.Corpus + 
  c.VerbSemantics +
  z.logWeightRatio +
  c.RecGivenness +
  c.ThemeGivenness +
  c.RecDefiniteness +
  c.ThemeDefiniteness +
  c.RecBinComplexity +
  c.ThemeBinComplexity +
  z.RecHeadFreq +
  z.ThemeHeadFreq +
  z.RecThematicity +
  z.ThemeThematicity +
  c.PrimeTypePruned +
  c.RecPron +
  c.ThemePron +
  c.RecBinAnimacy +
  c.ThemeBinAnimacy +
  z.TypeTokenRatio

forest.controls = cforest_unbiased(ntree=2000, mtry=5)
# add -x to make order reversed
t1 <- proc.time()
for (i in 1:9){
  
  d <- droplevels(subset(data2, data2$Variety == vars[i]))
  
  # rescale predictors
  d$z.logWeightRatio = z.(d[,6], factor=2)
  d$z.RecHeadFreq = z.(d[,13], factor=2)
  d$z.ThemeHeadFreq = z.(d[,14], factor=2)
  d$z.RecThematicity = z.(d[,15], factor=2)
  d$z.ThemeThematicity = z.(d[,16], factor=2)
  d$z.TypeTokenRatio = z.(d[,22], factor=2)
  d$c.Mode = c.(d[,23])
  d$c.Corpus = c.(d[,24])
  d$c.VerbSemantics = c.(d[,25])
  d$c.RecGivenness = c.(d[,26])
  d$c.ThemeGivenness = c.(d[,27])
  d$c.RecDefiniteness = c.(d[,28])
  d$c.ThemeDefiniteness = c.(d[,29])
  d$c.RecComplexity = c.(d[,30])
  d$c.ThemeComplexity = c.(d[,31])
  d$c.PrimeTypePruned = c.(d[,32])
  d$c.RecPron = c.(d[,33])
  d$c.ThemePron = c.(d[,34])
  d$c.RecAnimacy = c.(d[,35])
  d$c.ThemeAnimacy = c.(d[,36])
  
  
  # run model
  rf <- cforest(f2.crf, data = d, controls = forest.controls)
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
  varimp <- as.data.frame(varimp)
  names(varimp) <- c("Chisq")
  rank.df <- transform(varimp, 
                    Pred.rank = ave(Chisq,FUN = function(x) rank(-x, ties.method = "first")))
  rank <- c(as.vector(rank.df$Pred.rank))
  if (i == 1) {
    pars.crf.full <- data.frame(row.names = rownames(varimp))
    pars.crf.full <- cbind(pars.crf.full, parameters, rank)
    pars.crf.chisq <- data.frame(row.names = rownames(varimp))
    pars.crf.chisq <- cbind(pars.crf.chisq, parameters)
    pars.crf.rk <- data.frame(row.names = rownames(varimp))
    pars.crf.rk <- cbind(pars.crf.rk, rank)
  }
  else {
    pars.crf.full <- cbind(pars.crf.full, parameters, rank)
    pars.crf.chisq <- cbind(pars.crf.chisq, parameters)
    pars.crf.rk <- cbind(pars.crf.rk, rank)
    }
  rm("d")
  rm("varimp")
  rm("rf")
}

t2 <- proc.time() - t1

names(pars.crf.chisq) <- vars



# transform pars.crf to long data format

crf_long <- melt(t(pars.crf.chisq))
names(crf_long) <- c("Variety", "Factor", "Varimp")


# -- .. plot variable importances --------
# create dataframe with ranking included
test <- pars.crf.full
names(test) <- rep(vars,times=1, each=2)
names(pars.crf.full) <- rep(vars,times=1, each=2)

pars1 <- melt(t(pars.crf.full[,seq(from=1, to=ncol(pars.crf.full), by=2)]))
pars2 <- melt(t(pars.crf.full[,seq(from=2, to=ncol(pars.crf.full), by=2)]))
pars.crf.full_long <- cbind(pars1, pars2[,3])
names(pars.crf.full_long) <- c("Variety", "Factor", "Varimp", "Rank")


# order predictors by size
crf.sum <- transform(pars.crf.chisq, sum=rowSums(pars.crf.chisq))
pars.sorted <- pars.crf.chisq[order(crf.sum[,"sum"], decreasing = F),]

# rename factors
rownames(pars.sorted) <- c("ThemeAnimacy", "PrimeType", "TTR", "Corpus", "Mode", "RecDefiniteness", "RecThematicity", "ThemeThematicity", "ThemeGivenness", "VerbSemantics", "RecAnimacy", "ThemeDefiniteness", "RecComplexity", "ThemePron","ThemeHeadFreq" , "RecGivenness", "ThemeComplexity" ,  "RecHeadFreq", "RecPron", "WeightRatio")

pars_long <- melt(t(pars.sorted))
names(pars_long) <- c("Variety", "Factor", "Varimp")



(plot.crf.varimp.VAR <- ggplot(pars_long, aes(Factor,Varimp)) +
    geom_bar(fill="dimgray", stat="identity", color="black", width=.6) +
    coord_flip() +
    theme(axis.text = element_text(size=10, colour = "black")) +
    theme(title = element_text(size=12)) +
    labs(x="", y="Conditional variable importance") +
    facet_wrap(~Variety) +
    # add the Shih line
    geom_hline(yintercept=abs(min(varimp.all$vals)),linetype="dashed") )

pdf("Rplot_crf_varimp_byvar.pdf", width = 10, height=9)
plot.crf.varimp.VAR
dev.off()



# -- .. barplot of accuracy and C-stats -----------------------------
# every data frame df contains a column preds with predictions.
# create dataframe with all predictions/accuracy
require(party)

forest_list = list()
forest_list[[1]] <- CAN.rf
forest_list[[2]] <- GB.rf
forest_list[[3]] <- HK.rf
forest_list[[4]] <- IND.rf
forest_list[[5]] <- IRE.rf
forest_list[[6]] <- JA.rf
forest_list[[7]] <- NZ.rf
forest_list[[8]] <- PHI.rf
forest_list[[9]] <- SIN.rf

data_list = list()
data_list[[1]] <- CAN.df
data_list[[2]] <- GB.df
data_list[[3]] <- HK.df
data_list[[4]] <- IND.df
data_list[[5]] <- IRE.df
data_list[[6]] <- JA.df
data_list[[7]] <- NZ.df
data_list[[8]] <- PHI.df
data_list[[9]] <- SIN.df


for(i in 1:9){
  print(paste("Working on", vars[i], "..."))
  df <- data_list[[i]]
  probs <- treeresponse(forest_list[[i]])
  probs <- sapply(probs, FUN = function(x) return(x[1])) %>% as.vector
  somers <- somers2(1-probs, as.numeric(df$Resp)-1)
  somers <- somers[1]
  acc <- sum(df$Resp==df$preds)/nrow(df)
  
  if(i == 1){
  eval <- as.data.frame(row.names=vars[i], cbind(somers, acc))
  }
  else{
  row <- as.data.frame(row.names=vars[i], cbind(somers,acc))
  eval <- rbind(eval,row)
  }
  rm("probs")
}

eval$Variety <- vars

eval$Variety <- with(eval, reorder(Variety, acc)) # Sorting

require(ggthemes)
(ggplot(eval, aes(Variety, acc))
+ geom_bar(stat="identity")
+ geom_text(aes(label=Variety), hjust=1.5, size=4, color="white")
+ geom_text(aes(label=paste(round(100*acc,1), "%", sep=" ")), 
            hjust=-0.125, size=4, color="black")
+ labs(title="Accuracy of by-variety random forests")
+ theme(axis.ticks.y = element_blank(),
        axis.text.y  = element_blank(),
        plot.title = element_text(size=20, face="bold", hjust=.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),            
        axis.text.x  = element_text(size=12),
        strip.text   = element_text(size=15, face="bold"))
+ scale_y_continuous(label=function(x){paste(100 * x, "%", sep="")},
                     limits=c(0,1.15), 
                     breaks=c(0, 0.25, 0.5, 0.75, 1))
+ scale_fill_economist()
+ coord_flip()
)






# -- .. plot observed vs. predicted -----------------------------------
# plot proportion of observed vs predicted for prepositional variant

Varieties <- c("CAN", "GB", "HK", "IRE", "IND", "JA", "NZ", "PHI", "SIN")

predictions <- as.vector(c(sum(CAN.df$preds=="pd")/nrow(CAN.df), sum(GB.df$preds=="pd")/nrow(GB.df), sum(HK.df$preds=="pd")/nrow(HK.df), sum(IRE.df$preds=="pd")/nrow(IRE.df), sum(IND.df$preds=="pd")/nrow(IND.df), sum(JA.df$preds=="pd")/nrow(JA.df), sum(NZ.df$preds=="pd")/nrow(NZ.df), sum(PHI.df$preds=="pd")/nrow(PHI.df), sum(SIN.df$preds=="pd")/nrow(SIN.df)))

observed <- as.vector(c(sum(CAN.df$Resp=="pd")/nrow(CAN.df), sum(GB.df$Resp=="pd")/nrow(GB.df), sum(HK.df$Resp=="pd")/nrow(HK.df), sum(IRE.df$Resp=="pd")/nrow(IRE.df), sum(IND.df$Resp=="pd")/nrow(IND.df), sum(JA.df$Resp=="pd")/nrow(JA.df), sum(NZ.df$Resp=="pd")/nrow(NZ.df), sum(PHI.df$Resp=="pd")/nrow(PHI.df), sum(SIN.df$Resp=="pd")/nrow(SIN.df)))

obspreds <- as.data.frame(cbind(Varieties, round(predictions,5), round(observed,5)))
names(obspreds) <- c("Varieties", "predictions", "observations")
obspreds$predictions <- as.numeric(as.character(obspreds$predictions))
obspreds$observations <- as.numeric(as.character(obspreds$observations))

require(ggrepel)
set.seed(42)
(pred.pd <- ggplot(obspreds, aes(observations, predictions, label = Varieties)) +
  xlim(0.15,0.5) +
  ylim(0.15,0.5) +
  geom_point(fill="dimgray", size=3, color="black", pch=21) +
  geom_text_repel(size=7) +
  #geom_text_repel(aes(label=verb, size=4.5))
  geom_abline(intercept = 0, slope=1, linetype=2) +
  labs(x="proportion of observed prepositional datives", y="proportion of predicted prepositional datives") +
  theme(axis.title = element_text(size=18)))

pdf("Rplot_observedVSpredicted_crf.pdf", width=8, height=6)
pred.pd
dev.off()


# -- .. and for ditransitives --------------
predictions2 <- as.vector(c(sum(CAN.df$preds=="do")/nrow(CAN.df), sum(GB.df$preds=="do")/nrow(GB.df), sum(HK.df$preds=="do")/nrow(HK.df), sum(IRE.df$preds=="do")/nrow(IRE.df), sum(IND.df$preds=="do")/nrow(IND.df), sum(JA.df$preds=="do")/nrow(JA.df), sum(NZ.df$preds=="do")/nrow(NZ.df), sum(PHI.df$preds=="do")/nrow(PHI.df), sum(SIN.df$preds=="do")/nrow(SIN.df)))

observed2 <- as.vector(c(sum(CAN.df$Resp=="do")/nrow(CAN.df), sum(GB.df$Resp=="do")/nrow(GB.df), sum(HK.df$Resp=="do")/nrow(HK.df), sum(IRE.df$Resp=="do")/nrow(IRE.df), sum(IND.df$Resp=="do")/nrow(IND.df), sum(JA.df$Resp=="do")/nrow(JA.df), sum(NZ.df$Resp=="do")/nrow(NZ.df), sum(PHI.df$Resp=="do")/nrow(PHI.df), sum(SIN.df$Resp=="do")/nrow(SIN.df)))

obspreds2 <- as.data.frame(cbind(Varieties, round(predictions2,5), round(observed2,5)))
names(obspreds2) <- c("Varieties", "predictions", "observations")
obspreds2$predictions <- as.numeric(as.character(obspreds2$predictions))
obspreds2$observations <- as.numeric(as.character(obspreds2$observations))

(pred.do <- ggplot(obspreds2, aes(observations, predictions, label = Varieties)) +
  xlim(0.5,0.8) +
  ylim(0.5,0.8) +
  geom_point(color=qlvlblue, size=3) +
  geom_text(hjust=0.6, vjust=-1, size=7) +
  geom_abline(intercept = 0, slope=1, linetype=2) +
  labs(x="proportion of observed ditransitive datives", y="proportion of predicted ditransitive datives") +
  theme(axis.title = element_text(size=18)))

require(gridExtra)
grid.arrange(pred.pd, pred.do, ncol=2)

# -- .. confusion matrix ---------------------------------------------
library(caret)
require(caret)
x <- with(CAN.df, table(Resp, preds))
x <- with(GB.df, table(Resp, preds))
x <- with(HK.df, table(Resp, preds))
x <- with(IND.df, table(Resp, preds))
x <- with(IRE.df, table(Resp, preds))
x <- with(JA.df, table(Resp, preds))
x <- with(NZ.df, table(Resp, preds))
x <- with(PHI.df, table(Resp, preds))
x <- with(SIN.df, table(Resp, preds))


# matrix for the whole forest
y <- with(data, table(preds.rf3, Resp))
confusionMatrix(y)
confusionMatrix(x)


save.image("datives_CRF.RData")
