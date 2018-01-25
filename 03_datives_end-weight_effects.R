#=======================================================================
#   Script: 03_datives_end-weight_effects.R
#   Date last changed: 9 January 2018
#
#   Rversion: 3.3.2
#   Author: MR
#   This script provides the Rcode for the discussion of complexity and
#   end-weight effects in the dissertation.
#=======================================================================

#----SETUP: LIBRARIES-----------------------------------------
library(reshape2);library(plyr);library(dplyr);library(magrittr);library(lme4);library(ggplot2);library(gridExtra); library(JGmisc); library(beepr); library(effects); library(Hmisc); library(car)

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

cols <- c(do="dimgray",pd="gainsboro")

# -- Set the reference levels and the contrasts -------------------------------
# predictions are for the prepositional dative, i.e. all levels are set for the prototypical ditransitive dative
data$Resp <- relevel(data$Resp, ref = "do")
data$Nativity <- relevel(data$Nativity, ref = "L1")
data$Mode <- relevel(data$Mode, ref = "spoken")
data$Corpus <- relevel(data$Corpus, ref = "glowbe")
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

data$Variety.Sum <- data$Variety
contrasts(data$Variety.Sum) = contr.sum(9)


# -- Filter infrequent ranef levels -------------------------------
# go for the 90% threshold
data$RecHeadFilter <- filter.infrequent(data$RecHeadLemma, threshold = 4)
data$ThemeHeadFilter <- filter.infrequent(data$ThemeHeadLemma, threshold = 8)


# -- reduce dataset to ICE ---------------------------------------
ice <- droplevels(subset(data, data$Corpus=="ice"))
contrasts(ice$Variety.Sum) = contr.sum(9)

quantile(as.vector(table(ice$RecHeadLemma)), probs=seq(0,1,0.05))
quantile(as.vector(table(ice$ThemeHeadLemma)), probs=seq(0,1,0.05))
ice$RecHeadFilter <- filter.infrequent(ice$RecHeadLemma, threshold = 4)
ice$ThemeHeadFilter <- filter.infrequent(ice$ThemeHeadLemma, threshold = 9)




# -- EXPLORE LEVELS OF COMPLEXITY  ---------------------------------
# -- conflate levels where necessary due to sparseness --------------

with(ice, ftable(Variety, Resp, RecComplexity))
(plice1 = ggplot(ice, aes(RecComplexity)) + 
  geom_bar(aes(fill=Resp), position="fill", width = 0.8)+
  scale_fill_manual(values = cols)+
  geom_hline(aes(yintercept = 0.5), linetype=2))

# I should probably conflate single postmodifications into spp and svp but keep genitives (gn) as "s"
ice$RecComplexityConflated <- "s"
ice[ice$RecComplexity %in% c("co", "ge", "postad", "pp"), "RecComplexityConflated"] <- "spp"
ice[ice$RecComplexity %in% c("nonfin", "advc", "rc", "cp"), "RecComplexityConflated"] <- "svp"
ice[ice$RecComplexity %in% c("tpp", "mpp"), "RecComplexityConflated"] <- "mpp"
ice[ice$RecComplexity %in% c("tvp", "mvp"), "RecComplexityConflated"] <- "mvp"
ice$RecComplexityConflated <- factor(ice$RecComplexityConflated)

with(ice, ftable(Variety, Resp, RecComplexityConflated))

ice$ThemeComplexityConflated <- "s"
ice[ice$ThemeComplexity %in% c("co", "ge", "postad", "pp"), "ThemeComplexityConflated"] <- "spp"
ice[ice$ThemeComplexity %in% c("nonfin", "advc", "rc", "cp"), "ThemeComplexityConflated"] <- "svp"
ice[ice$ThemeComplexity %in% c("tpp", "mpp"), "ThemeComplexityConflated"] <- "mpp"
ice[ice$ThemeComplexity %in% c("tvp", "mvp"), "ThemeComplexityConflated"] <- "mvp"
ice$ThemeComplexityConflated <- factor(ice$ThemeComplexityConflated)
with(ice, ftable(Variety, Resp, ThemeComplexityConflated))

# -- set reference levels for new predictors ---------------
ice$ThemeComplexityConflated <- relevel(ice$ThemeComplexityConflated, ref = "s")
ice$RecComplexityConflated <- relevel(ice$RecComplexityConflated, ref = "s")


# -- Bar chart of frequencies ---------------------------
cols <- c(do='dimgray',pd='gainsboro')
ice$Variety <- factor(ice$Variety, levels=c("CAN", "GB", "HK", "IND", "IRE", "JA", "NZ", "PHI", "SIN"))

# -- .. all Complexity: Rec --------------
frequencies = prop.table(with(ice, table(RecComplexityConflated, Resp)),1) # create prop table
frequencies <- data.frame(frequencies)
tokens = with(ice, table(RecComplexityConflated, Resp)) # create raw freq table
tokens <- data.frame(tokens)
tokens_sorted <- arrange(tokens, RecComplexityConflated, Resp) # sort both using library(plyr)
freq_sorted <- arrange(frequencies, RecComplexityConflated, Resp)

# replace every second percentage by 1 to use as geom_text
freq_sorted[1,3] <- 1
freq_sorted[3,3] <- 1
freq_sorted[5,3] <- 1
freq_sorted[7,3] <- 1
freq_sorted[9,3] <- 1


freq_sorted$Tokens <- tokens_sorted$Freq  #add raw frequency to dataset
levels(ice$RecComplexityConflated)
ice$RecComplexityConflated <- factor(ice$RecComplexityConflated, levels = c("s", "spp", "svp", "mpp", "mvp"))



(plot.reccox <- ggplot(ice, aes(RecComplexityConflated, fill=Resp)) + 
  geom_bar(aes(fill=Resp), position="fill",width = 0.6, color="black") + 
  labs(title = "", y="proportion of tokens", x="RecComplexity") + 
    theme(axis.text.x = element_text(size=16)) +
    theme(axis.title.y = element_text(size=16)) +
    theme(axis.title.x = element_text(size=16)) +
  geom_text(data = freq_sorted, aes(x=RecComplexityConflated,y=Freq, label = Tokens), vjust=1.5, color=ifelse(freq_sorted$Resp=="pd", "black", "white"), size=6) + 
  scale_fill_manual(values = cols, name="Variant", labels=c("ditransitive", "prepositional")) +
  geom_hline(yintercept = 0.5, colour="black", linetype=2) +
    guides(fill=guide_legend(ncol=2,title.position = "top", title.hjust = 0.5))+
    theme(legend.position = "bottom",legend.title = element_text(size=16,face = "bold"), legend.text = element_text(size=16)))

pdf("Rplot_cox_barplot_reccox.pdf", width = 6, height = 8)
plot.reccox
dev.off()


# -- .. all Complexity: Theme --------------
frequencies = prop.table(with(ice, table(ThemeComplexityConflated, Resp)),1) # create prop table
frequencies <- data.frame(frequencies)
tokens = with(ice, table(ThemeComplexityConflated, Resp)) # create raw freq table
tokens <- data.frame(tokens)
tokens_sorted <- arrange(tokens, ThemeComplexityConflated, Resp) # sort both using library(plyr)
freq_sorted <- arrange(frequencies, ThemeComplexityConflated, Resp)

# replace every second percentage by 1 to use as geom_text
freq_sorted[1,3] <- 1
freq_sorted[3,3] <- 1
freq_sorted[5,3] <- 1
freq_sorted[7,3] <- 1
freq_sorted[9,3] <- 1

# make sure that 7 is printed in bar (svp, pd:
freq_sorted[6,3] <- 0.06

freq_sorted$Tokens <- tokens_sorted$Freq  #add raw frequency to dataset
ice$ThemeComplexityConflated <- factor(ice$ThemeComplexityConflated, levels = c("s", "spp", "svp", "mpp", "mvp"))

(plot.thcox <- ggplot(ice, aes(ThemeComplexityConflated, fill=Resp)) + 
  geom_bar(aes(fill=Resp), position="fill",width = 0.6, color="black") + 
  labs(title = "", y="proportion of tokens", x="ThemeComplexity") + 
  theme(axis.text.x = element_text(size=16)) +
  theme(axis.title.y = element_text(size=16)) +
  theme(axis.title.x = element_text(size=16)) +
  geom_text(data = freq_sorted, aes(x=ThemeComplexityConflated,y=Freq, label = Tokens), vjust=1.5, color=ifelse(freq_sorted$Resp=="pd", "black", "white"), size=6) + 
  scale_fill_manual(values = cols, name="Variant", labels=c("ditransitive", "prepositional")) +
  geom_hline(yintercept = 0.5, colour="black", linetype=2) +
  guides(fill=guide_legend(ncol=2,title.position = "top", title.hjust = 0.5))+
    theme(legend.position = "bottom",legend.title = element_text(size=16,face = "bold"), legend.text = element_text(size=16)))

pdf("Rplot_cox_barplot_thcox.pdf", width = 6, height = 8)
plot.thcox
dev.off()

# -- .. all Complexity: together --------------
require(gridExtra)
grid.arrange(plot.reccox, plot.thcox, ncol=2)

my_legend <- g_legend(plot.reccox)


pdf("Rplot_cox_barplot_thandreccox.pdf", width = 12, height = 8)
grid.arrange(arrangeGrob(plot.reccox + theme(legend.position = "none"),
                         plot.thcox + theme(legend.position = "none") + ylab(""),
                         nrow=1), my_legend, nrow=2, heights=c(10,1))
dev.off()




# -- ..by variety and variant --------------------------------
(plot.reccox.var <- ggplot(ice, aes(RecComplexityConflated)) + 
  facet_wrap(~Variety) +
  geom_hline(yintercept = 0.5, colour="black", linetype=2, size=.3) +
  geom_bar(aes(fill=Resp), position="fill", width=0.6, color="black") +
  scale_fill_manual(values = cols, labels=c("ditransitive", "prepositional"), name="Variant") +
  labs(title = "", y="proportion of tokens", x="RecComplexity") + 
  guides(fill=guide_legend(ncol=2,title.position = "top", title.hjust = 0.5))+
  theme(legend.position = "bottom",legend.title = element_text(size=16,face = "bold"), legend.text = element_text(size=16)) +
  theme(axis.text.x = element_text(size=14)) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.title.x = element_text(size=14)) )
  

pdf("Rplot_cox_barplot_reccoxbyvar.pdf", width = 9, height = 6)
plot.reccox.var
dev.off()



(plot.thcox.var <- ggplot(ice, aes(ThemeComplexityConflated)) + 
  facet_wrap(~Variety) +
  geom_hline(yintercept = 0.5, colour="black", linetype=2, size=.3) +
  geom_bar(aes(fill=Resp), position="fill", width=0.6, color="black") +
  scale_fill_manual(values = cols, labels=c("ditransitive", "prepositional"), name="Variant") +
  labs(title = "", y="proportion of tokens", x="ThemeComplexity") + 
    guides(fill=guide_legend(ncol=2,title.position = "top", title.hjust = 0.5))+
    theme(legend.position = "bottom",legend.title = element_text(size=16,face = "bold"), legend.text = element_text(size=16)) +
  theme(axis.text.x = element_text(size=14)) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.title.x = element_text(size=14)) )

pdf("Rplot_cox_barplot_thcoxbyvar.pdf", width = 9, height = 6)
plot.thcox.var
dev.off()  

# -- CORRELATION BETWEEN LENGTH AND COMPLEXITY ---------------------------

# -- make complexity numerical -------------
ice <- cbind(ice, as.numeric2(ice[,86:87], as.char=F))
colnames(ice)[88:89] = c("NumRecComplexity", "NumThemeComplexity")
summary(ice$NumRecComplexity)
ice[ice$NumRecComplexity== 5, "RecComplexityConflated"]

ice$RecComplexity <- factor(ice$RecComplexity, levels=c("s", "postad", "gn", "pp", "ge", "appnom","co", "nonfin", "cp", "rc", "advc", "nc", "tpp", "tvp", "mpp", "mvp"))
ice$ThemeComplexity <- factor(ice$ThemeComplexity, levels=c("s", "postad", "gn", "pp", "ge", "appnom", "co", "nonfin", "cp", "rc", "advc", "nc", "tpp", "tvp", "mpp", "mvp"))

ice<- cbind(ice, as.numeric2(ice[,43:44], as.char=F))
colnames(ice)[90:91] = c("NumFullRecComplexity", "NumFullThemeComplexity")
summary(ice$NumFullRecComplexity)


# -- Correlations ------
cor.test(ice$NumRecComplexity, ice$RecLetterLth)
cor.test(ice$NumThemeComplexity, ice$ThemeLetterLth)

cor(ice$NumFullRecComplexity, ice$RecLetterLth)
cor(ice$NumFullThemeComplexity, ice$ThemeLetterLth)

cor.test(ice$NumRecComplexity, ice$RecWordLth)
cor.test(ice$ThemeWordLth, ice$NumThemeComplexity)



# -- Linear plots --------
ggplot(ice, aes(NumRecComplexity, RecLetterLth)) +
  geom_smooth(method='lm', colour="black", size=.5) +
  facet_wrap(~Variety)
  
ggplot(ice, aes(NumThemeComplexity, ThemeLetterLth)) +
  geom_smooth(method='lm', colour="black") +
  facet_wrap(~Variety)

corr.df.th <- as.data.frame(cbind(as.character(ice$Variety), ice$NumThemeComplexity, ice$ThemeLetterLth, rep("theme", nrow(ice))))
corr.df.rec <- as.data.frame(cbind(as.character(ice$Variety), ice$NumRecComplexity, ice$RecLetterLth, rep("recipient", nrow(ice))))
corr.df <- as.data.frame(rbind(corr.df.th, corr.df.rec))
names(corr.df) <- c("variety","complexity", "length", "constituent")
corr.df$complexity <- as.numeric2(corr.df$complexity)
corr.df$length <- as.numeric2(corr.df$length)

# Both in one
shapes = c(theme=1, recipient=2)
(plot.corr.lm <- ggplot(corr.df, aes(complexity, length)) +
    geom_point(aes(shape=constituent), size=1.5, position = position_dodge(width=0.5)) +
    geom_smooth(aes(linetype=constituent), method="lm", size=.5, color="black") +
    geom_smooth(aes(linetype=constituent), method="lm", size=.5, color="black", fill=NA) +
    facet_wrap(~variety) +
    theme(axis.title = element_text(size=14)) +
    scale_shape_manual(values = shapes, name="Constituent") +
    scale_linetype_manual(values=c("solid", "dashed"), name="") +
    theme(legend.position = "bottom", legend.title=element_text(size=13, face="bold")) +
    guides(linetype=guide_legend(override.aes=list(fill=NA), ncol=1, byrow=T, title.position = "top", title.hjust = 0.5)) +
    guides(shape=guide_legend(ncol=1, title.position="top", byrow=T, title.hjust = 0.5)) +
    theme(legend.text=element_text(size=13)) + 
    labs(y="length in letters", x="level of complexity") )


pdf("Rplot_cox_smooth_correlation.pdf", width = 8, height = 10)
plot.corr.lm
dev.off()  

# -- long but simple recipients ----------------
summary(ice[ice$RecComplexityConflated=="s", "RecLetterLth"])
quantile(ice[ice$RecComplexityConflated=="s", "RecLetterLth"], probs=seq(0,1,0.005))
ice[ice$RecComplexityConflated=="s" & ice$RecLetterLth > 33, "RecipientShort"]
ice[ice$RecComplexityConflated=="s" & ice$RecLetterLth > 33, "FileID"]

ice[ice$ThemeComplexityConflated=="s" & ice$ThemeLetterLth > 33, "FileID"]
ice[ice$ThemeComplexityConflated=="s" & ice$ThemeLetterLth > 33, "Resp"]
quantile(ice[ice$ThemeComplexityConflated=="s", "ThemeLetterLth"], probs=seq(0,1,0.005))

# -- ..plot simple constituents ---------------
cox.df <- data[data$RecBinComplexity=="simple",c("RecLetterLth", "Resp")]
cox.df$Constituent <- as.factor(rep("recipient", nrow(cox.df)))
names(cox.df) <- c("Length", "Resp", "Constituent")

cox.df.th <- data[data$ThemeBinComplexity=="simple",c("ThemeLetterLth", "Resp")]
cox.df.th$Constituent <- as.factor(rep("theme", nrow(cox.df.th)))
names(cox.df.th) <- c("Length", "Resp", "Constituent")

cox.df.all <- rbind(cox.df, cox.df.th)


mean.all <- c(mean(data[data$Resp=="do", "RecHeadFreq"]), mean(data[data$Resp=="do", "ThemeHeadFreq"]), mean(data[data$Resp=="pd", "RecHeadFreq"]), mean(data[data$Resp=="pd", "ThemeHeadFreq"]))
mean.all <- as.data.frame(mean.all)
mean.all$Resp <- as.factor(c("do", "do", "pd", "pd"))
mean.all$Constituent <- as.factor(c("recipient", "theme", "recipient", "theme"))
mean.all$mean.all <- round(mean.all$mean.all, 1)

cols3 <- c(recipient="dimgray", theme="gainsboro")

(plot.simple <- ggplot(cox.df.all, aes(Constituent, Length)) +
    geom_boxplot(position=position_dodge(0.8), width=0.7, outlier.alpha = 0.3, outlier.size = 0.8) +
    labs(title = "", y="length in letters", x="") + 
    theme(legend.position="bottom", legend.title=element_text(size=16, face="bold")) + 
    theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16)) +
    theme(legend.text=element_text(size=16)) +
    theme(axis.title = element_text(size=14)) +
    guides(fill=guide_legend(ncol=2,title.position = "top", title.hjust = 0.5)))

pdf("Rplot_simple_constituents_length.pdf", width = 6, height = 5)
plot.simple
dev.off()



# -- Random forest -------------------------------------------------------
#detach("package:partykit", unload=TRUE)
# ---- ..generate complexity ratio ------------
ice$ComplexityRatio <- ice$NumRecComplexity/ice$NumThemeComplexity
summary(ice$ComplexityRatio)

# calculate correlation between ComplexityRatio and WeightRatio
cor.test(ice$WeightRatio, ice$ComplexityRatio) # rather low

# -- ..fit the forest ------------
library(party)
forest.controls = cforest_unbiased(ntree=5000, mtry=5)
set.seed(3452341)


rf1a <- cforest(Resp ~
                  Variety +  
                  Mode +
                  VerbSemantics +
                  z.logWeightRatio +
                  RecGivenness +
                  ThemeGivenness +
                  RecDefiniteness +
                  ThemeDefiniteness +
                  ComplexityRatio +
                  z.RecHeadFreq +
                  z.ThemeHeadFreq +
                  z.RecThematicity +
                  z.ThemeThematicity +
                  PrimeType +
                  RecPron +
                  ThemePron +
                  RecBinAnimacy +
                  ThemeBinAnimacy +
                  z.TypeTokenRatio,
                data = ice, controls = forest.controls)

varimp.rf1a = party::varimpAUC(rf1a)


# -- ..verify forest with different seed ------------------
forest.controls = cforest_unbiased(ntree=5000, mtry=5)
set.seed(67852)


rf1b <- cforest(Resp ~
                  Variety +  
                  Mode +
                  VerbSemantics +
                  z.logWeightRatio +
                  RecGivenness +
                  ThemeGivenness +
                  RecDefiniteness +
                  ThemeDefiniteness +
                  ComplexityRatio +
                  z.RecHeadFreq +
                  z.ThemeHeadFreq +
                  z.RecThematicity +
                  z.ThemeThematicity +
                  PrimeType +
                  RecPron +
                  ThemePron +
                  RecBinAnimacy +
                  ThemeBinAnimacy +
                  z.TypeTokenRatio,
                data = ice, controls = forest.controls)

varimp.rf1b = party::varimpAUC(rf1b)
dot1b <- dotplot(sort(varimp.rf1b))
dot1a <- dotplot(sort(varimp.rf1a))

require(gridExtra)
grid.arrange(dot1a, dot1b, ncol=2)


# -- ..evaluation of the random forest --------------------
# predictive accuracy of forest
probs.ratio <- treeresponse(rf1b)
# Add predictions for the prepositional dative variant to dataframe
ice$probs.rf1b <- sapply(probs.ratio, FUN = function(x) return(x[1])) %>%
  as.vector
ice$preds.rf1b <- predict(rf1b)
(acc.rf1b <- sum(ice$Resp==ice$preds.rf1b)/nrow(ice)) # 90.1%


somers2(1-ice$probs.rf1b, as.numeric(ice$Resp)-1)


# --..Plotting of varimp-----------------------
varimpAUC.rf1b <- sort(varimp.rf1b)
varimp.all <- data.frame(vals = varimpAUC.rf1b, 
                         preds = names(varimpAUC.rf1b),
                         levels = names(sort(varimpAUC.rf1b)),
                         var = rep("Variable Importance", 19))

# change names of predictors:
varimp.all$preds <- c("TTR", "ThemeAnimacy", "RecDefiniteness", "RecThematicity", "ThemeThematicity", "PrimeType", "ThemeGivenness", "RecAnimacy", "Mode", "VerbSemantics", "Variety", "RecHeadFreq", "RecGivenness", "ThemeHeadFreq", "ThemeDefiniteness", "ThemePron", "ComplexityRatio", "RecPron", "WeightRatio")
varimp.all$levels <- c("TTR", "ThemeAnimacy", "RecDefiniteness", "RecThematicity", "ThemeThematicity", "PrimeType", "ThemeGivenness", "RecAnimacy", "Mode", "VerbSemantics", "Variety", "RecHeadFreq", "RecGivenness", "ThemeHeadFreq", "ThemeDefiniteness", "ThemePron", "ComplexityRatio", "RecPron", "WeightRatio")

varimp.all$preds <- factor(varimp.all$preds, as.character(varimp.all$preds))

(plot.rf.varimp <- ggplot(varimp.all, aes(preds,vals)) +
  geom_bar(fill="dimgray", stat="identity", color="black", width=.6) +
  coord_flip() +
  theme(axis.text = element_text(size=15, colour = "black")) +
  theme(title = element_text(size=18)) +
  labs(x="", y="Conditional variable importance") +
  # add the Shih line
  geom_hline(yintercept=abs(min(varimp.all$vals)),linetype="dashed") )

pdf("Rplot_cox_crf_varimp_all.pdf", width = 7, height=7)
plot.rf.varimp
dev.off()






# -- REGIONAL VARIABILITY OF COMPLEXITY --------------------------
# no glmer model with separate Complexity due to low levels

# -- Random forest by variety -----------
# (only with ratio data)

vars <- levels(ice$Variety)

# formula
f1 <- Resp ~
  Mode +
  VerbSemantics +
  z.logWeightRatio +
  RecGivenness +
  ThemeGivenness +
  RecDefiniteness +
  ThemeDefiniteness +
  ComplexityRatio +
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

forest.controls = cforest_unbiased(ntree=5000, mtry=5)

t1 <- proc.time()
for (i in 1:9){
  cat("Getting", vars[i], "...\n")
  d <- droplevels(subset(ice, ice$Variety == vars[i]))
  
  # rescale predictors
  d$z.logWeightRatio = z.(d[,"logWeightRatio"], factor=2)
  d$z.RecHeadFreq = z.(d[,"RecHeadFreq"], factor=2)
  d$z.ThemeHeadFreq = z.(d[,"ThemeHeadFreq"], factor=2)
  d$z.RecThematicity = z.(d[,"RecThematicity"], factor=2)
  d$z.ThemeThematicity = z.(d[,"ThemeThematicity"], factor=2)
  d$z.TypeTokenRatio = z.(d[,"TypeTokenRatio"], factor=2)
  
  # run model
  rf <- cforest(f1, data = d, controls = forest.controls)
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
rownames(pars.sorted) <- c("ThemeAnimacy", "TTR", "PrimeType", "RecDefiniteness", "RecThematicity", "Mode", "ThemeGivenness", "VerbSemantics","ThemeThematicity", "RecAnimacy", "ThemeDefiniteness", "ThemeHeadFreq" , "ThemePron", "RecGivenness", "RecHeadFreq" ,  "ComplexityRatio", "RecPron", "WeightRatio")

pars_long <- melt(t(pars.sorted))
names(pars_long) <- c("Variety", "Factor", "Varimp")

pars_long$Variety <- factor(pars_long$Variety, levels=c("CAN", "GB", "HK", "IND", "IRE", "JA", "NZ", "PHI", "SIN"))


(plot.rf.varimp.VAR <- ggplot(pars_long, aes(Factor,Varimp)) +
    geom_bar(fill="dimgray", stat="identity", color="black", width=.6) +
    coord_flip() +
    theme(axis.text = element_text(size=10, colour = "black")) +
    theme(title = element_text(size=12)) +
    labs(x="", y="Conditional variable importance") +
    facet_wrap(~Variety) +
    # add the Shih line
    geom_hline(yintercept=abs(min(varimp.all$vals)),linetype="dashed") )

pdf("Rplot_cox_crf_varimp_byvar.pdf", width = 9, height=10)
plot.rf.varimp.VAR
dev.off()



# -- GLMER models by variety ------------
contrasts(ice$Variety.Sum)

m01 <- glmer(Resp ~ (1|GenreCoarse) + (1|SpeakerID) + (1|Verb) +
               Variety.Sum * z.logWeightRatio   
                , data=ice, family=binomial, 
             glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1e6)))

summary(m01)
somers.mer(m01) 
m01.fits <- fitted(m01)
sum(ifelse(m01.fits > 0.5, 1, 0) == as.numeric(ice$Resp)-1)/nrow(ice) 
MuMIn::r.squaredGLMM(m01) 

ggLogit.plot(m01, ice)
collin.fnc.mer(m01)$cnumber
overdisp.mer(m01)
m01.resid.bin <- as.data.frame(binned.resids(m01.fits, (as.numeric(ice$Resp)-1)-m01.fits, nclass=40)$binned)
ggplot(m01.resid.bin, aes(xbar, ybar)) +
  geom_ribbon(aes(ymin=ybar - 2*se, ymax = ybar + 2*se), alpha=0.1) +
  geom_point() + labs(x="Estimated Pr(pd-dative)", y="mean residual") +
  geom_smooth(method="loess", se=F)


# Check out SinE
ice$Variety.Sum2 <- factor(ice$Variety, levels=c("SIN","GB","HK","IND","IRE", "JA","NZ","PHI", "CAN"))
contrasts(ice$Variety.Sum2) <- contr.sum(9)

m01.sg <- glmer(Resp ~ (1|GenreCoarse) + (1|SpeakerID) + (1|Verb) +
               Variety.Sum2 * z.logWeightRatio   
             , data=ice, family=binomial, 
             glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1e6)))
summary(m01.sg)
summary(m01)


m02 <- glmer(Resp ~ (1 | GenreCoarse) + (1|SpeakerID) + (1|Verb) + 
               Variety.Sum * ComplexityRatio   
             , data=ice, family=binomial, 
             glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1e6)))

summary(m02)
somers.mer(m02)  
m02.fits <- fitted(m02)
sum(ifelse(m02.fits > 0.5, 1, 0) == as.numeric(ice$Resp)-1)/nrow(ice) 
srd.m02 <- MuMIn::r.squaredGLMM(m02)

ggLogit.plot(m02, ice)
collin.fnc.mer(m02)$cnumber
overdisp.mer(m02) # yes, overdispersion, i.e. there is more variation in the data than predicted by the model
m02.resid.bin <- as.data.frame(binned.resids(m02.fits, (as.numeric(ice$Resp)-1)-m02.fits, nclass=40)$binned)
ggplot(m02.resid.bin, aes(xbar, ybar)) +
  geom_ribbon(aes(ymin=ybar - 2*se, ymax = ybar + 2*se), alpha=0.1) +
  geom_point() + labs(x="Estimated Pr(pd-dative)", y="mean residual") +
  geom_smooth(method="loess", se=F)


# check out SinE
m02.sg <- glmer(Resp ~ (1 | GenreCoarse) + (1|SpeakerID) + (1|Verb) + 
               Variety.Sum2 * ComplexityRatio   
             , data=ice, family=binomial, 
             glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1e6)))
summary(m02.sg)


# -- NON-PRONOMINAL MODELS --------------
# -- forest without pronouns -----------------------------------
nouns <- droplevels(subset(ice, ice$RecPron=="non-pron" & ice$ThemePron=="non-pron"))

forest.controls = cforest_unbiased(ntree=3000, mtry=5)
set.seed(3452341)


rf2a <- cforest(Resp ~
                  Variety +  
                  Mode +
                  VerbSemantics +
                  z.logWeightRatio +
                  RecGivenness +
                  ThemeGivenness +
                  RecDefiniteness +
                  ThemeDefiniteness +
                  ComplexityRatio +
                  z.RecHeadFreq +
                  z.ThemeHeadFreq +
                  z.RecThematicity +
                  z.ThemeThematicity +
                  PrimeTypePruned +
                  #RecPron +
                  #ThemePron +
                  RecBinAnimacy +
                  ThemeBinAnimacy +
                  z.TypeTokenRatio,
                data = nouns, controls = forest.controls)

varimp.rf2a = party::varimpAUC(rf2a)


# -- ..verify forest with different seed ------------------

set.seed(67852)


rf2b <- cforest(Resp ~
                  Variety +  
                  Mode +
                  VerbSemantics +
                  z.logWeightRatio +
                  RecGivenness +
                  ThemeGivenness +
                  RecDefiniteness +
                  ThemeDefiniteness +
                  ComplexityRatio +
                  z.RecHeadFreq +
                  z.ThemeHeadFreq +
                  z.RecThematicity +
                  z.ThemeThematicity +
                  PrimeTypePruned +
                  # RecPron +
                  #ThemePron +
                  RecBinAnimacy +
                  ThemeBinAnimacy +
                  z.TypeTokenRatio,
                data = nouns, controls = forest.controls)

varimp.rf2b = party::varimpAUC(rf2b)
dot2b <- dotplot(sort(varimp.rf2b))
dot2a <- dotplot(sort(varimp.rf2a))

require(gridExtra)
grid.arrange(dot2a, dot2b, ncol=2)


varimp.rf2a-varimp.rf2b # hardly any difference



# --..Plotting of varimp without pronouns-----------------------
varimpAUC.rf2b <- sort(varimp.rf2b)
varimp.all <- data.frame(vals = varimpAUC.rf2b, 
                         preds = names(varimpAUC.rf2b),
                         levels = names(sort(varimpAUC.rf2b)),
                         var = rep("Variable Importance", 17))



# change names of predictors:
varimp.all$preds <- c("ThemeAnimacy", "ThemeGivenness", "Mode","RecHeadFreq", "RecGivenness", "ThemeThematicity", "RecDefiniteness", "TTR", "PrimeType", "RecThematicity", "ThemeDefiniteness", "ThemeHeadFreq", "RecAnimacy", "VerbSemantics", "Variety", "ComplexityRatio", "WeightRatio")
varimp.all$levels <- c("ThemeAnimacy", "ThemeGivenness", "Mode","RecHeadFreq", "RecGivenness", "ThemeThematicity", "RecDefiniteness", "TTR", "PrimeType", "RecThematicity", "ThemeDefiniteness", "ThemeHeadFreq", "RecAnimacy", "VerbSemantics", "Variety", "ComplexityRatio", "WeightRatio")

varimp.all$preds <- factor(varimp.all$preds, as.character(varimp.all$preds))

(plot.rf.varimp.wo.pron <- ggplot(varimp.all, aes(preds,vals)) +
    geom_bar(fill="dimgray", stat="identity", color="black", width=.6) +
    coord_flip() +
    theme(axis.text = element_text(size=15, colour = "black")) +
    theme(title = element_text(size=18)) +
    labs(x="", y="Conditional variable importance") +
    # add the Shih line
    geom_hline(yintercept=abs(min(varimp.all$vals)),linetype="dashed") )

pdf("Rplot_cox_crf_varimp_wo_pron.pdf", width = 7, height=7)
plot.rf.varimp.wo.pron
dev.off()



# -- glmer model of interactions without pronouns ----------------
contrasts(nouns$Variety.Sum) <- contr.sum(9)

m.nom.length <- glmer(Resp ~ (1 | GenreCoarse) + (1|SpeakerID) + (1|Verb) +
                        WeightRatio * Variety.Sum
                      , data=nouns, family=binomial, 
                      glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1e6)))



m.nom.cox <- glmer(Resp ~ (1 | GenreCoarse) + (1|SpeakerID) + (1|Verb) +
                     ComplexityRatio * Variety.Sum
                   , data=nouns, family=binomial, 
                   glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1e6)))


somers.mer(m.nom.length) # C: 0.9095071
somers.mer(m.nom.cox) # C: 0.9157313
m.nom.length.fits <- fitted(m.nom.length)
m.nom.cox.fits <- fitted(m.nom.cox)
sum(ifelse(m.nom.length.fits > 0.5, 1, 0) == as.numeric(nouns$Resp)-1)/nrow(nouns) # 
sum(ifelse(m.nom.cox.fits > 0.5, 1, 0) == as.numeric(nouns$Resp)-1)/nrow(nouns) #

summary(m.nom.length)
summary(m.nom.cox)


# -- IS A 5-level PREDICTOR REALLY NECESSARY -------------------------------
# compare model with binary to that with 5-level predictor (excluding or including length) and see which one scores better

# --  random forest -------------------------
forest.controls = cforest_unbiased(ntree=5000, mtry=5)
set.seed(3452341)

# -- ..Random forest with two levels----------
rf3a <- cforest(Resp ~
                  Variety +  
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
                data = ice, controls = forest.controls)

varimp.rf3a = party::varimpAUC(rf3a)

# -- ..verify crf with different seed -------------------

forest.controls = cforest_unbiased(ntree=5000, mtry=5)
set.seed(896713)

rf3b <- cforest(Resp ~    Variety +  
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
                        data = ice, controls = forest.controls)

varimp.rf3b = party::varimpAUC(rf3b)


# -- ..Random forest with 5 levels and separate constituents -------------
set.seed(2367)
rf4a <- cforest(Resp ~
                  Variety +  
                  Mode +
                  VerbSemantics +
                  z.logWeightRatio +
                  RecGivenness +
                  ThemeGivenness +
                  RecDefiniteness +
                  ThemeDefiniteness +
                  RecComplexityConflated +
                  ThemeComplexityConflated +
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
                data = ice, controls = forest.controls)
varimp.rf4a = party::varimpAUC(rf4a)

# -- ..verify with different seed -------------
set.seed(457994)
rf4b <- cforest(Resp ~
                  Variety +  
                  Mode +
                  VerbSemantics +
                  z.logWeightRatio +
                  RecGivenness +
                  ThemeGivenness +
                  RecDefiniteness +
                  ThemeDefiniteness +
                  RecComplexityConflated +
                  ThemeComplexityConflated +
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
                data = ice, controls = forest.controls)
varimp.rf4b = party::varimpAUC(rf4b)

# -- .. check predictive accuracy--------
# predictive accuracy of forest (is it better than rf1a?)
probs.bin <- treeresponse(rf3b)
# Add predictions for the prepositional dative variant to dataframe
ice$probs.rf3b <- sapply(probs.bin, FUN = function(x) return(x[1])) %>%
  as.vector
ice$preds.rf3b <- predict(rf3b)
(acc.rf3b <- sum(ice$Resp==ice$preds.rf3b)/nrow(ice)) # %


somers2(1-ice$probs.rf3b, as.numeric(ice$Resp)-1)


# -- glmer -------------------

m02b <- glmer(Resp ~ (1 | GenreCoarse) + (1|SpeakerID) + (1|Verb) +
               RecBinComplexity + ThemeBinComplexity
             , data=ice, family=binomial, 
             glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1e6)))

m02c <- glmer(Resp ~ (1 | GenreCoarse) + (1|SpeakerID) + (1|Verb) +
                RecComplexityConflated + ThemeComplexityConflated
              , data=ice, family=binomial, 
              glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1e6)))

summary(m02b)
summary(m02c)
somers.mer(m02b) # C: 0.9296211
somers.mer(m02c) # C: 0.9316025
m02b.fits <- fitted(m02b)
m02c.fits <- fitted(m02c)
sum(ifelse(m02b.fits > 0.5, 1, 0) == as.numeric(ice$Resp)-1)/nrow(ice) # 0.8563701
sum(ifelse(m02c.fits > 0.5, 1, 0) == as.numeric(ice$Resp)-1)/nrow(ice) # 0.8572533
srd.m02b <- MuMIn::r.squaredGLMM(m02b)
srd.m02c <- MuMIn::r.squaredGLMM(m02c)

anova(m02, m02b)
AIC(m02b)
AIC(m02c)

ggLogit.plot(m02b, ice)
collin.fnc.mer(m02b)$cnumber
overdisp.mer(m02b) # no overdisperson
overdisp.mer(m02c) # no overdisperson
m02b.resid.bin <- as.data.frame(binned.resids(m02b.fits, (as.numeric(ice$Resp)-1)-m02b.fits, nclass=40)$binned)
ggplot(m02b.resid.bin, aes(xbar, ybar)) +
  geom_ribbon(aes(ymin=ybar - 2*se, ymax = ybar + 2*se), alpha=0.1) +
  geom_point() + labs(x="Estimated Pr(pd-dative)", y="mean residual") +
  geom_smooth(method="loess", se=F)

save.image("datives_complexity.RData")