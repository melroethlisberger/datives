#=======================================================================
#   Script: 00b_datives_descriptive_stats.R
#   Date last changed: 5 January 2018
#
#   Rversion: 3.3.2
#   Author: MR
#   This script creates the plots and calculates the X-squares for
#   the descriptive statistics of the dissertation.
#=======================================================================

#----SETUP: LIBRARIES-----------------------------------------
library(reshape2);library(plyr);library(dplyr);library(magrittr);library(lme4);library(ggplot2);library(gridExtra); library(JGmisc); library(beepr); library(effects); library(Hmisc)

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



#-----SETUP: DATA------------------------
data <- read.delim("OA/datives.txt")


#-- PLOT: Proportional distributions ---------------------------------
cols <- c(do="dimgray",pd="gainsboro")
data$Variety <- factor(data$Variety, levels=c("GB", "CAN", "IRE", "NZ", "JA", "SIN", "IND", "HK", "PHI"))

frequencies = prop.table(with(data, table(Variety, Resp)),1) # create prop table
frequencies <- data.frame(frequencies)
tokens = with(data, table(Variety, Resp)) # create raw freq table
tokens <- data.frame(tokens)
tokens_sorted <- arrange(tokens, Variety, Resp) # sort both using library(plyr)
freq_sorted <- arrange(frequencies, Variety, Resp)

# replace every second percentage by 1 to use as geom_text
freq_sorted[1,3] <- 1
freq_sorted[3,3] <- 1
freq_sorted[5,3] <- 1
freq_sorted[7,3] <- 1
freq_sorted[9,3] <- 1
freq_sorted[11,3] <- 1
freq_sorted[13,3] <- 1
freq_sorted[15,3] <- 1
freq_sorted[17,3] <- 1

freq_sorted$Tokens <- tokens_sorted$Freq  #add raw frequency to dataset


# --> add raw numbers into cells:
p1 = ggplot(data, aes(Variety, fill=Resp)) + 
  geom_bar(aes(fill=Resp), position="fill",width = 0.7, color="black") + 
  labs(title = "", y="proportion of tokens", x="") + 
  theme(legend.position="bottom", legend.title=element_text(size=16, face="bold")) + 
  theme(axis.text.x = element_text(size=16)) +
  theme(legend.text=element_text(size=16)) +
  geom_text(data = freq_sorted, aes(x=Variety,y=Freq, label = Tokens), vjust=2, color=ifelse(freq_sorted$Resp=="pd", "black", "white"), size=8) + 
  scale_fill_manual(values = cols, labels=c("ditransitive", "prepositional"), name="Variant") +
  guides(fill=guide_legend(ncol=2,title.position = "top", title.hjust = 0.5));p1

pdf("Rplot_distribution_all_datives.pdf", width = 11, height = 7)
p1
dev.off()

# Chisquare test on proportions
tab.vars <- with(data, table(Variety, Resp))
prop.test(tab.vars, correct=FALSE) 


#-- PLOT: Proportional distributions: glowbe only -----------------------------
glowbe <- droplevels(subset(data, data$Corpus=="glowbe"))

frequencies = prop.table(with(glowbe, table(Variety, Resp)),1) # create prop table
frequencies <- data.frame(frequencies)
tokens = with(glowbe, table(Variety, Resp)) # create raw freq table
tokens <- data.frame(tokens)
tokens_sorted <- arrange(tokens, Variety, Resp) # sort both using library(plyr)
freq_sorted <- arrange(frequencies, Variety, Resp)

# replace every second percentage by 1 to use as geom_text
freq_sorted[1,3] <- 1
freq_sorted[3,3] <- 1
freq_sorted[5,3] <- 1
freq_sorted[7,3] <- 1
freq_sorted[9,3] <- 1
freq_sorted[11,3] <- 1
freq_sorted[13,3] <- 1
freq_sorted[15,3] <- 1
freq_sorted[17,3] <- 1

freq_sorted$Tokens <- tokens_sorted$Freq  #add raw frequency to dataset


# --> add raw numbers into cells:
p2 = ggplot(glowbe, aes(Variety, fill=Resp)) + 
  geom_bar(aes(fill=Resp), position="fill",width = 0.7, color="black") + 
  labs(title = "", y="proportion of tokens", x="GloWbE") + 
  theme(legend.position="bottom", legend.title=element_text(size=16, face="bold")) + 
  theme(axis.text.x = element_text(size=16)) +
  theme(axis.title = element_text(size=14)) +
  theme(legend.text=element_text(size=16)) +
  geom_text(data = freq_sorted, aes(x=Variety,y=Freq, label = Tokens), vjust=2, color=ifelse(freq_sorted$Resp=="pd", "black", "white"), size=6.5) + 
  scale_fill_manual(values = cols, labels=c("ditransitive", "prepositional"), name="Variant") +
  guides(fill=guide_legend(ncol=2,title.position = "top", title.hjust = 0.5));p2

#-- PLOT: Proportional distributions: ice only -----------------------------
ice <- droplevels(subset(data, data$Corpus=="ice"))

frequencies = prop.table(with(ice, table(Variety, Resp)),1) # create prop table
frequencies <- data.frame(frequencies)
tokens = with(ice, table(Variety, Resp)) # create raw freq table
tokens <- data.frame(tokens)
tokens_sorted <- arrange(tokens, Variety, Resp) # sort both using library(plyr)
freq_sorted <- arrange(frequencies, Variety, Resp)

# replace every second percentage by 1 to use as geom_text
freq_sorted[1,3] <- 1
freq_sorted[3,3] <- 1
freq_sorted[5,3] <- 1
freq_sorted[7,3] <- 1
freq_sorted[9,3] <- 1
freq_sorted[11,3] <- 1
freq_sorted[13,3] <- 1
freq_sorted[15,3] <- 1
freq_sorted[17,3] <- 1

freq_sorted$Tokens <- tokens_sorted$Freq  #add raw frequency to dataset


# --> add raw numbers into cells:
p4 = ggplot(ice, aes(Variety, fill=Resp)) + 
  geom_bar(aes(fill=Resp), position="fill",width = 0.7, color="black") + 
  labs(title = "", y="proportion of tokens", x="ICE") + 
  theme(legend.position="bottom", legend.title=element_text(size=20, face="bold")) + 
  theme(axis.text.x = element_text(size=16)) +
  theme(axis.title = element_text(size=14)) +
  theme(legend.text=element_text(size=20)) +
  geom_text(data = freq_sorted, aes(x=Variety,y=Freq, label = Tokens), vjust=2, color=ifelse(freq_sorted$Resp=="pd", "black", "white"), size=6.5) + 
  scale_fill_manual(values = cols, labels=c("ditransitive", "prepositional"), name="Variant") +
  guides(fill=guide_legend(ncol=2,title.position = "top", title.hjust = 0.5));p4

require(gridExtra)
mylegend <- g_legend(p4)
bothplots.distr <- grid.arrange(arrangeGrob(p2 + theme(legend.position = "none"),
                                            p4 + theme(legend.position = "none") + ylab(""),
                                            nrow=1), mylegend, nrow=2, heights=c(10,1))

pdf("Rplot_distribution_iceVSglowbe.pdf", width = 15, height = 7)
grid.arrange(arrangeGrob(p4 + theme(legend.position = "none"),
                         p2 + theme(legend.position = "none") + ylab(""),
                         nrow=1), mylegend, nrow=2, heights=c(10,1))
dev.off()

# Chisquare test on proportions
tab.vars.ice <- with(data[data$Corpus=="ice",], table(Variety, Resp))
prop.test(tab.vars.ice, correct=FALSE) 

tab.vars.glowbe <- with(data[data$Corpus=="glowbe",], table(Variety, Resp))
prop.test(tab.vars.glowbe, correct=FALSE) 


#-- PLOT: proportional distributions by Genre -------------------------

frequencies = prop.table(with(data, table(GenreCoarse, Resp)),1) # create prop table
frequencies <- data.frame(frequencies)
tokens = with(data, table(GenreCoarse, Resp)) # create raw freq table
tokens <- data.frame(tokens)
tokens_sorted <- arrange(tokens, GenreCoarse, Resp) # sort both using library(plyr)
freq_sorted <- arrange(frequencies, GenreCoarse, Resp)

# replace every second percentage by 1 to use as geom_text
freq_sorted[1,3] <- 1
freq_sorted[3,3] <- 1
freq_sorted[5,3] <- 1
freq_sorted[7,3] <- 1
freq_sorted[9,3] <- 1


freq_sorted$Tokens <- tokens_sorted$Freq  #add raw frequency to dataset


# --> add raw numbers into cells:
plot.genre = ggplot(data, aes(GenreCoarse, fill=Resp)) + 
  geom_bar(aes(fill=Resp), position="fill",width = 0.7, color="black") + 
  labs(title = "", y="proportion of tokens", x="") + 
  theme(legend.position="bottom", legend.title=element_text(size=16, face="bold")) + 
  theme(axis.text.x = element_text(size=16)) +
  theme(legend.text=element_text(size=16)) +
  geom_text(data = freq_sorted, aes(x=GenreCoarse,y=Freq, label = Tokens), vjust=2, color=ifelse(freq_sorted$Resp=="pd", "black", "white"), size=6.5) + 
  scale_fill_manual(values = cols, labels=c("ditransitive", "prepositional"), name="Variant") +
  guides(fill=guide_legend(ncol=2,title.position = "top", title.hjust = 0.5));plot.genre


pdf("Rplot_distribution_GenreCoarse.pdf", width = 8, height = 5)
plot.genre
dev.off()

# Chisquare test
tab.reg <- with(data, table(GenreCoarse, Resp))
prop.test(tab.reg, correct=FALSE) 




#-- PLOT: proportional distributions by Animacy -------------------------
# Recipient animacy
frequencies = prop.table(with(data, table(RecBinAnimacy, Resp)),1) # create prop table
frequencies <- data.frame(frequencies)
tokens = with(data, table(RecBinAnimacy, Resp)) # create raw freq table
tokens <- data.frame(tokens)
tokens_sorted <- arrange(tokens, RecBinAnimacy, Resp) # sort both using library(plyr)
freq_sorted <- arrange(frequencies, RecBinAnimacy, Resp)

# replace every second percentage by 1 to use as geom_text
freq_sorted[1,3] <- 1
freq_sorted[3,3] <- 1



freq_sorted$Tokens <- tokens_sorted$Freq  #add raw frequency to dataset


# --> add raw numbers into cells:
plot.recani = ggplot(data, aes(RecBinAnimacy, fill=Resp)) + 
  geom_bar(aes(fill=Resp), position="fill",width = 0.7, color="black") + 
  labs(title = "", y="proportion of tokens", x="RecAnimacy") + 
  theme(legend.position="bottom", legend.title=element_text(size=16, face="bold")) + 
  theme(axis.text.x = element_text(size=16)) +
  theme(legend.text=element_text(size=16)) +
  theme(axis.title = element_text(size=14)) +
  geom_text(data = freq_sorted, aes(x=RecBinAnimacy,y=Freq, label = Tokens), vjust=2, color=ifelse(freq_sorted$Resp=="pd", "black", "white"), size=6.5) + 
  scale_fill_manual(values = cols, labels=c("ditransitive", "prepositional"), name="Variant") +
  guides(fill=guide_legend(ncol=2,title.position = "top", title.hjust = 0.5));plot.recani


pdf("Rplot_distribution_RecAnimacy.pdf", width = 4, height = 5)
plot.recani
dev.off()


# Chisquare test
tab.recani <- with(data, table(RecBinAnimacy, Resp))
prop.test(t(tab.recani), correct=FALSE) 

# Theme animacy
# reverse levels
data$ThemeBinAnimacy <- factor(data$ThemeBinAnimacy, levels=c("animate", "inanimate"))
frequencies = prop.table(with(data, table(ThemeBinAnimacy, Resp)),1) # create prop table
frequencies <- data.frame(frequencies)
tokens = with(data, table(ThemeBinAnimacy, Resp)) # create raw freq table
tokens <- data.frame(tokens)
tokens_sorted <- arrange(tokens, ThemeBinAnimacy, Resp) # sort both using library(plyr)
freq_sorted <- arrange(frequencies, ThemeBinAnimacy, Resp)

# replace every second percentage by 1 to use as geom_text
freq_sorted[1,3] <- 1
freq_sorted[3,3] <- 1



freq_sorted$Tokens <- tokens_sorted$Freq  #add raw frequency to dataset


# --> add raw numbers into cells:
plot.thani = ggplot(data, aes(ThemeBinAnimacy, fill=Resp)) + 
  geom_bar(aes(fill=Resp), position="fill",width = 0.7, color="black") + 
  labs(title = "", y="proportion of tokens", x="ThemeAnimacy") + 
  theme(legend.position="bottom", legend.title=element_text(size=16, face="bold")) + 
  theme(axis.text.x = element_text(size=16)) +
  theme(legend.text=element_text(size=16)) +
  theme(axis.title = element_text(size=14)) +
  geom_text(data = freq_sorted, aes(x=ThemeBinAnimacy,y=Freq, label = Tokens), vjust=2, color=ifelse(freq_sorted$Resp=="pd", "black", "white"), size=6.5) + 
  scale_fill_manual(values = cols, labels=c("ditransitive", "prepositional"), name="Variant") +
  guides(fill=guide_legend(ncol=2,title.position = "top", title.hjust = 0.5));plot.thani


pdf("Rplot_distribution_ThemeAnimacy.pdf", width = 4, height = 5)
plot.thani
dev.off()

# Chisquare test
tab.theani <- with(data, table(ThemeBinAnimacy, Resp))
prop.test(t(tab.theani), correct=FALSE) 


require(gridExtra)
mylegend <- g_legend(plot.thani)


pdf("Rplot_distribution_animacy.pdf", width = 8, height = 6)
grid.arrange(arrangeGrob(plot.recani + theme(legend.position = "none"),
                         plot.thani + theme(legend.position = "none") + ylab(""), nrow=1), mylegend, nrow=2, heights=c(10,1))
dev.off()

# reverse ThemeAnimacy again
data$ThemeBinAnimacy <- factor(data$ThemeBinAnimacy, levels=c("inanimate", "animate"))



#-- PLOT: proportional distributions by Length -------------------------

library(reshape2)
mydata.summary = melt(xtabs(~Resp+RecLetterLth,data=data))
mydata.summary2 = mydata.summary[mydata.summary$Resp=="do",]
mydata.summary2$percent <- mydata.summary2$value/melt(xtabs(~RecLetterLth,data=data))$value
mydata.summary2$constituent <- rep(as.factor("recipient"), nrow(mydata.summary2))
names(mydata.summary2) <- c("Resp", "LetterLth", "value", "percent", "constituent")

mydata.summary.th = melt(xtabs(~Resp+ThemeLetterLth,data=data))
mydata.summary.th2 = mydata.summary.th[mydata.summary.th$Resp=="do",]
mydata.summary.th2$percent <- mydata.summary.th2$value/melt(xtabs(~ThemeLetterLth,data=data))$value
mydata.summary.th2$constituent <- rep(as.factor("theme"), nrow(mydata.summary.th2))
names(mydata.summary.th2) <- c("Resp", "LetterLth", "value", "percent", "constituent")

d.length <- rbind(mydata.summary2, mydata.summary.th2)


plot.length <- ggplot(d.length, aes(LetterLth, percent, group = constituent, colour=constituent)) +
  geom_smooth(aes(linetype=constituent), method="lm", size=.5, color="black") +
  geom_smooth(aes(linetype=constituent), method="lm", size=.5, color="black", fill=NA) +
  labs(title = "", y="percentage of ditransitive datives", x="length in letters") + 
  scale_linetype_manual(values=c("solid", "dashed"), name="Constituent") +
  theme(legend.position="bottom", legend.title=element_text(size=16, face="bold")) + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=14)) +
  theme(legend.text=element_text(size=16)) +
  theme(axis.title = element_text(size=16)) +
  guides(linetype=guide_legend(override.aes=list(fill=NA), ncol=2, byrow=T, title.position = "top", title.hjust = 0.5)); plot.length
  
pdf("Rplot_distribution_Length.pdf", width = 8, height = 6)
plot.length
dev.off()

#-- PLOT: proportional distributions by Complexity -------------------------
# Recipient complexity
frequencies = prop.table(with(data, table(RecBinComplexity, Resp)),1) # create prop table
frequencies <- data.frame(frequencies)
tokens = with(data, table(RecBinComplexity, Resp)) # create raw freq table
tokens <- data.frame(tokens)
tokens_sorted <- arrange(tokens, RecBinComplexity, Resp) # sort both using library(plyr)
freq_sorted <- arrange(frequencies, RecBinComplexity, Resp)

# replace every second percentage by 1 to use as geom_text
freq_sorted[1,3] <- 1
freq_sorted[3,3] <- 1

freq_sorted$Tokens <- tokens_sorted$Freq  #add raw frequency to dataset

plot.reccox = ggplot(data, aes(RecBinComplexity, fill=Resp)) + 
  geom_bar(aes(fill=Resp), position="fill",width = 0.7, color="black") + 
  labs(title = "", y="proportion of tokens", x="RecComplexity") + 
  theme(legend.position="bottom", legend.title=element_text(size=16, face="bold")) + 
  theme(axis.text.x = element_text(size=16)) +
  theme(legend.text=element_text(size=16)) +
  theme(axis.title = element_text(size=14)) +
  geom_text(data = freq_sorted, aes(x=RecBinComplexity,y=Freq, label = Tokens), vjust=2, color=ifelse(freq_sorted$Resp=="pd", "black", "white"), size=6.5) + 
  scale_fill_manual(values = cols, labels=c("ditransitive", "prepositional"), name="Variant") +
  guides(fill=guide_legend(ncol=2,title.position = "top", title.hjust = 0.5));plot.reccox


pdf("Rplot_distribution_RecComplexity.pdf", width = 4, height = 5)
plot.reccox
dev.off()

# Chisquare test
tab.reccox
prop.test(t(tab.reccox), correct=FALSE) 


# Theme animacy
# reverse levels
data$ThemeBinComplexity <- factor(data$ThemeBinComplexity, levels=c("simple", "complex"))
frequencies = prop.table(with(data, table(ThemeBinComplexity, Resp)),1) # create prop table
frequencies <- data.frame(frequencies)
tokens = with(data, table(ThemeBinComplexity, Resp)) # create raw freq table
tokens <- data.frame(tokens)
tokens_sorted <- arrange(tokens, ThemeBinComplexity, Resp) # sort both using library(plyr)
freq_sorted <- arrange(frequencies, ThemeBinComplexity, Resp)

# replace every second percentage by 1 to use as geom_text
freq_sorted[1,3] <- 1
freq_sorted[3,3] <- 1

freq_sorted$Tokens <- tokens_sorted$Freq  #add raw frequency to dataset

plot.thcox = ggplot(data, aes(ThemeBinComplexity, fill=Resp)) + 
  geom_bar(aes(fill=Resp), position="fill",width = 0.7, color="black") + 
  labs(title = "", y="proportion of tokens", x="ThemeComplexity") + 
  theme(legend.position="bottom", legend.title=element_text(size=16, face="bold")) + 
  theme(axis.text.x = element_text(size=16)) +
  theme(legend.text=element_text(size=16)) +
  theme(axis.title = element_text(size=14)) +
  geom_text(data = freq_sorted, aes(x=ThemeBinComplexity,y=Freq, label = Tokens), vjust=2, color=ifelse(freq_sorted$Resp=="pd", "black", "white"), size=6.5) + 
  scale_fill_manual(values = cols, labels=c("ditransitive", "prepositional"), name="Variant") +
  guides(fill=guide_legend(ncol=2,title.position = "top", title.hjust = 0.5));plot.thcox


pdf("Rplot_distribution_ThemeComplexity.pdf", width = 4, height = 6)
plot.thcox
dev.off()

# Chisquare test
tab.thecox
prop.test(t(tab.thecox), correct=FALSE) 

require(gridExtra)
mylegend <- g_legend(plot.thcox)


pdf("Rplot_distribution_complexity.pdf", width = 8, height = 6)
grid.arrange(arrangeGrob(plot.reccox + theme(legend.position = "none"),
                         plot.thcox + theme(legend.position = "none") + ylab(""), nrow=1), mylegend, nrow=2, heights=c(10,1))
dev.off()

# reverse ThemeAnimacy again
data$ThemeBinComplexity <- factor(data$ThemeBinComplexity, levels=c("complex", "simple"))


#-- PLOT: proportional distributions by Definiteness -------------------------
# Recipient definiteness
frequencies = prop.table(with(data, table(RecDefiniteness, Resp)),1) # create prop table
frequencies <- data.frame(frequencies)
tokens = with(data, table(RecDefiniteness, Resp)) # create raw freq table
tokens <- data.frame(tokens)
tokens_sorted <- arrange(tokens, RecDefiniteness, Resp) # sort both using library(plyr)
freq_sorted <- arrange(frequencies, RecDefiniteness, Resp)

# replace every second percentage by 1 to use as geom_text
freq_sorted[1,3] <- 1
freq_sorted[3,3] <- 1

freq_sorted$Tokens <- tokens_sorted$Freq  #add raw frequency to dataset

plot.recdef = ggplot(data, aes(RecDefiniteness, fill=Resp)) + 
  geom_bar(aes(fill=Resp), position="fill",width = 0.7, color="black") + 
  labs(title = "", y="proportion of tokens", x="RecDefiniteness") + 
  theme(legend.position="bottom", legend.title=element_text(size=16, face="bold")) + 
  theme(axis.text.x = element_text(size=16)) +
  theme(legend.text=element_text(size=16)) +
  theme(axis.title = element_text(size=14)) +
  geom_text(data = freq_sorted, aes(x=RecDefiniteness,y=Freq, label = Tokens), vjust=2, color=ifelse(freq_sorted$Resp=="pd", "black", "white"), size=6.5) + 
  scale_fill_manual(values = cols, labels=c("ditransitive", "prepositional"), name="Variant") +
  guides(fill=guide_legend(ncol=2,title.position = "top", title.hjust = 0.5));plot.recdef


pdf("Rplot_distribution_RecDefiniteness.pdf", width = 4, height = 5)
plot.recdef
dev.off()

# Chisquare test
tab.recdef
prop.test(t(tab.recdef), correct=FALSE) 


# Theme definiteness
# reverse levels
data$ThemeDefiniteness <- factor(data$ThemeDefiniteness, levels=c("def", "indef"))
frequencies = prop.table(with(data, table(ThemeDefiniteness, Resp)),1) # create prop table
frequencies <- data.frame(frequencies)
tokens = with(data, table(ThemeDefiniteness, Resp)) # create raw freq table
tokens <- data.frame(tokens)
tokens_sorted <- arrange(tokens, ThemeDefiniteness, Resp) # sort both using library(plyr)
freq_sorted <- arrange(frequencies, ThemeDefiniteness, Resp)

# replace every second percentage by 1 to use as geom_text
freq_sorted[1,3] <- 1
freq_sorted[3,3] <- 1

freq_sorted$Tokens <- tokens_sorted$Freq  #add raw frequency to dataset

plot.thdef = ggplot(data, aes(ThemeDefiniteness, fill=Resp)) + 
  geom_bar(aes(fill=Resp), position="fill",width = 0.7, color="black") + 
  labs(title = "", y="proportion of tokens", x="ThemeDefiniteness") + 
  theme(legend.position="bottom", legend.title=element_text(size=16, face="bold")) + 
  theme(axis.text.x = element_text(size=16)) +
  theme(legend.text=element_text(size=16)) +
  theme(axis.title = element_text(size=14)) +
  geom_text(data = freq_sorted, aes(x=ThemeDefiniteness,y=Freq, label = Tokens), vjust=2, color=ifelse(freq_sorted$Resp=="pd", "black", "white"), size=6.5) + 
  scale_fill_manual(values = cols, labels=c("ditransitive", "prepositional"), name="Variant") +
  guides(fill=guide_legend(ncol=2,title.position = "top", title.hjust = 0.5));plot.thdef


pdf("Rplot_distribution_ThemeDefiniteness.pdf", width = 4, height = 6)
plot.thdef
dev.off()

# Chisquare test
tab.thedef
prop.test(t(tab.thedef), correct=FALSE) 

require(gridExtra)
mylegend <- g_legend(plot.thdef)


pdf("Rplot_distribution_definiteness.pdf", width = 8, height = 6)
grid.arrange(arrangeGrob(plot.recdef + theme(legend.position = "none"),
                         plot.thdef + theme(legend.position = "none") + ylab(""), nrow=1), mylegend, nrow=2, heights=c(10,1))
dev.off()

# reverse ThemeDefiniteness again
data$ThemeDefiniteness <- factor(data$ThemeDefiniteness, levels=c("indef", "def"))

#-- PLOT: proportional distributions by Pronominality -------------------------
# Recipient pronominality
frequencies = prop.table(with(data, table(RecPron, Resp)),1) # create prop table
frequencies <- data.frame(frequencies)
tokens = with(data, table(RecPron, Resp)) # create raw freq table
tokens <- data.frame(tokens)
tokens_sorted <- arrange(tokens, RecPron, Resp) # sort both using library(plyr)
freq_sorted <- arrange(frequencies, RecPron, Resp)

# replace every second percentage by 1 to use as geom_text
freq_sorted[1,3] <- 1
freq_sorted[3,3] <- 1

freq_sorted$Tokens <- tokens_sorted$Freq  #add raw frequency to dataset

plot.recpron = ggplot(data, aes(RecPron, fill=Resp)) + 
  geom_bar(aes(fill=Resp), position="fill",width = 0.7, color="black") + 
  labs(title = "", y="proportion of tokens", x="RecPron") + 
  theme(legend.position="bottom", legend.title=element_text(size=16, face="bold")) + 
  theme(axis.text.x = element_text(size=16)) +
  theme(legend.text=element_text(size=16)) +
  theme(axis.title = element_text(size=14)) +
  geom_text(data = freq_sorted, aes(x=RecPron,y=Freq, label = Tokens), vjust=2, color=ifelse(freq_sorted$Resp=="pd", "black", "white"), size=6.5) + 
  scale_fill_manual(values = cols, labels=c("ditransitive", "prepositional"), name="Variant") +
  guides(fill=guide_legend(ncol=2,title.position = "top", title.hjust = 0.5));plot.recpron


pdf("Rplot_distribution_RecPron.pdf", width = 4, height = 6)
plot.recpron
dev.off()

# Chisquare test
tab.recpron.all <- with(data, table(RecPron, Resp))
prop.test(t(tab.recpron.all), correct=FALSE) 

# Theme pronominality
# reverse levels
data$ThemePron <- factor(data$ThemePron, levels=c("pron", "non-pron"))
frequencies = prop.table(with(data, table(ThemePron, Resp)),1) # create prop table
frequencies <- data.frame(frequencies)
tokens = with(data, table(ThemePron, Resp)) # create raw freq table
tokens <- data.frame(tokens)
tokens_sorted <- arrange(tokens, ThemePron, Resp) # sort both using library(plyr)
freq_sorted <- arrange(frequencies, ThemePron, Resp)

# replace every second percentage by 1 to use as geom_text
freq_sorted[1,3] <- 1
freq_sorted[3,3] <- 1

freq_sorted$Tokens <- tokens_sorted$Freq  #add raw frequency to dataset

plot.thpron = ggplot(data, aes(ThemePron, fill=Resp)) + 
  geom_bar(aes(fill=Resp), position="fill",width = 0.7, color="black") + 
  labs(title = "", y="proportion of tokens", x="ThemePron") + 
  theme(legend.position="bottom", legend.title=element_text(size=16, face="bold")) + 
  theme(axis.text.x = element_text(size=16)) +
  theme(legend.text=element_text(size=16)) +
  theme(axis.title = element_text(size=14)) +
  geom_text(data = freq_sorted, aes(x=ThemePron,y=Freq, label = Tokens), vjust=2, color=ifelse(freq_sorted$Resp=="pd", "black", "white"), size=6.5) + 
  scale_fill_manual(values = cols, labels=c("ditransitive", "prepositional"), name="Variant") +
  guides(fill=guide_legend(ncol=2,title.position = "top", title.hjust = 0.5));plot.thpron


pdf("Rplot_distribution_ThemePron.pdf", width = 4, height = 6)
plot.thpron
dev.off()

# Chisquare test
tab.thepron <- with(data, table(ThemePron, Resp))
prop.test(t(tab.thepron), correct=FALSE) 


require(gridExtra)
mylegend <- g_legend(plot.thpron)


pdf("Rplot_distribution_pronominality.pdf", width = 8, height = 6)
grid.arrange(arrangeGrob(plot.recpron + theme(legend.position = "none"),
                         plot.thpron + theme(legend.position = "none") + ylab(""), nrow=1), mylegend, nrow=2, heights=c(10,1))
dev.off()

# reverse ThemePron again
data$ThemePron <- factor(data$ThemePron, levels=c("non-pron", "pron"))

# Pronominality across varieties
frequencies = prop.table(with(data, table(Variety, RecPron)),1) # create prop table
frequencies <- data.frame(frequencies)
tokens = with(data, table(Variety, RecPron)) # create raw freq table
tokens <- data.frame(tokens)
tokens_sorted <- arrange(tokens, Variety, RecPron) # sort both using library(plyr)
freq_sorted <- arrange(frequencies, Variety, RecPron)

# replace every second percentage by 1 to use as geom_text
freq_sorted[1,3] <- 1
freq_sorted[3,3] <- 1
freq_sorted[5,3] <- 1
freq_sorted[7,3] <- 1
freq_sorted[9,3] <- 1
freq_sorted[11,3] <- 1
freq_sorted[13,3] <- 1
freq_sorted[15,3] <- 1
freq_sorted[17,3] <- 1


freq_sorted$Tokens <- tokens_sorted$Freq  #add raw frequency to dataset
freq_sorted$RecPron <- as.character(freq_sorted$RecPron)
freq_sorted$RecPron[freq_sorted$RecPron=="non-pron"] <- "non_pron" 
freq_sorted$RecPron <- as.factor(freq_sorted$RecPron) %>% droplevels

# change level non-pron because plotting doesn't work otherwise
datives.pron <- data
datives.pron$RecPron <- as.character(datives.pron$RecPron)
datives.pron$RecPron[datives.pron$RecPron=="non-pron"] <- "non_pron"
datives.pron$RecPron <- as.factor(datives.pron$RecPron)

cols2 <- c(pron="dimgray", non_pron="gainsboro")


plot.recpronvar = ggplot(datives.pron, aes(Variety, fill=RecPron)) +
  geom_hline(yintercept = 0.5, linetype="dashed") +
  geom_bar(aes(fill=RecPron), position="fill",width = 0.7, color="black") + 
  labs(title = "", y="proportion of tokens", x="Variety") + 
  theme(legend.position="bottom", legend.title=element_text(size=16, face="bold")) + 
  theme(axis.text.x = element_text(size=16)) +
  theme(legend.text=element_text(size=16)) +
  theme(axis.title = element_text(size=14)) +
  geom_text(data = freq_sorted, aes(x=Variety,y=rev(Freq), label = Tokens), vjust=2, color=ifelse(freq_sorted$RecPron=="non_pron", "black", "white"), size=6.5) + 
  scale_fill_manual(values = cols2, labels=c("non-pron", "pron"), name="recipient pronominality") +
  guides(fill=guide_legend(ncol=2,title.position = "top", title.hjust = 0.5));plot.recpronvar


# Chisquare test on proportions
tab.recpron <- with(data, table(RecPron, Variety))
prop.test(t(tab.recpron), correct=FALSE) 

pdf("Rplot_distribution_RecPronByVar.pdf", width = 8, height = 5)
plot.recpronvar
dev.off()

#-- PLOT: proportional distributions by Givenness -------------------------
# Recipient givenness
frequencies = prop.table(with(data, table(RecGivenness, Resp)),1) # create prop table
frequencies <- data.frame(frequencies)
tokens = with(data, table(RecGivenness, Resp)) # create raw freq table
tokens <- data.frame(tokens)
tokens_sorted <- arrange(tokens, RecGivenness, Resp) # sort both using library(plyr)
freq_sorted <- arrange(frequencies, RecGivenness, Resp)

# replace every second percentage by 1 to use as geom_text
freq_sorted[1,3] <- 1
freq_sorted[3,3] <- 1

freq_sorted$Tokens <- tokens_sorted$Freq  #add raw frequency to dataset

plot.recgiv = ggplot(data, aes(RecGivenness, fill=Resp)) + 
  geom_bar(aes(fill=Resp), position="fill",width = 0.7, color="black") + 
  labs(title = "", y="proportion of tokens", x="RecGivenness") + 
  theme(legend.position="bottom", legend.title=element_text(size=16, face="bold")) + 
  theme(axis.text.x = element_text(size=16)) +
  theme(legend.text=element_text(size=16)) +
  theme(axis.title = element_text(size=14)) +
  geom_text(data = freq_sorted, aes(x=RecGivenness,y=Freq, label = Tokens), vjust=2, color=ifelse(freq_sorted$Resp=="pd", "black", "white"), size=6.5) + 
  scale_fill_manual(values = cols, labels=c("ditransitive", "prepositional"), name="Variant") +
  guides(fill=guide_legend(ncol=2,title.position = "top", title.hjust = 0.5));plot.recgiv


pdf("Rplot_distribution_RecGivenness.pdf", width = 4, height = 5)
plot.recgiv
dev.off()

# Chisquare test
tab.recgiv
prop.test(t(tab.recgiv), correct=FALSE) 

# Theme givenness
# reverse levels
data$ThemeGivenness <- factor(data$ThemeGivenness, levels=c("given", "new"))
frequencies = prop.table(with(data, table(ThemeGivenness, Resp)),1) # create prop table
frequencies <- data.frame(frequencies)
tokens = with(data, table(ThemeGivenness, Resp)) # create raw freq table
tokens <- data.frame(tokens)
tokens_sorted <- arrange(tokens, ThemeGivenness, Resp) # sort both using library(plyr)
freq_sorted <- arrange(frequencies, ThemeGivenness, Resp)

# replace every second percentage by 1 to use as geom_text
freq_sorted[1,3] <- 1
freq_sorted[3,3] <- 1

freq_sorted$Tokens <- tokens_sorted$Freq  #add raw frequency to dataset

plot.thgiv = ggplot(data, aes(ThemeGivenness, fill=Resp)) + 
  geom_bar(aes(fill=Resp), position="fill",width = 0.7, color="black") + 
  labs(title = "", y="proportion of tokens", x="ThemeGivenness") + 
  theme(legend.position="bottom", legend.title=element_text(size=16, face="bold")) + 
  theme(axis.text.x = element_text(size=16)) +
  theme(legend.text=element_text(size=16)) +
  theme(axis.title = element_text(size=14)) +
  geom_text(data = freq_sorted, aes(x=ThemeGivenness,y=Freq, label = Tokens), vjust=2, color=ifelse(freq_sorted$Resp=="pd", "black", "white"), size=6.5) + 
  scale_fill_manual(values = cols, labels=c("ditransitive", "prepositional"), name="Variant") +
  guides(fill=guide_legend(ncol=2,title.position = "top", title.hjust = 0.5));plot.thgiv


pdf("Rplot_distribution_ThemeGivenness.pdf", width = 4, height = 5)
plot.thgiv
dev.off()

# Chisquare test
tab.thegiv
prop.test(t(tab.thegiv), correct=FALSE) 

require(gridExtra)
mylegend <- g_legend(plot.thgiv)


pdf("Rplot_distribution_givenness.pdf", width = 8, height = 6)
grid.arrange(arrangeGrob(plot.recgiv + theme(legend.position = "none"),
                         plot.thgiv + theme(legend.position = "none") + ylab(""), nrow=1), mylegend, nrow=2, heights=c(10,1))
dev.off()

# reverse ThemePron again
data$ThemeGivenness <- factor(data$ThemeGivenness, levels=c("new", "given"))

#-- PLOT: proportional distributions by VerbSemantics -------------------------
# VerbSemantics
# change level
data$VerbSemantics <- factor(data$VerbSemantics, levels=c("t", "f", "c", "p", "a"))
data$VerbSemantics_sp <- data$VerbSemantics

levels(data$VerbSemantics_sp) <- c("transfer", "future t.", "communication", "prevention", "abstract")

frequencies = prop.table(with(data, table(VerbSemantics_sp, Resp)),1) # create prop table
frequencies <- data.frame(frequencies)
tokens = with(data, table(VerbSemantics_sp, Resp)) # create raw freq table
tokens <- data.frame(tokens)
tokens_sorted <- arrange(tokens, VerbSemantics_sp, Resp) # sort both using library(plyr)
freq_sorted <- arrange(frequencies, VerbSemantics_sp, Resp)

# replace every second percentage by 1 to use as geom_text
freq_sorted[1,3] <- 1
freq_sorted[3,3] <- 1
freq_sorted[5,3] <- 1
freq_sorted[7,3] <- 1
freq_sorted[9,3] <- 1

freq_sorted$Tokens <- tokens_sorted$Freq  #add raw frequency to dataset

plot.verbsense = ggplot(data, aes(VerbSemantics_sp, fill=Resp)) + 
  geom_bar(aes(fill=Resp), position="fill",width = 0.7, color="black") + 
  labs(title = "", y="proportion of tokens", x="") + 
  theme(legend.position="bottom", legend.title=element_text(size=16, face="bold")) + 
  theme(axis.text.x = element_text(size=16)) +
  theme(legend.text=element_text(size=16)) +
  theme(axis.title = element_text(size=14)) +
  geom_text(data = freq_sorted, aes(x=VerbSemantics_sp,y=Freq, label = Tokens), vjust=2, color=ifelse(freq_sorted$Resp=="pd", "black", "white"), size=6.5) + 
  scale_fill_manual(values = cols, labels=c("ditransitive", "prepositional"), name="Variant") +
  guides(fill=guide_legend(ncol=2,title.position = "top", title.hjust = 0.5));plot.verbsense


pdf("Rplot_distribution_VerbSemantics.pdf", width = 9, height = 5)
plot.verbsense
dev.off()

data$VerbSemantics <- factor(data$VerbSemantics, levels=c("a", "c", "f", "p", "t"))

# test for X-square
tab.verbsem <- with(data, table(VerbSemantics, Resp))
prop.test(tab.verbsem, correct=FALSE) 



#-- PLOT: proportional distributions by PrimeType -------------------------
frequencies = prop.table(with(data, table(PrimeTypePruned, Resp)),1) # create prop table
frequencies <- data.frame(frequencies)
tokens = with(data, table(PrimeTypePruned, Resp)) # create raw freq table
tokens <- data.frame(tokens)
tokens_sorted <- arrange(tokens, PrimeTypePruned, Resp) # sort both using library(plyr)
freq_sorted <- arrange(frequencies, PrimeTypePruned, Resp)

# replace every second percentage by 1 to use as geom_text
freq_sorted[1,3] <- 1
freq_sorted[3,3] <- 1
freq_sorted[5,3] <- 1

freq_sorted$Tokens <- tokens_sorted$Freq  #add raw frequency to dataset

plot.prime = ggplot(data, aes(PrimeTypePruned, fill=Resp)) + 
  geom_bar(aes(fill=Resp), position="fill",width = 0.7, color="black") + 
  labs(title = "", y="proportion of tokens", x="") + 
  theme(legend.position="bottom", legend.title=element_text(size=16, face="bold")) + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16)) +
  theme(legend.text=element_text(size=16)) +
  theme(axis.title = element_text(size=14)) +
  geom_text(data = freq_sorted, aes(x=PrimeTypePruned,y=Freq, label = Tokens), vjust=2, color=ifelse(freq_sorted$Resp=="pd", "black", "white"), size=6.5) + 
  scale_fill_manual(values = cols, labels=c("ditransitive", "prepositional"), name="Variant") +
  guides(fill=guide_legend(ncol=2,title.position = "top", title.hjust = 0.5));plot.prime


pdf("Rplot_distribution_PrimeType.pdf", width = 6, height = 5)
plot.prime
dev.off()

# test for X-square
tab.prime <- with(data, table(PrimeTypePruned, Resp))
prop.test(tab.prime, correct=FALSE) 



#-- PLOT: proportional distributions by Frequency -------------------------

freq.df <- data[,c("RecHeadFreq", "Resp")]
freq.df$Constituent <- as.factor(rep("recipient", nrow(freq.df)))
names(freq.df) <- c("Frequency", "Resp", "Constituent")

freq.df.th <- data[,c("ThemeHeadFreq", "Resp")]
freq.df.th$Constituent <- as.factor(rep("theme", nrow(freq.df.th)))
names(freq.df.th) <- c("Frequency", "Resp", "Constituent")

freq.df.all <- rbind(freq.df, freq.df.th)

# get sd in order of boxplots
sd.all <- c(sd(data[data$Resp=="do", "RecHeadFreq"]), sd(data[data$Resp=="do", "ThemeHeadFreq"]), sd(data[data$Resp=="pd", "RecHeadFreq"]), sd(data[data$Resp=="pd", "ThemeHeadFreq"]))
sd.all <- as.data.frame(sd.all)
sd.all$Resp <- as.factor(c("do", "do", "pd", "pd"))
sd.all$Constituent <- as.factor(c("recipient", "theme", "recipient", "theme"))
sd.all$sd.all <- round(sd.all$sd.all, 1)


mean.all <- c(mean(data[data$Resp=="do", "RecHeadFreq"]), mean(data[data$Resp=="do", "ThemeHeadFreq"]), mean(data[data$Resp=="pd", "RecHeadFreq"]), mean(data[data$Resp=="pd", "ThemeHeadFreq"]))
mean.all <- as.data.frame(mean.all)
mean.all$Resp <- as.factor(c("do", "do", "pd", "pd"))
mean.all$Constituent <- as.factor(c("recipient", "theme", "recipient", "theme"))
mean.all$mean.all <- round(mean.all$mean.all, 1)

cols3 <- c(recipient="dimgray", theme="gainsboro")

plot.freq <- ggplot(freq.df.all, aes(Resp, Frequency, fill=Constituent)) +
  geom_boxplot(position=position_dodge(0.8), width=0.7, outlier.alpha = 0.3, outlier.size = 0.8) +
  geom_text(data=sd.all, aes(label=paste("SD =", sd.all), y=max(freq.df.all$Frequency) + 1000),position=position_dodge(0.8), size=4, fontface="italic") +
  geom_text(data=mean.all, aes(label=paste("mean =", mean.all), y=max(freq.df.all$Frequency) + 2000),position=position_dodge(0.8), size=4, fontface="italic") +
  labs(title = "", y="normalized frequency", x="") + 
  theme(legend.position="bottom", legend.title=element_text(size=16, face="bold")) + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16)) +
  theme(legend.text=element_text(size=16)) +
  theme(axis.title = element_text(size=14)) +
  scale_x_discrete(breaks=c("do","pd"),
                   labels=c("ditransitive", "prepositional")) +
  scale_fill_manual(values = cols3, labels=c("recipient", "theme"), name="Constituent") +
  guides(fill=guide_legend(ncol=2,title.position = "top", title.hjust = 0.5));plot.freq

pdf("Rplot_distribution_Freq.pdf", width = 9, height = 5)
plot.freq
dev.off()

# t-test to compare means
th.test <- t.test(freq.df.th[freq.df.th$Resp=="pd", "Frequency"], freq.df.th[freq.df.th$Resp=="do", "Frequency"])
str(th.test)
p.th <- th.test$p.value

rec.test <- t.test(freq.df[freq.df$Resp=="do", "Frequency"], freq.df[freq.df$Resp=="pd", "Frequency"])
str(rec.test)
p.rec <- rec.test$p.value

# make same plot as with length
library(reshape2)
mydata.summary = melt(xtabs(~Resp+RecHeadFreq,data=data))
mydata.summary2 = mydata.summary[mydata.summary$Resp=="do",]
mydata.summary2$percent <- mydata.summary2$value/melt(xtabs(~RecHeadFreq,data=data))$value
mydata.summary2$constituent <- rep(as.factor("recipient"), nrow(mydata.summary2))
names(mydata.summary2) <- c("Resp", "HeadFreq", "value", "percent", "constituent")

mydata.summary.th = melt(xtabs(~Resp+ThemeHeadFreq,data=data))
mydata.summary.th2 = mydata.summary.th[mydata.summary.th$Resp=="do",]
mydata.summary.th2$percent <- mydata.summary.th2$value/melt(xtabs(~ThemeHeadFreq,data=data))$value
mydata.summary.th2$constituent <- rep(as.factor("theme"), nrow(mydata.summary.th2))
names(mydata.summary.th2) <- c("Resp", "HeadFreq", "value", "percent", "constituent")

d.freq <- rbind(mydata.summary2, mydata.summary.th2)


plot.freq2 <- ggplot(d.freq, aes(HeadFreq, percent, group = constituent, colour=constituent)) +
  geom_smooth(aes(linetype=constituent), method="lm", size=.5, color="black") +
  geom_smooth(aes(linetype=constituent), method="lm", size=.5, color="black", fill=NA) +
  labs(title = "", y="percentage of ditransitive datives", x="normalized frequency") + 
  scale_linetype_manual(values=c("solid", "dashed"), name="constituent") +
  theme(legend.position="bottom", legend.title=element_text(size=16, face="bold")) + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=14)) +
  theme(legend.text=element_text(size=16)) +
  theme(axis.title = element_text(size=16)) +
  guides(linetype=guide_legend(override.aes=list(fill=NA), ncol=2, byrow=T, title.position = "top", title.hjust = 0.5)); plot.freq2

pdf("Rplot_distribution_freq2.pdf", width = 8, height = 6)
plot.freq2
dev.off()


#-- PLOT: proportional distributions by Thematicity -------------------------

thema <- data[,c("RecThematicity", "Resp")]
thema$Constituent <- as.factor(rep("recipient", nrow(thema)))
names(thema) <- c("Thematicity", "Resp", "Constituent")

thema.th <- data[,c("ThemeThematicity", "Resp")]
thema.th$Constituent <- as.factor(rep("theme", nrow(thema.th)))
names(thema.th) <- c("Thematicity", "Resp", "Constituent")

thema.all <- rbind(thema, thema.th)


min.mean.sd.max <- function(x) {
  r <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

# get sd in order of boxplots
sd.all <- c(sd(data[data$Resp=="do", "RecThematicity"]), sd(data[data$Resp=="do", "ThemeThematicity"]), sd(data[data$Resp=="pd", "RecThematicity"]), sd(data[data$Resp=="pd", "ThemeThematicity"]))
sd.all <- as.data.frame(sd.all)
sd.all$Resp <- as.factor(c("do", "do", "pd", "pd"))
sd.all$Constituent <- as.factor(c("recipient", "theme", "recipient", "theme"))
sd.all$sd.all <- round(sd.all$sd.all, 3)


mean.all <- c(mean(data[data$Resp=="do", "RecThematicity"]), mean(data[data$Resp=="do", "ThemeThematicity"]), mean(data[data$Resp=="pd", "RecThematicity"]), mean(data[data$Resp=="pd", "ThemeThematicity"]))
mean.all <- as.data.frame(mean.all)
mean.all$Resp <- as.factor(c("do", "do", "pd", "pd"))
mean.all$Constituent <- as.factor(c("recipient", "theme", "recipient", "theme"))
mean.all$mean.all <- round(mean.all$mean.all, 3)

cols3 <- c(recipient="dimgray", theme="gainsboro")

plot.thematicity <- ggplot(thema.all, aes(Resp, Thematicity, fill=Constituent)) +
  geom_boxplot(position=position_dodge(0.8), width=0.7, outlier.alpha = 0.3, outlier.size = 0.8) +
  geom_text(data=sd.all, aes(label=paste("SD =", sd.all), y=max(thema.all$Thematicity) - 0.03),position=position_dodge(0.8), size=4, fontface="italic") +
  geom_text(data=mean.all, aes(label=paste("mean =", mean.all), y=max(thema.all$Thematicity) - 0.02),position=position_dodge(0.8), size=4, fontface="italic") +
  labs(title = "", y="normalized thematicity", x="") + 
  theme(legend.position="bottom", legend.title=element_text(size=16, face="bold")) + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16)) +
  theme(legend.text=element_text(size=16)) +
  theme(axis.title = element_text(size=14)) +
  scale_x_discrete(breaks=c("do","pd"),
                   labels=c("ditransitive", "prepositional")) +
  scale_fill_manual(values = cols3, labels=c("recipient", "theme"), name="Constituent") +
  guides(fill=guide_legend(ncol=2,title.position = "top", title.hjust = 0.5));plot.thematicity

pdf("Rplot_distribution_Thematicity.pdf", width = 7, height = 5)
plot.thematicity
dev.off()


# t-test to check for significance
plot(density(thema[,1])) # long tail but lots of data

rec.test <- t.test(thema[thema$Resp=="do", "Thematicity"], thema[thema$Resp=="pd", "Thematicity"])
str(rec.test)
p.rec <- rec.test$p.value

th.test <- t.test(thema.th[thema.th$Resp=="pd", "Thematicity"], thema.th[thema.th$Resp=="do", "Thematicity"])
str(th.test)
p.th <- th.test$p.value

#-- PLOT: proportional distributions by TTR -------------------------
# make same plot as with length
library(reshape2)
mydata.summary = melt(xtabs(~Resp+TypeTokenRatio,data=data))
d.ttr = mydata.summary[mydata.summary$Resp=="do",]
d.ttr$percent <- d.ttr$value/melt(xtabs(~TypeTokenRatio,data=data))$value
names(d.ttr) <- c("Resp", "TTR", "value", "percent")


plot.ttr <- ggplot(d.ttr, aes(TTR, percent)) +
  geom_smooth(method="lm", size=.5, color="black") +
  geom_smooth(method="lm", size=.5, color="black", fill=NA) +
  labs(title = "", y="percentage of ditransitive datives", x="type-token ratio") + 
  theme(legend.position="bottom", legend.title=element_text(size=16, face="bold")) + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=14)) +
  theme(legend.text=element_text(size=16)) +
  theme(axis.title = element_text(size=16)) +
  guides(linetype=guide_legend(override.aes=list(fill=NA), title.position = "top", title.hjust = 0.5)); plot.ttr

pdf("Rplot_distribution_ttr.pdf", width = 8, height = 6)
plot.ttr
dev.off()


# -- create pic for Manhattan und Euclidean distance------

df <- data.frame(x=1:2, y=1:2)

plot.euclid <- ggplot(df, aes(x,y)) +
  geom_point(size=4, col="black", fill="white", shape=21) +
  geom_text(label=c("Var A", "Var B"), nudge_y =0.12, nudge_x=-0.06, fontface="bold", size=5) +
  ylim(0.5,2.5) +
  xlim(0.5,2.5) +
  geom_path(arrow=arrow(length=unit(0.3, "cm"), type="closed", ends="first"), size=1) +
  geom_segment(aes(x = 2, y = 1, xend = 2, yend = 2), arrow=arrow(length=unit(0.3, "cm"), type = "closed"), linetype="dashed", size=1) +
  geom_segment(aes(x = 1, y = 1, xend = 2, yend = 1), linetype="dashed", size=1) +
  annotate("text", x=1.8, y=0.9, label="Manhattan distance", size=5) +
  annotate("text", x=1.2, y=1.6, label="Euclidean \n distance", size=5);plot.euclid

pdf("Rplot_euclidean_distance.pdf", width = 5, height = 4)
plot.euclid
dev.off()




# -- find super-token for example -------------
dotest <- data[data$Resp=="do", c(21,29,32)]
pdtest <- data[data$Resp=="pd", c(21,29,32)]
rownames(dotest) <- NULL
rownames(pdtest) <- NULL

dotest$Resp <- "DO"
pdtest$Resp <- "PD"

test <- rbind(dotest, pdtest)
dupRows <- dupsBetweenGroups(test, "Resp")
test <- cbind(test, dup=dupRows)

test[test$dup=="TRUE", ]

data[data$RecHeadLemma=="police" & data$ThemeHeadLemma == 'statement' & data$Verb=="give", "SpeakerID"] # no. 6 & 7, JA:S1B-069



save.image("datives_descriptive.RData")
