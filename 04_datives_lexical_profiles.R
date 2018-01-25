#=======================================================================
#   Script: 04_datives_lexical_profiles.R
#   Date last changed: 10 January 2018
#
#   Rversion: 3.3.2
#   Author: MR
#   This script calculates association strength between collexemes
#   (covarying collexeme analysis) and between lexical items and the
#   constructions (DO vs. PD) they occur in (distinctive collexeme
#   analysis). It makes use of Gries' 2014 Rscript.
#=======================================================================

#---SETUP: LIBRARIES -----------------------------------------

library(reshape2);library(plyr);library(dplyr);library(magrittr);library(lme4);library(ggplot2);library(gridExtra); library(JGmisc); library(beepr); library(effects); library(Hmisc); library(car); library(Rmisc)

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
cols <- c(do="dimgray",pd="gainsboro")

#-----SETUP: DATA------------------------
data <- read.delim("OA/datives.txt")
dats <- data

# get the varieties
vars <- levels(dats$Variety)

# rename labels
levels(dats$Resp) = c("ditransitive", "prepositional") 

#------- COVARYING COLLEXEME ANALYSIS -----------------------------------------

# create list of dataframes for each Variety and response
data.list = list()
for (i in seq(2,18,2)){
	v <- vars[(i/2)]
	sub.data <- droplevels(subset(dats, Variety == v)) # drop unused levels
	sub.data <- sub.data[, c("Verb", "VerbThemeLemma", "RecHeadLemma", "ThemeHeadLemma", "Resp")] # keep dataframes small -> leave out unnecessary cols
	sub.data.d <- droplevels(subset(sub.data, Resp == "ditransitive"))
	sub.data.p <- droplevels(subset(sub.data, Resp == "prepositional"))
	# assign name, e.g. "CAN.data"
	data.list[[i - 1]] <- sub.data.d; names(data.list)[i - 1] <- paste(v, "d.data", sep = '.')
	data.list[[i]] <- sub.data.p; names(data.list)[i] <- paste(v, "p.data", sep = '.')
	rm(list = c("sub.data.d", "sub.data.p", "v")) # remove temporary stuff
}

# Rerun collexeme list separately by constituent combination and then plot, before calculating the next collexeme list since the same variable name (collexeme.list) is used.

# -- .. recipient + theme --------------------
# make list of collexeme associations for each var and resp with the rec and th
collexeme.list <- list()
for (i in 1:18){
	cur.data <- data.list[[i]]
	out.df <- data.frame() # blank output dataframe
	for (j in 1:nrow(cur.data)){
		r <- cur.data[j, "RecHeadLemma"]; freq.r <- length(which(cur.data$RecHeadLemma == r))
		t <- cur.data[j, "ThemeHeadLemma"]; freq.t <- length(which(cur.data$ThemeHeadLemma == t))
		freq.rt <- length(which(cur.data$RecHeadLemma == r & cur.data$ThemeHeadLemma == t))
		# create frequency table from Stefanowitsch & Gries (2005: 9)
		tab1 <- matrix(c(freq.rt,
										 freq.t - freq.rt,
										 freq.r - freq.rt,
										 nrow(cur.data) - (freq.r + freq.t - freq.rt)),
									 ncol = 2)
		# calculate p-value and log10 p-value
		Ft.pvalue <- fisher.test(tab1)$p.value
		log.pvalue <- log10(Ft.pvalue)
		# adjust sign. If observed freq is greater than expected then make log10p-value positive
		expected <- round((addmargins(tab1)[3,1]/sum(tab1)) * (addmargins(tab1)[1,3]/sum(tab1)) * sum(tab1))
		if (tab1[1, 1] > expected) {
			log.pvalue <- abs(log.pvalue) 
		}
		# add to output df
		out.df <- rbind(out.df,
										data.frame(pair = paste(r, t),  # left out trim() in "trim(paste(v,h))"
															 p.value = Ft.pvalue, 
															 log.pvalue = log.pvalue))
	}
	out.df <- unique(out.df) # remove duplicates
	out.df <- out.df[order(out.df$log.pvalue, decreasing = T), ] # sort
	collexeme.list[[i]] <- out.df # add to list and fix the name
	names(collexeme.list)[i] <- paste(unlist(strsplit(names(data.list)[i], split = '[.]'))[1:2], collapse = ".")
	rm(list = c("out.df", "cur.data", "r", "j", "tab1"))
}



# -- .. verb + theme --------------------
# make list of collexeme associations for each var and resp with verb + theme
collexeme.list <- list()
for (i in 1:18){
  cur.data <- data.list[[i]]
  out.df <- data.frame() # blank output dataframe
  for (j in 1:nrow(cur.data)){
    v <- cur.data[j, "Verb"]; freq.v <- length(which(cur.data$Verb == v))
    t <- cur.data[j, "ThemeHeadLemma"]; freq.t <- length(which(cur.data$ThemeHeadLemma == t))
    freq.vt <- length(which(cur.data$Verb == v & cur.data$ThemeHeadLemma == t))
    # create frequency table from Stefanowitsch & Gries (2005: 9)
    tab1 <- matrix(c(freq.vt,
                     freq.t - freq.vt,
                     freq.v - freq.vt,
                     nrow(cur.data) - (freq.v + freq.t - freq.vt)),
                   ncol = 2)
    # calculate p-value and log10 p-value
    Ft.pvalue <- fisher.test(tab1)$p.value
    log.pvalue <- log10(Ft.pvalue)
    # adjust sign. If observed freq is greater than expected then make log10p-value positive
    expected <- round((addmargins(tab1)[3,1]/sum(tab1)) * (addmargins(tab1)[1,3]/sum(tab1)) * sum(tab1))
    if (tab1[1, 1] > expected) {
      log.pvalue <- abs(log.pvalue) 
    }
    # add to output df
    out.df <- rbind(out.df,
                    data.frame(pair = paste(v, t),  # left out trim() in "trim(paste(v,h))"
                               p.value = Ft.pvalue, 
                               log.pvalue = log.pvalue))
  }
  out.df <- unique(out.df) # remove duplicates
  out.df <- out.df[order(out.df$log.pvalue, decreasing = T), ] # sort
  collexeme.list[[i]] <- out.df # add to list and fix the name
  names(collexeme.list)[i] <- paste(unlist(strsplit(names(data.list)[i], split = '[.]'))[1:2], collapse = ".")
  rm(list = c("out.df", "cur.data", "t", "j", "tab1"))
}

# -- .. verb + recipient --------------------
# make list of collexeme associations for each var and resp with verb + recipient
collexeme.list <- list()
for (i in 1:18){
  cur.data <- data.list[[i]]
  out.df <- data.frame() # blank output dataframe
  for (j in 1:nrow(cur.data)){
    v <- cur.data[j, "Verb"]; freq.v <- length(which(cur.data$Verb == v))
    r <- cur.data[j, "RecHeadLemma"]; freq.r <- length(which(cur.data$RecHeadLemma == r))
    freq.vr <- length(which(cur.data$Verb == v & cur.data$RecHeadLemma == r))
    # create frequency table from Stefanowitsch & Gries (2005: 9)
    tab1 <- matrix(c(freq.vr,
                     freq.r - freq.vr,
                     freq.v - freq.vr,
                     nrow(cur.data) - (freq.v + freq.r - freq.vr)),
                   ncol = 2)
    # calculate p-value and log10 p-value
    Ft.pvalue <- fisher.test(tab1)$p.value
    log.pvalue <- log10(Ft.pvalue)
    # adjust sign. If observed freq is greater than expected then make log10p-value positive
    expected <- round((addmargins(tab1)[3,1]/sum(tab1)) * (addmargins(tab1)[1,3]/sum(tab1)) * sum(tab1))
    if (tab1[1, 1] > expected) {
      log.pvalue <- abs(log.pvalue) 
    }
    # add to output df
    out.df <- rbind(out.df,
                    data.frame(pair = paste(v, r),  # left out trim() in "trim(paste(v,h))"
                               p.value = Ft.pvalue, 
                               log.pvalue = log.pvalue))
  }
  out.df <- unique(out.df) # remove duplicates
  out.df <- out.df[order(out.df$log.pvalue, decreasing = T), ] # sort
  collexeme.list[[i]] <- out.df # add to list and fix the name
  names(collexeme.list)[i] <- paste(unlist(strsplit(names(data.list)[i], split = '[.]'))[1:2], collapse = ".")
  rm(list = c("out.df", "cur.data", "r", "j", "tab1"))
}



#------------ results ---------------------------------------------------------------

# look at the top 10 strongest attracted collexemes
for (i in 1:18){
	print(names(collexeme.list)[i])
	print(head(collexeme.list[[i]], 10))
	cat("\n\n")
}

# look at the distribution of collexeme strengths across VoEs
for (i in 1:18){
	print(names(collexeme.list)[i])
	print(summary(collexeme.list[[i]]["log.pvalue"]))
	cat("\n\n")
}

# -- extract and combine the dataframes ------------------
# ditransitive Cxs
jd <- data.frame()
for (i in seq(1, 17, 2)){
  var <- unlist(strsplit(names(collexeme.list)[i], split = '[.]'))[1]
  cur.df <- collexeme.list[[i]]
  cur.df$Variety <- var
  jd <- rbind(jd, cur.df)
}	

# prepositional Cxs
spld <- data.frame()
for (i in seq(2, 18, 2)){
  var <- unlist(strsplit(names(collexeme.list)[i], split = '[.]'))[1]
  cur.df <- collexeme.list[[i]]
  cur.df$Variety <- var
  spld <- rbind(spld, cur.df)
}	

# combine datasets
combined <- rbind(jd, spld)
combined$Cx <- c(rep("ditransitive", nrow(jd)), rep("prepositional", nrow(spld)))


#--- plot the means ---------------------------------------------------

# test:
colMeans(collexeme.list[[1]]["log.pvalue"])

# calculate the means
#vars <- levels(as.factor(combined$Variety))
dtlist.mean <- matrix(ncol=2, nrow=9)
Resp <- factor(c("ditransitive", "prepositional"))
for (i in seq(2,18,2)){
  mean.p <- round(colMeans(collexeme.list[[i]]["log.pvalue"]),2)
  mean.d <- round(colMeans(collexeme.list[[i-1]]["log.pvalue"]),2)
  dtlist.mean[i/2,2] <- mean.p
  dtlist.mean[i/2,1] <- mean.d
  
}
rownames(dtlist.mean) <- vars
colnames(dtlist.mean) <- Resp
dtlist.mean <- as.data.frame(dtlist.mean)
dtlist.mean$Variety <- rownames(dtlist.mean)
dtlist.mean.melted <- melt(dtlist.mean)
dtlist.mean.melted$Variety <- factor(dtlist.mean.melted$Variety, levels=c("CAN", "GB", "IRE", "NZ","JA", "SIN", "HK", "IND", "PHI"))

# extract mean with other stats as well (ci, se)
dtlist.mean.all <- summarySE(combined, measurevar = "log.pvalue", groupvars = c("Variety", "Cx"))
dtlist.mean.all$Variety <- factor(dtlist.mean.all$Variety, levels=c("CAN", "GB", "IRE", "NZ","JA", "SIN", "HK", "IND", "PHI"))


# barplot of mean coll strength (df.list2)
cols <- c(ditransitive="dimgray",prepositional="gainsboro")

# Make barplot: use barplot since boxplot is not much informative
# meanplot.verbrec: verb-recipient
# meanplot.verbth: verb-theme
# meanplot.recth: recipient-theme
(meanplot.verbrec <- ggplot(dtlist.mean.all, aes(Variety, log.pvalue, fill=Cx)) +
  geom_bar(stat="identity", position=position_dodge(), color="black", width=.8) +
  geom_text(aes(y=log.pvalue+ci+0.1, label=round(log.pvalue, 2)), position=position_dodge(0.8), size=5. ) +
  geom_errorbar(aes(ymin=log.pvalue-ci, ymax=log.pvalue+ci), colour="black", width=.1, position=position_dodge(width=0.8), linetype=1) +
  scale_fill_manual(values=cols, name="Variant")+
  labs(title="verb-recipient associations", y="mean collstr. strength", x="") +
  theme(axis.title.x = element_text(size=12), axis.title.y = element_text(size=15)) +
  theme(axis.text.x = element_text(size=16), plot.title = element_text(size=16, face="bold", hjust = 0.5)) +
  guides(fill=guide_legend(ncol=2,title.position = "top", title.hjust = 0.5))+
  theme(legend.position = "bottom",legend.title = element_text(size=16,face = "bold"), legend.text = element_text(size=16)))

pdf("Rplot_lexical_collexeme_mean_recth.pdf", width = 12, height = 7)
meanplot.recth
dev.off()

pdf("Rplot_lexical_collexeme_mean_verbth.pdf", width = 12, height = 7)
meanplot.verbth
dev.off()

pdf("Rplot_lexical_collexeme_mean_verbrec.pdf", width = 12, height = 7)
meanplot.verbrec
dev.off()


require(gridExtra)
mylegend <- g_legend(meanplot.verbth)

grid.arrange(arrangeGrob(meanplot.verbrec + theme(legend.position = "none"),
            meanplot.verbth + theme(legend.position = "none"),
            meanplot.recth + theme(legend.position = "none"),
            mylegend,
            heights = c(0.3, 0.3, 0.3, 0.1), ncol=1))



pdf("Rplot_lexical_collexeme_mean_allplots.pdf", width = 12, height = 20)
grid.arrange(arrangeGrob(meanplot.verbrec + theme(legend.position = "none"),
                         meanplot.verbth + theme(legend.position = "none"),
                         meanplot.recth + theme(legend.position = "none"),
                         mylegend,
                         heights = c(0.28, 0.28, 0.28, 0.07), ncol=1))
dev.off()





# -- t-test to compare means ------------------------------------------
head(combined)
# create subsets by variant and variety to compare
length(as.vector(collexeme.list$PHI.p[,3])) # PhiE, prep, log.pvalue

# is the data normally distributed?
plot(density(combined$log.pvalue))
plot(density(collexeme.list$CAN.d[,3]))
plot(density(collexeme.list$CAN.p[,3]))
plot(density(collexeme.list$GB.d[,3]))
plot(density(collexeme.list$GB.p[,3]))
plot(density(collexeme.list$HK.d[,3]))
plot(density(collexeme.list$HK.p[,3]))
plot(density(collexeme.list$IND.d[,3]))
plot(density(collexeme.list$IND.p[,3]))
plot(density(collexeme.list$IRE.d[,3]))
plot(density(collexeme.list$IRE.p[,3]))
plot(density(collexeme.list$JA.d[,3]))
plot(density(collexeme.list$JA.p[,3]))
plot(density(collexeme.list$NZ.p[,3]))
plot(density(collexeme.list$NZ.d[,3]))
plot(density(collexeme.list$PHI.d[,3]))
plot(density(collexeme.list$PHI.p[,3]))
plot(density(collexeme.list$SIN.p[,3]))
plot(density(collexeme.list$SIN.d[,3]))

qqnorm(collexeme.list$CAN.d[,3])
qqline(collexeme.list$CAN.d[,3])
qqnorm(collexeme.list$IND.d[,3])
qqline(collexeme.list$IND.d[,3])

shapiro.test(collexeme.list$CAN.d[,3])
shapiro.test(collexeme.list$CAN.p[,3])
shapiro.test(collexeme.list$GB.d[,3])
shapiro.test(collexeme.list$GB.p[,3])

# verb-rec: largely normally distributed but enough data points for t.test
# verb-th: to some extent normally, but looks a bit binomial and very long tails
# rec-th: long tails and bumps in the curve; ditrans and prep pattern similar across varieties

# run t-tests on all separately by variety, correct for Bonferroni
test.can <- t.test(collexeme.list$CAN.d[,3], collexeme.list$CAN.p[,3])
str(test.can)
pCAN <- test.can$p.value
test.gb <- t.test(collexeme.list$GB.d[,3], collexeme.list$GB.p[,3])
pGB <- test.gb$p.value
test.hk <- t.test(collexeme.list$HK.d[,3], collexeme.list$HK.p[,3])
pHK <- test.hk$p.value
test.ind <- t.test(collexeme.list$IND.d[,3], collexeme.list$IND.p[,3])
pIND <- test.ind$p.value
test.ire <- t.test(collexeme.list$IRE.d[,3], collexeme.list$IRE.p[,3])
pIRE <- test.ire$p.value
test.ja <- t.test(collexeme.list$JA.d[,3], collexeme.list$JA.p[,3])
pJA <- test.ja$p.value
test.nz <- t.test(collexeme.list$NZ.d[,3], collexeme.list$NZ.p[,3])
pNZ <- test.nz$p.value
test.phi <- t.test(collexeme.list$PHI.d[,3], collexeme.list$PHI.p[,3])
pPHI <- test.phi$p.value
test.sin <- t.test(collexeme.list$SIN.d[,3], collexeme.list$SIN.p[,3])
pSIN <- test.sin$p.value

p.adjust(c(pCAN, pGB, pHK, pIND, pIRE, pJA, pNZ, pPHI, pSIN), method = "bonferroni")

# compare mean ditrans IND vs global average (for verb-rec)
# get global subset -IND
vs.ind <- droplevels(subset(combined, combined$Variety != "IND" & combined$Cx=="ditransitive"))
head(vs.ind)
t.test(collexeme.list$IND.d[,3], vs.ind$log.pvalue)


# -- find examples of max collexeme strength -------------------
# -- .. verb - recipients----------------------------
head(collexeme.list$CAN.d)
head(as.data.frame(collexeme.list$CAN.d[,1]))
data[data$Variety=="CAN" & data$Resp=="do" & data$Verb=="give" & data$RecHeadPlain=="it", "FileID"] # 20
head(collexeme.list$CAN.p)
data[data$Variety=="CAN" & data$Resp=="pd" & data$Verb=="assign" & data$RecHeadLemma=="Association", "FileID"] # 3
head(as.data.frame(collexeme.list$CAN.p[,1]))

head(collexeme.list$GB.d)
head(as.data.frame(collexeme.list$GB.d[,1]))
data[data$Variety=="GB" & data$Resp=="do" & data$Verb=="cause" & data$RecHeadLemma=="defender", "FileID"] # 1
head(collexeme.list$GB.p)
data[data$Variety=="GB" & data$Resp=="pd" & data$Verb=="charge" & data$RecHeadLemma=="customer", "FileID"] #1
head(as.data.frame(collexeme.list$GB.p[,1]))

head(collexeme.list$HK.d)
head(as.data.frame(collexeme.list$HK.d[,1]))
data[data$Variety=="HK" & data$Resp=="do" & data$Verb=="drop" & data$RecHeadLemma=="Carolyn", "FileID"] #1
head(collexeme.list$HK.p)
data[data$Variety=="HK" & data$Resp=="pd" & data$Verb=="recommend" & data$RecHeadLemma=="other", "FileID"] #9
head(as.data.frame(collexeme.list$HK.p[,1]))

head(collexeme.list$IND.d)
head(as.data.frame(collexeme.list$IND.d[,1]))
data[data$Variety=="IND" & data$Resp=="do" & data$Verb=="promise" & data$RecHeadLemma=="Louise", "FileID"] # 2
head(collexeme.list$IND.p)
data[data$Variety=="IND" & data$Resp=="pd" & data$Verb=="pay" & data$RecHeadLemma=="company", "FileID"] # 4
head(as.data.frame(collexeme.list$IND.p[,1]))

head(collexeme.list$IRE.d)
head(as.data.frame(collexeme.list$IRE.d[,1]))
data[data$Variety=="IRE" & data$Resp=="do" & data$Verb=="give" & data$RecHeadLemma=="it", "FileID"] # 25
head(collexeme.list$IRE.p)
data[data$Variety=="IRE" & data$Resp=="pd" & data$Verb=="extend" & data$RecHeadLemma=="family", "FileID"] # 2
head(as.data.frame(collexeme.list$IRE.p[,1]))

head(collexeme.list$JA.d)
head(as.data.frame(collexeme.list$JA.d[,1]))
data[data$Variety=="JA" & data$Resp=="do" & data$Verb=="grant" & data$RecHeadLemma=="colony", "FileID"] # 1
head(collexeme.list$JA.p)
data[data$Variety=="JA" & data$Resp=="pd" & data$Verb=="hand" & data$RecHeadLemma=="Lodge", "FileID"] #2
head(as.data.frame(collexeme.list$JA.p[,1]))

head(collexeme.list$NZ.d)
head(as.data.frame(collexeme.list$NZ.d[,1]))
data[data$Variety=="NZ" & data$Resp=="do" & data$Verb=="give" & data$RecHeadLemma=="it", "FileID"] # 30
head(collexeme.list$NZ.p)
data[data$Variety=="NZ" & data$Resp=="pd" & data$Verb=="recommend" & data$RecHeadLemma=="friend", "FileID"] # 2
head(as.data.frame(collexeme.list$NZ.p[,1]))

head(collexeme.list$PHI.d)
head(as.data.frame(collexeme.list$PHI.d[,1]))
data[data$Variety=="PHI" & data$Resp=="do" & data$Verb=="show" & data$RecHeadLemma=="student", "FileID"] #4
head(collexeme.list$PHI.p)
data[data$Variety=="PHI" & data$Resp=="pd" & data$Verb=="deliver" & data$RecHeadLemma=="government", "FileID"] #1
head(as.data.frame(collexeme.list$PHI.p[,1]))

head(collexeme.list$SIN.d)
head(as.data.frame(collexeme.list$SIN.d[,1]))
data[data$Variety=="SIN" & data$Resp=="do" & data$Verb=="give" & data$RecHeadLemma=="it", "FileID"] # 28
head(collexeme.list$SIN.p)
data[data$Variety=="SIN" & data$Resp=="pd" & data$Verb=="pass" & data$RecHeadLemma=="Kalthom", "FileID"] # 1
head(as.data.frame(collexeme.list$SIN.p[,1]))


# -- .. verb - themes --------------
head(collexeme.list$CAN.d)
head(as.data.frame(collexeme.list$CAN.d[,1]))
data[data$Variety=="CAN" & data$Resp=="do" & data$Verb=="tell" & data$ThemeHeadLemma=="this", "FileID"] # 7
head(collexeme.list$CAN.p)
data[data$Variety=="CAN" & data$Resp=="pd" & data$Verb=="pay" & data$ThemeHeadLemma=="attention", "FileID"] #9
head(as.data.frame(collexeme.list$CAN.p[,1]))

head(collexeme.list$GB.d)
head(as.data.frame(collexeme.list$GB.d[,1]))
data[data$Variety=="GB" & data$Resp=="do" & data$Verb=="tell" & data$ThemeHeadLemma=="story", "FileID"] # 5
head(collexeme.list$GB.p)
data[data$Variety=="GB" & data$Resp=="pd" & data$Verb=="write" & data$ThemeHeadLemma=="letter", "FileID"] # 3
head(as.data.frame(collexeme.list$GB.p[,1]))

head(collexeme.list$HK.d)
head(as.data.frame(collexeme.list$HK.d[,1]))
data[data$Variety=="HK" & data$Resp=="do" & data$Verb=="send" & data$ThemeHeadLemma=="mail", "FileID"] # 10
head(collexeme.list$HK.p)
data[data$Variety=="HK" & data$Resp=="pd" & data$Verb=="pay" & data$ThemeHeadLemma=="attention", "FileID"] # 46
head(as.data.frame(collexeme.list$HK.p[,1]))

head(collexeme.list$IND.d)
head(as.data.frame(collexeme.list$IND.d[,1]))
data[data$Variety=="IND" & data$Resp=="do" & data$Verb=="drop" & data$ThemeHeadLemma=="line", "FileID"] # 3
head(collexeme.list$IND.p)
data[data$Variety=="IND" & data$Resp=="pd" & data$Verb=="pay" & data$ThemeHeadLemma=="attention", "FileID"] # 3
head(as.data.frame(collexeme.list$IND.p[,1]))

head(collexeme.list$IRE.d)
head(as.data.frame(collexeme.list$IRE.d[,1]))
data[data$Variety=="IRE" & data$Resp=="do" & data$Verb=="tell" & data$ThemeHeadLemma=="this", "FileID"] # 17
head(collexeme.list$IRE.p)
data[data$Variety=="IRE" & data$Resp=="pd" & data$Verb=="cause" & data$ThemeHeadLemma=="damage", "FileID"] # 2
head(as.data.frame(collexeme.list$IRE.p[,1]))

head(collexeme.list$JA.d)
head(as.data.frame(collexeme.list$JA.d[,1]))
data[data$Variety=="JA" & data$Resp=="do" & data$Verb=="tell" & data$ThemeHeadLemma=="something", "FileID"] # 25
head(collexeme.list$JA.p)
data[data$Variety=="JA" & data$Resp=="pd" & data$Verb=="pay" & data$ThemeHeadLemma=="attention", "FileID"] # 11
head(as.data.frame(collexeme.list$JA.p[,1]))

head(collexeme.list$NZ.d)
head(as.data.frame(collexeme.list$NZ.d[,1]))
data[data$Variety=="NZ" & data$Resp=="do" & data$Verb=="write" & data$ThemeHeadLemma=="letter", "FileID"] # 2
head(collexeme.list$NZ.p)
data[data$Variety=="NZ" & data$Resp=="pd" & data$Verb=="pay" & data$ThemeHeadLemma=="attention", "FileID"] # 15
head(as.data.frame(collexeme.list$NZ.p[,1]))

head(collexeme.list$PHI.d)
head(as.data.frame(collexeme.list$PHI.d[,1]))
data[data$Variety=="PHI" & data$Resp=="do" & data$Verb=="write" & data$ThemeHeadLemma=="letter", "FileID"] # 3
head(collexeme.list$PHI.p)
data[data$Variety=="PHI" & data$Resp=="pd" & data$Verb=="pay" & data$ThemeHeadLemma=="attention", "FileID"] #6
head(as.data.frame(collexeme.list$PHI.p[,1]))

head(collexeme.list$SIN.d)
head(as.data.frame(collexeme.list$SIN.d[,1]))
data[data$Variety=="SIN" & data$Resp=="do" & data$Verb=="tell" & data$ThemeHeadLemma=="story", "FileID"] # 10
head(collexeme.list$SIN.p)
data[data$Variety=="SIN" & data$Resp=="pd" & data$Verb=="pay" & data$ThemeHeadLemma=="attention", "FileID"] # 18
head(as.data.frame(collexeme.list$SIN.p[,1]))


# -- .. recipient - themes ------------
head(collexeme.list$CAN.d)
head(as.data.frame(collexeme.list$CAN.d[,1]))
data[data$Variety=="CAN" & data$Resp=="do" & data$RecHeadLemma=="application" & data$ThemeHeadLemma=="consideration", "FileID"] # 4
head(collexeme.list$CAN.p)
data[data$Variety=="CAN" & data$Resp=="pd" & data$RecHeadLemma=="Association" & data$ThemeHeadLemma=="remissions", "FileID"] # 1
head(as.data.frame(collexeme.list$CAN.p[,1]))

head(collexeme.list$GB.d)
head(as.data.frame(collexeme.list$GB.d[,1]))
data[data$Variety=="GB" & data$Resp=="do" & data$RecHeadLemma=="defender" & data$ThemeHeadLemma=="problem", "FileID"] # 1
head(collexeme.list$GB.p)
data[data$Variety=="GB" & data$Resp=="pd" & data$RecHeadLemma=="me" & data$ThemeHeadLemma=="it", "FileID"] # 3
head(as.data.frame(collexeme.list$GB.p[,1]))

head(collexeme.list$HK.d)
head(as.data.frame(collexeme.list$HK.d[,1]))
data[data$Variety=="HK" & data$Resp=="do" & data$RecHeadLemma=="her" & data$ThemeHeadLemma=="dollar", "FileID"] # 5
head(collexeme.list$HK.p)
data[data$Variety=="HK" & data$Resp=="pd" & data$RecHeadLemma=="other" & data$ThemeHeadLemma=="book", "FileID"] # 5
head(as.data.frame(collexeme.list$HK.p[,1]))

head(collexeme.list$IND.d)
head(as.data.frame(collexeme.list$IND.d[,1]))
data[data$Variety=="IND" & data$Resp=="do" & data$RecHeadLemma=="Louise" & data$ThemeHeadLemma=="happiness", "FileID"] # 2
head(collexeme.list$IND.p)
data[data$Variety=="IND" & data$Resp=="pd" & data$RecHeadLemma=="company" & data$ThemeHeadLemma=="amount", "FileID"] # 7
head(as.data.frame(collexeme.list$IND.p[,1]))

head(collexeme.list$IRE.d)
head(as.data.frame(collexeme.list$IRE.d[,1]))
data[data$Variety=="IRE" & data$Resp=="do" & data$RecHeadLemma=="it" & data$ThemeHeadLemma=="try", "FileID"] # 3
head(collexeme.list$IRE.p)
data[data$Variety=="IRE" & data$Resp=="pd" & data$RecHeadLemma=="node" & data$ThemeHeadLemma=="colour", "FileID"] # 1
head(as.data.frame(collexeme.list$IRE.p[,1]))

head(collexeme.list$JA.d)
head(as.data.frame(collexeme.list$JA.d[,1]))
data[data$Variety=="JA" & data$Resp=="do" & data$RecHeadLemma=="himself" & data$ThemeHeadLemma=="room", "FileID"] # 3
head(collexeme.list$JA.p)
data[data$Variety=="JA" & data$Resp=="pd" & data$RecHeadLemma=="police" & data$ThemeHeadLemma=="name", "FileID"] # 1
head(as.data.frame(collexeme.list$JA.p[,1]))

head(collexeme.list$NZ.d)
head(as.data.frame(collexeme.list$NZ.d[,1]))
data[data$Variety=="NZ" & data$Resp=="do" & data$RecHeadLemma=="God" & data$ThemeHeadLemma=="authority", "FileID"] # 1
head(collexeme.list$NZ.p)
data[data$Variety=="NZ" & data$Resp=="pd" & data$RecHeadLemma=="owner" & data$ThemeHeadLemma=="price", "FileID"] # 2
head(as.data.frame(collexeme.list$NZ.p[,1]))

head(collexeme.list$PHI.d)
head(as.data.frame(collexeme.list$PHI.d[,1]))
data[data$Variety=="PHI" & data$Resp=="do" & data$RecHeadLemma=="u" & data$ThemeHeadLemma=="room", "FileID"] # 1
head(collexeme.list$PHI.p)
data[data$Variety=="PHI" & data$Resp=="pd" & data$RecHeadLemma=="principle" & data$ThemeHeadLemma=="life", "FileID"] #1
head(as.data.frame(collexeme.list$PHI.p[,1]))

head(collexeme.list$SIN.d)
head(as.data.frame(collexeme.list$SIN.d[,1]))
data[data$Variety=="SIN" & data$Resp=="do" & data$RecHeadLemma=="Sato" & data$ThemeHeadLemma=="train", "FileID"] # 1
head(collexeme.list$SIN.p)
data[data$Variety=="SIN" & data$Resp=="pd" & data$RecHeadLemma=="winner" & data$ThemeHeadLemma=="medal", "FileID"] # 1
head(as.data.frame(collexeme.list$SIN.p[,1]))


# -- find examples ----------------------------------------------------------
# give it in IRE DOs
collexeme.list$IRE.d[collexeme.list$IRE.d$pair=="give it", ]
crdt <- data.list$IRE.d.data

for (j in 1:nrow(crdt)){
v <- crdt[j, "Verb"];  freq.v <- length(which(crdt$Verb == v))
r <- crdt[j, "RecHeadLemma"]; freq.r <- length(which(crdt$RecHeadLemma == r))
freq.vr <- length(which(crdt$Verb == v & crdt$RecHeadLemma == r))
}

# find number of occurrences in IRE.d.data
length(which(crdt$RecHeadLemma == "it"))
length(which(crdt$Verb == "give"))
length(which(crdt$Verb != "give" & crdt$RecHeadLemma!= "it"))
length(which(crdt$Verb == "give" & crdt$RecHeadLemma== "it"))


#------------------ DISTINCTIVE COLLEXEME ANALYSIS (STG) --------------------------
# -- Create files with separate columns for verb, theme, recipient --------------

#--- verb ----
coll.analysis_verb <- dats[,c("Variety", "Resp", "Verb")]
coll.analysis_verbCAN <- data.frame(droplevels(subset(coll.analysis_verb, Variety == "CAN")))
coll.analysis_verbCAN <- coll.analysis_verbCAN[,c("Resp", "Verb")]
write.table(data.frame(coll.analysis_verbCAN), "coll_verb/coll.analysis_verbCAN.txt", sep="\t", row.names=F, quote=F)

coll.analysis_verb <- dats[,c("Variety", "Resp", "Verb")]
coll.analysis_verbGB <- data.frame(droplevels(subset(coll.analysis_verb, Variety == "GB")))
coll.analysis_verbGB <- coll.analysis_verbGB[,c("Resp", "Verb")]
write.table(data.frame(coll.analysis_verbGB), "coll_verb/coll.analysis_verbGB.txt", sep="\t", row.names=F, quote=F)

coll.analysis_verb <- dats[,c("Variety", "Resp", "Verb")]
coll.analysis_verbHK <- data.frame(droplevels(subset(coll.analysis_verb, Variety == "HK")))
coll.analysis_verbHK <- coll.analysis_verbHK[,c("Resp", "Verb")]
write.table(data.frame(coll.analysis_verbHK), "coll_verb/coll.analysis_verbHK.txt", sep="\t", row.names=F, quote=F)

coll.analysis_verb <- dats[,c("Variety", "Resp", "Verb")]
coll.analysis_verbIND <- data.frame(droplevels(subset(coll.analysis_verb, Variety == "IND")))
coll.analysis_verbIND <- coll.analysis_verbIND[,c("Resp", "Verb")]
write.table(data.frame(coll.analysis_verbIND), "coll_verb/coll.analysis_verbIND.txt", sep="\t", row.names=F, quote=F)

coll.analysis_verb <- dats[,c("Variety", "Resp", "Verb")]
coll.analysis_verbIRE <- data.frame(droplevels(subset(coll.analysis_verb, Variety == "IRE")))
coll.analysis_verbIRE <- coll.analysis_verbIRE[,c("Resp", "Verb")]
write.table(data.frame(coll.analysis_verbIRE), "coll_verb/coll.analysis_verbIRE.txt", sep="\t", row.names=F, quote=F)

coll.analysis_verb <- dats[,c("Variety", "Resp", "Verb")]
coll.analysis_verbJA <- data.frame(droplevels(subset(coll.analysis_verb, Variety == "JA")))
coll.analysis_verbJA <- coll.analysis_verbJA[,c("Resp", "Verb")]
write.table(data.frame(coll.analysis_verbJA), "coll_verb/coll.analysis_verbJA.txt", sep="\t", row.names=F, quote=F)

coll.analysis_verb <- dats[,c("Variety", "Resp", "Verb")]
coll.analysis_verbNZ <- data.frame(droplevels(subset(coll.analysis_verb, Variety == "NZ")))
coll.analysis_verbNZ <- coll.analysis_verbNZ[,c("Resp", "Verb")]
write.table(data.frame(coll.analysis_verbNZ), "coll_verb/coll.analysis_verbNZ.txt", sep="\t", row.names=F, quote=F)

coll.analysis_verb <- dats[,c("Variety", "Resp", "Verb")]
coll.analysis_verbPHI <- data.frame(droplevels(subset(coll.analysis_verb, Variety == "PHI")))
coll.analysis_verbPHI <- coll.analysis_verbPHI[,c("Resp", "Verb")]
write.table(data.frame(coll.analysis_verbPHI), "coll_verb/coll.analysis_verbPHI.txt", sep="\t", row.names=F, quote=F)

coll.analysis_verb <- dats[,c("Variety", "Resp", "Verb")]
coll.analysis_verbSIN <- data.frame(droplevels(subset(coll.analysis_verb, Variety == "SIN")))
coll.analysis_verbSIN <- coll.analysis_verbSIN[,c("Resp", "Verb")]
write.table(data.frame(coll.analysis_verbSIN), "coll_verb/coll.analysis_verbSIN.txt", sep="\t", row.names=F, quote=F)

#---- theme ----
coll.analysis_theme <- dats[,c("Variety", "Resp", "ThemeHeadLemma")]
coll.analysis_themeCAN <- data.frame(droplevels(subset(coll.analysis_theme, Variety == "CAN")))
coll.analysis_themeCAN <- coll.analysis_themeCAN[,c("Resp", "ThemeHeadLemma")]
write.table(data.frame(coll.analysis_themeCAN), "coll_theme/coll.analysis_themeCAN.txt", sep="\t", row.names=F, quote=F)

coll.analysis_theme <- dats[,c("Variety", "Resp", "ThemeHeadLemma")]
coll.analysis_themeGB <- data.frame(droplevels(subset(coll.analysis_theme, Variety == "GB")))
coll.analysis_themeGB <- coll.analysis_themeGB[,c("Resp", "ThemeHeadLemma")]
write.table(data.frame(coll.analysis_themeGB), "coll_theme/coll.analysis_themeGB.txt", sep="\t", row.names=F, quote=F)

coll.analysis_theme <- dats[,c("Variety", "Resp", "ThemeHeadLemma")]
coll.analysis_themeHK <- data.frame(droplevels(subset(coll.analysis_theme, Variety == "HK")))
coll.analysis_themeHK <- coll.analysis_themeHK[,c("Resp", "ThemeHeadLemma")]
write.table(data.frame(coll.analysis_themeHK), "coll_theme/coll.analysis_themeHK.txt", sep="\t", row.names=F, quote=F)

coll.analysis_theme <- dats[,c("Variety", "Resp", "ThemeHeadLemma")]
coll.analysis_themeIND <- data.frame(droplevels(subset(coll.analysis_theme, Variety == "IND")))
coll.analysis_themeIND <- coll.analysis_themeIND[,c("Resp", "ThemeHeadLemma")]
write.table(data.frame(coll.analysis_themeIND), "coll_theme/coll.analysis_themeIND.txt", sep="\t", row.names=F, quote=F)

coll.analysis_theme <- dats[,c("Variety", "Resp", "ThemeHeadLemma")]
coll.analysis_themeIRE <- data.frame(droplevels(subset(coll.analysis_theme, Variety == "IRE")))
coll.analysis_themeIRE <- coll.analysis_themeIRE[,c("Resp", "ThemeHeadLemma")]
write.table(data.frame(coll.analysis_themeIRE), "coll_theme/coll.analysis_themeIRE.txt", sep="\t", row.names=F, quote=F)

coll.analysis_theme <- dats[,c("Variety", "Resp", "ThemeHeadLemma")]
coll.analysis_themeJA <- data.frame(droplevels(subset(coll.analysis_theme, Variety == "JA")))
coll.analysis_themeJA <- coll.analysis_themeJA[,c("Resp", "ThemeHeadLemma")]
write.table(data.frame(coll.analysis_themeJA), "coll_theme/coll.analysis_themeJA.txt", sep="\t", row.names=F, quote=F)

coll.analysis_theme <- dats[,c("Variety", "Resp", "ThemeHeadLemma")]
coll.analysis_themeNZ <- data.frame(droplevels(subset(coll.analysis_theme, Variety == "NZ")))
coll.analysis_themeNZ <- coll.analysis_themeNZ[,c("Resp", "ThemeHeadLemma")]
write.table(data.frame(coll.analysis_themeNZ), "coll_theme/coll.analysis_themeNZ.txt", sep="\t", row.names=F, quote=F)

coll.analysis_theme <- dats[,c("Variety", "Resp", "ThemeHeadLemma")]
coll.analysis_themePHI <- data.frame(droplevels(subset(coll.analysis_theme, Variety == "PHI")))
coll.analysis_themePHI <- coll.analysis_themePHI[,c("Resp", "ThemeHeadLemma")]
write.table(data.frame(coll.analysis_themePHI), "coll_theme/coll.analysis_themePHI.txt", sep="\t", row.names=F, quote=F)

coll.analysis_theme <- dats[,c("Variety", "Resp", "ThemeHeadLemma")]
coll.analysis_themeSIN <- data.frame(droplevels(subset(coll.analysis_theme, Variety == "SIN")))
coll.analysis_themeSIN <- coll.analysis_themeSIN[,c("Resp", "ThemeHeadLemma")]
write.table(data.frame(coll.analysis_themeSIN), "coll_theme/coll.analysis_themeSIN.txt", sep="\t", row.names=F, quote=F)

#--- recipient------
coll.analysis_rec <- dats[,c("Variety", "Resp", "RecHeadLemma")]
coll.analysis_recCAN <- data.frame(droplevels(subset(coll.analysis_rec, Variety == "CAN")))
coll.analysis_recCAN <- coll.analysis_recCAN[,c("Resp", "RecHeadLemma")]
write.table(data.frame(coll.analysis_recCAN), "coll_rec/coll.analysis_recCAN.txt", sep="\t", row.names=F, quote=F)

coll.analysis_rec <- dats[,c("Variety", "Resp", "RecHeadLemma")]
coll.analysis_recGB <- data.frame(droplevels(subset(coll.analysis_rec, Variety == "GB")))
coll.analysis_recGB <- coll.analysis_recGB[,c("Resp", "RecHeadLemma")]
write.table(data.frame(coll.analysis_recGB), "coll_rec/coll.analysis_recGB.txt", sep="\t", row.names=F, quote=F)

coll.analysis_rec <- dats[,c("Variety", "Resp", "RecHeadLemma")]
coll.analysis_recHK <- data.frame(droplevels(subset(coll.analysis_rec, Variety == "HK")))
coll.analysis_recHK <- coll.analysis_recHK[,c("Resp", "RecHeadLemma")]
write.table(data.frame(coll.analysis_recHK), "coll_rec/coll.analysis_recHK.txt", sep="\t", row.names=F, quote=F)

coll.analysis_rec <- dats[,c("Variety", "Resp", "RecHeadLemma")]
coll.analysis_recIND <- data.frame(droplevels(subset(coll.analysis_rec, Variety == "IND")))
coll.analysis_recIND <- coll.analysis_recIND[,c("Resp", "RecHeadLemma")]
write.table(data.frame(coll.analysis_recIND), "coll_rec/coll.analysis_recIND.txt", sep="\t", row.names=F, quote=F)

coll.analysis_rec <- dats[,c("Variety", "Resp", "RecHeadLemma")]
coll.analysis_recIRE <- data.frame(droplevels(subset(coll.analysis_rec, Variety == "IRE")))
coll.analysis_recIRE <- coll.analysis_recIRE[,c("Resp", "RecHeadLemma")]
write.table(data.frame(coll.analysis_recIRE), "coll_rec/coll.analysis_recIRE.txt", sep="\t", row.names=F, quote=F)

coll.analysis_rec <- dats[,c("Variety", "Resp", "RecHeadLemma")]
coll.analysis_recJA <- data.frame(droplevels(subset(coll.analysis_rec, Variety == "JA")))
coll.analysis_recJA <- coll.analysis_recJA[,c("Resp", "RecHeadLemma")]
write.table(data.frame(coll.analysis_recJA), "coll_rec/coll.analysis_recJA.txt", sep="\t", row.names=F, quote=F)

coll.analysis_rec <- dats[,c("Variety", "Resp", "RecHeadLemma")]
coll.analysis_recPHI <- data.frame(droplevels(subset(coll.analysis_rec, Variety == "PHI")))
coll.analysis_recPHI <- coll.analysis_recPHI[,c("Resp", "RecHeadLemma")]
write.table(data.frame(coll.analysis_recPHI), "coll_rec/coll.analysis_recPHI.txt", sep="\t", row.names=F, quote=F)

coll.analysis_rec <- dats[,c("Variety", "Resp", "RecHeadLemma")]
coll.analysis_recNZ <- data.frame(droplevels(subset(coll.analysis_rec, Variety == "NZ")))
coll.analysis_recNZ <- coll.analysis_recNZ[,c("Resp", "RecHeadLemma")]
write.table(data.frame(coll.analysis_recNZ), "coll_rec/coll.analysis_recNZ.txt", sep="\t", row.names=F, quote=F)

coll.analysis_rec <- dats[,c("Variety", "Resp", "RecHeadLemma")]
coll.analysis_recSIN <- data.frame(droplevels(subset(coll.analysis_rec, Variety == "SIN")))
coll.analysis_recSIN <- coll.analysis_recSIN[,c("Resp", "RecHeadLemma")]
write.table(data.frame(coll.analysis_recSIN), "coll_rec/coll.analysis_recSIN.txt", sep="\t", row.names=F, quote=F)

#----------verb_theme------------------
coll.analysis_vt <- dats[,c("Variety", "Resp", "VerbThemeLemma")]
coll.analysis_vtCAN <- data.frame(droplevels(subset(coll.analysis_vt, Variety == "CAN")))
coll.analysis_vtCAN <- coll.analysis_vtCAN[,c("Resp", "VerbThemeLemma")]
write.table(data.frame(coll.analysis_vtCAN), "coll.analysis_vtCAN.txt", sep="\t", row.names=F, quote=F)

coll.analysis_vt <- dats[,c("Variety", "Resp", "VerbThemeLemma")]
coll.analysis_vtGB <- data.frame(droplevels(subset(coll.analysis_vt, Variety == "GB")))
coll.analysis_vtGB <- coll.analysis_vtGB[,c("Resp", "VerbThemeLemma")]
write.table(data.frame(coll.analysis_vtGB), "coll.analysis_vtGB.txt", sep="\t", row.names=F, quote=F)

coll.analysis_vt <- dats[,c("Variety", "Resp", "VerbThemeLemma")]
coll.analysis_vtHK <- data.frame(droplevels(subset(coll.analysis_vt, Variety == "HK")))
coll.analysis_vtHK <- coll.analysis_vtHK[,c("Resp", "VerbThemeLemma")]
write.table(data.frame(coll.analysis_vtHK), "coll.analysis_vtHK.txt", sep="\t", row.names=F, quote=F)

coll.analysis_vt <- dats[,c("Variety", "Resp", "VerbThemeLemma")]
coll.analysis_vtIND <- data.frame(droplevels(subset(coll.analysis_vt, Variety == "IND")))
coll.analysis_vtIND <- coll.analysis_vtIND[,c("Resp", "VerbThemeLemma")]
write.table(data.frame(coll.analysis_vtIND), "coll.analysis_vtIND.txt", sep="\t", row.names=F, quote=F)

coll.analysis_vt <- dats[,c("Variety", "Resp", "VerbThemeLemma")]
coll.analysis_vtIRE <- data.frame(droplevels(subset(coll.analysis_vt, Variety == "IRE")))
coll.analysis_vtIRE <- coll.analysis_vtIRE[,c("Resp", "VerbThemeLemma")]
write.table(data.frame(coll.analysis_vtIRE), "coll.analysis_vtIRE.txt", sep="\t", row.names=F, quote=F)

coll.analysis_vt <- dats[,c("Variety", "Resp", "VerbThemeLemma")]
coll.analysis_vtJA <- data.frame(droplevels(subset(coll.analysis_vt, Variety == "JA")))
coll.analysis_vtJA <- coll.analysis_vtJA[,c("Resp", "VerbThemeLemma")]
write.table(data.frame(coll.analysis_vtJA), "coll.analysis_vtJA.txt", sep="\t", row.names=F, quote=F)

coll.analysis_vt <- dats[,c("Variety", "Resp", "VerbThemeLemma")]
coll.analysis_vtNZ <- data.frame(droplevels(subset(coll.analysis_vt, Variety == "NZ")))
coll.analysis_vtNZ <- coll.analysis_vtNZ[,c("Resp", "VerbThemeLemma")]
write.table(data.frame(coll.analysis_vtNZ), "coll.analysis_vtNZ.txt", sep="\t", row.names=F, quote=F)

coll.analysis_vt <- dats[,c("Variety", "Resp", "VerbThemeLemma")]
coll.analysis_vtPHI <- data.frame(droplevels(subset(coll.analysis_vt, Variety == "PHI")))
coll.analysis_vtPHI <- coll.analysis_vtPHI[,c("Resp", "VerbThemeLemma")]
write.table(data.frame(coll.analysis_vtPHI), "coll.analysis_vtPHI.txt", sep="\t", row.names=F, quote=F)

coll.analysis_vt <- dats[,c("Variety", "Resp", "VerbThemeLemma")]
coll.analysis_vtSIN <- data.frame(droplevels(subset(coll.analysis_vt, Variety == "SIN")))
coll.analysis_vtSIN <- coll.analysis_vtSIN[,c("Resp", "VerbThemeLemma")]
write.table(data.frame(coll.analysis_vtSIN), "coll.analysis_vtSIN.txt", sep="\t", row.names=F, quote=F)

#----------verb_sense-type------------------

coll.analysis_vs <- dats[,c("Variety", "Resp", "VerbSense")]
coll.analysis_vsCAN <- data.frame(droplevels(subset(coll.analysis_vs, Variety == "CAN")))
coll.analysis_vsCAN <- coll.analysis_vsCAN[,c("Resp", "VerbSense")]
write.table(data.frame(coll.analysis_vsCAN), "coll.analysis_vsCAN.txt", sep="\t", row.names=F, quote=F)

coll.analysis_vs <- dats[,c("Variety", "Resp", "VerbSense")]
coll.analysis_vsGB <- data.frame(droplevels(subset(coll.analysis_vs, Variety == "GB")))
coll.analysis_vsGB <- coll.analysis_vsGB[,c("Resp", "VerbSense")]
write.table(data.frame(coll.analysis_vsGB), "coll.analysis_vsGB.txt", sep="\t", row.names=F, quote=F)

coll.analysis_vs <- dats[,c("Variety", "Resp", "VerbSense")]
coll.analysis_vsHK <- data.frame(droplevels(subset(coll.analysis_vs, Variety == "HK")))
coll.analysis_vsHK <- coll.analysis_vsHK[,c("Resp", "VerbSense")]
write.table(data.frame(coll.analysis_vsHK), "coll.analysis_vsHK.txt", sep="\t", row.names=F, quote=F)

coll.analysis_vs <- dats[,c("Variety", "Resp", "VerbSense")]
coll.analysis_vsIND <- data.frame(droplevels(subset(coll.analysis_vs, Variety == "IND")))
coll.analysis_vsIND <- coll.analysis_vsIND[,c("Resp", "VerbSense")]
write.table(data.frame(coll.analysis_vsIND), "coll.analysis_vsIND.txt", sep="\t", row.names=F, quote=F)

coll.analysis_vs <- dats[,c("Variety", "Resp", "VerbSense")]
coll.analysis_vsIRE <- data.frame(droplevels(subset(coll.analysis_vs, Variety == "IRE")))
coll.analysis_vsIRE <- coll.analysis_vsIRE[,c("Resp", "VerbSense")]
write.table(data.frame(coll.analysis_vsIRE), "coll.analysis_vsIRE.txt", sep="\t", row.names=F, quote=F)

coll.analysis_vs <- dats[,c("Variety", "Resp", "VerbSense")]
coll.analysis_vsJA <- data.frame(droplevels(subset(coll.analysis_vs, Variety == "JA")))
coll.analysis_vsJA <- coll.analysis_vsJA[,c("Resp", "VerbSense")]
write.table(data.frame(coll.analysis_vsJA), "coll.analysis_vsJA.txt", sep="\t", row.names=F, quote=F)

coll.analysis_vs <- dats[,c("Variety", "Resp", "VerbSense")]
coll.analysis_vsNZ <- data.frame(droplevels(subset(coll.analysis_vs, Variety == "NZ")))
coll.analysis_vsNZ <- coll.analysis_vsNZ[,c("Resp", "VerbSense")]
write.table(data.frame(coll.analysis_vsNZ), "coll.analysis_vsNZ.txt", sep="\t", row.names=F, quote=F)

coll.analysis_vs <- dats[,c("Variety", "Resp", "VerbSense")]
coll.analysis_vsPHI <- data.frame(droplevels(subset(coll.analysis_vs, Variety == "PHI")))
coll.analysis_vsPHI <- coll.analysis_vsPHI[,c("Resp", "VerbSense")]
write.table(data.frame(coll.analysis_vsPHI), "coll.analysis_vsPHI.txt", sep="\t", row.names=F, quote=F)

coll.analysis_vs <- dats[,c("Variety", "Resp", "VerbSense")]
coll.analysis_vsSIN <- data.frame(droplevels(subset(coll.analysis_vs, Variety == "SIN")))
coll.analysis_vsSIN <- coll.analysis_vsSIN[,c("Resp", "VerbSense")]
write.table(data.frame(coll.analysis_vsSIN), "coll.analysis_vsSIN.txt", sep="\t", row.names=F, quote=F)

# -- Analysis: external script ------------------------------
rm(list=ls(all=TRUE))
source("coll.analysis_mpfr.r") # script by Stefan Gries

# In order to go through the analysis, type the below and the programm will guide you through
coll.analysis()

#--------- RESULTS verb-------------
#-------...mean coll.strength verb-----------
# Look at mean coll.strength of constituents in all varieties by preferred variant
verbs <- read.delim("/Users/melanierothlisberger/Dropbox/Universitaet/04_Statistics/statistical analyses/20170209_THESIS_full_model/Lexical/coll_verb/verb_all_results2.csv", strip.white = T, header=T)


# separate by variety and variant
verbs = verbs[, c("variety", "words", "pref.occur", "coll.strength")]
vars <- levels(verbs$variety)
data.list3 = list()
for (i in seq(2,18,2)){
  v <- vars[(i/2)]
  sub.data <- droplevels(subset(verbs, variety==v)) # drop unused levels
  sub.data.d <- droplevels(subset(sub.data, pref.occur=="ditransitive"))
  sub.data.p <- droplevels(subset(sub.data, pref.occur=="prepositional"))
  data.list3[[i - 1]] <- sub.data.d; names(data.list3)[i - 1] <- paste(v, "d.data", sep = '.')
  data.list3[[i]] <- sub.data.p; names(data.list3)[i] <- paste(v, "p.data", sep = '.')
  rm(list = c("sub.data.d", "sub.data.p")) # remove temporary stuff
}


# barplot of mean coll strength
cols <- c(ditransitive="dimgray",prepositional="gainsboro")

collstr.mean.verb <- summarySE(verbs, measurevar = "coll.strength", groupvars = c("variety", "pref.occur"))

# Make barplot: 
# mean.collstr.verb: verb collostructions

collstr.mean.verb$variety <- factor(collstr.mean.verb$variety, levels=c("CAN", "GB", "IRE", "NZ","JA", "SIN", "HK", "IND", "PHI"))

(mean.collstr.verb <- ggplot(collstr.mean.verb, aes(variety, coll.strength, fill=pref.occur)) +
    geom_bar(stat="identity", position=position_dodge(), color="black", width=.8) +
    geom_text(aes(y=coll.strength+0.4, label=round(coll.strength, 2)), position=position_dodge(0.8), size=5. ) +
    geom_errorbar(aes(ymin=coll.strength-ci, ymax=coll.strength+ci), colour="black", width=.1, position=position_dodge(width=0.8), linetype=1) +
    scale_fill_manual(values=cols, name = "Variant")+
    labs(title="verb collostructions", y="mean collstr. strength", x="") +
    theme(axis.title.y = element_text(size=20)) +
    theme(axis.text.x = element_text(size=20), plot.title = element_text(size=20, face="bold", hjust = 0.5)) +
    guides(fill=guide_legend(ncol=2,title.position = "top", title.hjust = 0.5))+
    theme(legend.position = "bottom",legend.title = element_text(size=20,face = "bold"), legend.text = element_text(size=20)))

pdf("Rplot_lexical_collostruction_mean_verb.pdf", width = 12, height = 7)
mean.collstr.verb
dev.off()




#-- ..Top 7 verbs per variety and variant----
collstr.v.do <- as.data.frame(data.list3$CAN.d.data[1:3,]) %>% droplevels
collstr.v.do <- rbind(collstr.v.do, as.data.frame(data.list3$GB.d.data[1:3,]), as.data.frame(data.list3$HK.d.data[1:3,]), as.data.frame(data.list3$IND.d.data[1:3,]), as.data.frame(data.list3$IRE.d.data[1:3,]), as.data.frame(data.list3$JA.d.data[1:3,]), as.data.frame(data.list3$NZ.d.data[1:3,]), as.data.frame(data.list3$PHI.d.data[1:3,]), as.data.frame(data.list3$SIN.d.data[1:3,])) %>% droplevels


row.names(collstr.v.do) <- NULL
collstr.v.do$variety <- factor(collstr.v.do$variety, levels=c("CAN", "GB", "IRE", "NZ","JA", "SIN", "HK", "IND", "PHI"))

# reorder to plot from highest to lowest collstrength
collstr.v.do$words <- reorder(collstr.v.do$words, -collstr.v.do$coll.strength)

unique(collstr.v.do$words)

(ex.collstr.verb <- ggplot(collstr.v.do, aes(variety, coll.strength, fill=words)) +
  geom_bar(stat="identity", position=position_dodge(), color="black", width=.8) +
  geom_text(aes(label=round(coll.strength,1)), vjust=-0.4, position=position_dodge(0.8), size=4.5) +
  scale_fill_brewer(palette = "BrBG", direction=-1, name="verb") +
  labs(title="ditransitive dative", y="collostructional strength", x="") +
    theme(axis.title.x = element_text(size=12), axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=16), plot.title = element_text(size=16, face="bold", hjust = 0.5)) +
    guides(fill=guide_legend(ncol=3,title.position = "top", title.hjust = 0.5))+
    theme(legend.position = "bottom",legend.title = element_text(size=16,face = "bold"), legend.text = element_text(size=16)))


pdf("Rplot_lexical_collostruction_examples_verb_do.pdf", width = 10, height = 6)
ex.collstr.verb
dev.off()

# Summary of top 10 verbs and in which varieties they occur
addmargins(xtabs(~ words + variety, data=collstr.v.do))  
nrow(xtabs(~words+variety, data=collstr.v.do))


# PREPOSITIONAL DATIVE
collstr.v.pd <- as.data.frame(data.list3$CAN.p.data[1:3,]) %>% droplevels
collstr.v.pd <- rbind(collstr.v.pd, as.data.frame(data.list3$GB.p.data[1:3,]), as.data.frame(data.list3$HK.p.data[1:3,]), as.data.frame(data.list3$IND.p.data[1:3,]), as.data.frame(data.list3$IRE.p.data[1:3,]), as.data.frame(data.list3$JA.p.data[1:3,]), as.data.frame(data.list3$NZ.p.data[1:3,]), as.data.frame(data.list3$PHI.p.data[1:3,]), as.data.frame(data.list3$SIN.p.data[1:3,])) %>% droplevels


row.names(collstr.v.pd) <- NULL
collstr.v.pd$variety <- factor(collstr.v.pd$variety, levels=c("CAN", "GB", "IRE", "NZ","JA", "SIN", "HK", "IND", "PHI"))

# reorder to plot from highest to lowest collstrength
collstr.v.pd$words <- reorder(collstr.v.pd$words, -collstr.v.pd$coll.strength)

unique(collstr.v.pd$words) # 10 verbs

(ex.collstr.verb.pd <- ggplot(collstr.v.pd, aes(variety, coll.strength, fill=words)) +
  geom_bar(stat="identity", position=position_dodge(), color="black", width=.8) +
  geom_text(aes(label=round(coll.strength,1)), vjust=-0.4, position=position_dodge(0.8), size=4.5) +
  scale_fill_brewer(palette = "BrBG", direction=-1, name="verb") +
    labs(title="prepositional dative", y="collostructional strength", x="") +
    theme(axis.title.x = element_text(size=12), axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=16), plot.title = element_text(size=16, face="bold", hjust = 0.5)) +
    guides(fill=guide_legend(ncol=5,title.position = "top", title.hjust = 0.5))+
    theme(legend.position = "bottom",legend.title = element_text(size=16,face = "bold"), legend.text = element_text(size=16)))


pdf("Rplot_lexical_collostruction_examples_verb_pd.pdf", width = 12, height = 7)
ex.collstr.verb.pd
dev.off()

# Summary of top 10 verbs and in which varieties they occur
addmargins(xtabs(~ words + variety, data=collstr.th.do))  
nrow(xtabs(~words+variety, data=collstr.th.do))

pdf("Rplot_lexical_collostruction_examples_verb_all.pdf", width = 12, height = 15)
grid.arrange(arrangeGrob(ex.collstr.verb, ex.collstr.verb.pd))
dev.off()


#---------- RESULTS Recipient-------------------
#---------...mean coll.strength recipients----------
# Look at mean coll.strength of constituents in all varieties by preferred variant
recipients <- read.delim("/Users/melanierothlisberger/Dropbox/Universitaet/04_Statistics/statistical analyses/20170209_THESIS_full_model/Lexical/coll_rec/rec_all_results2.csv", strip.white = T, header=T)

# separate by variety and variant
recipients = recipients[, c("variety", "words", "pref.occur", "coll.strength")]
vars <- levels(recipients$variety)
data.lst = list()
for (i in seq(2,18,2)){
  v <- vars[(i/2)]
  sub.data <- droplevels(subset(recipients, variety==v)) # drop unused levels
  sub.data.d <- droplevels(subset(sub.data, pref.occur=="ditransitive"))
  sub.data.p <- droplevels(subset(sub.data, pref.occur=="prepositional"))
  data.lst[[i - 1]] <- sub.data.d; names(data.lst)[i - 1] <- paste(v, "d.data", sep = '.')
  data.lst[[i]] <- sub.data.p; names(data.lst)[i] <- paste(v, "p.data", sep = '.')
  rm(list = c("sub.data.d", "sub.data.p")) # remove temporary stuff
}

# matrix table
vars <- levels(recipients$variety)


# barplot of mean coll strength (df.list)
cols <- c(ditransitive="dimgray",prepositional="gainsboro")

collstr.mean.rec <- summarySE(recipients, measurevar = "coll.strength", groupvars = c("variety", "pref.occur"))

# change order of factors to be displayed for variety:
collstr.mean.rec$variety <- factor(collstr.mean.rec$variety, levels=c("CAN", "GB", "IRE", "NZE","JA", "SIN", "HK", "IND", "PHI"))

# Make barplot: 

(mean.collstr.recip <- ggplot(collstr.mean.rec, aes(variety,coll.strength, fill=pref.occur)) +
  geom_bar(stat="identity", position=position_dodge(), color="black", width=.8) +
    geom_text(aes(y=0.08, label=round(coll.strength, 2)), position=position_dodge(0.8), size=5. ) +
    geom_errorbar(aes(ymin=coll.strength-ci, ymax=coll.strength+ci), colour="black", width=.1, position=position_dodge(width=0.8), linetype=1) +
    scale_fill_manual(values=cols, name = "Variant")+
    labs(title="recipient collostructions", y="mean collstr. strength", x="") +
    theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=20)) +
    theme(axis.text.x = element_text(size=20), plot.title = element_text(size=20, face="bold", hjust = 0.5)) +
    guides(fill=guide_legend(ncol=2,title.position = "top", title.hjust = 0.5))+
    theme(legend.position = "bottom",legend.title = element_text(size=20,face = "bold"), legend.text = element_text(size=20)))

pdf("Rplot_lexical_collostruction_mean_recipient.pdf", width = 12, height = 7)
mean.collstr.recip
dev.off()


#-- ..Top 7 recipients (=7 pronouns) per variety and variant----
collstr.rec.do <- as.data.frame(data.lst$CAN.d.data[1:7,])
collstr.rec.do <- rbind(collstr.rec.do, as.data.frame(data.lst$GB.d.data[1:7,]))
collstr.rec.do <- rbind(collstr.rec.do, as.data.frame(data.lst$HK.d.data[1:7,]))
collstr.rec.do <- rbind(collstr.rec.do, as.data.frame(data.lst$IND.d.data[1:7,]))
collstr.rec.do <- rbind(collstr.rec.do, as.data.frame(data.lst$IRE.d.data[1:7,]))
collstr.rec.do <- rbind(collstr.rec.do, as.data.frame(data.lst$JA.d.data[1:7,]))
collstr.rec.do <- rbind(collstr.rec.do, as.data.frame(data.lst$NZE.d.data[1:7,]))
collstr.rec.do <- rbind(collstr.rec.do, as.data.frame(data.lst$PHI.d.data[1:7,]))
collstr.rec.do <- rbind(collstr.rec.do, as.data.frame(data.lst$SIN.d.data[1:7,]))
row.names(collstr.rec.do) <- NULL
collstr.rec.do$variety <- factor(collstr.rec.do$variety, levels=c("CAN", "GB", "IRE", "NZE","JA", "SIN", "HK", "IND", "PHI"))

# reorder to plot from highest to lowest collstrength
collstr.rec.do$words <- reorder(collstr.rec.do$words, -collstr.rec.do$coll.strength)

(ex.collstr.rec.do <- ggplot(collstr.rec.do, aes(variety, coll.strength, fill=words)) +
  geom_bar(stat="identity", position=position_dodge(), color="black", width=.8) +
  geom_text(aes(label=round(coll.strength,0)), vjust=-0.4, position=position_dodge(0.8), size=4) +
    scale_fill_brewer(palette = "BrBG", direction=-1, name="recipient") +
    labs(title="ditransitive dative", y="collostructional strength", x="") +
    theme(axis.title.x = element_text(size=12), axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=16), plot.title = element_text(size=16, face="bold", hjust = 0.5)) +
    guides(fill=guide_legend(ncol=5,title.position = "top", title.hjust = 0.5))+
    theme(legend.position = "bottom",legend.title = element_text(size=16,face = "bold"), legend.text = element_text(size=16)))

pdf("Rplot_lexical_collostruction_examples_rec_do.pdf", width = 12, height = 7)
ex.collstr.rec.do
dev.off()


# PREPOSITIONAL DATIVE
as.data.frame(data.lst$CAN.p.data[1:7,])
as.data.frame(data.lst$GB.p.data[1:7,])
as.data.frame(data.lst$IRE.p.data[1:7,])
as.data.frame(data.lst$NZE.p.data[1:7,])
as.data.frame(data.lst$JA.p.data[1:7,])
as.data.frame(data.lst$SIN.p.data[1:7,])
as.data.frame(data.lst$HK.p.data[1:7,])
as.data.frame(data.lst$IND.p.data[1:7,])
as.data.frame(data.lst$PHI.p.data[1:7,])


#--------------RESULTS Themes-----------------
#---------...mean coll.strength themes----------
# Look at mean coll.strength of constituents in all varieties by preferred variant
themes <- read.delim("/Users/melanierothlisberger/Dropbox/Universitaet/04_Statistics/statistical analyses/20170209_THESIS_full_model/Lexical/coll_theme/th_all_results2.csv", strip.white = T, header=T)

# separate by variety and variant
themes = themes[, c("variety", "words", "pref.occur", "coll.strength")]
vars <- levels(themes$variety)
data.list2 = list()
for (i in seq(2,18,2)){
  v <- vars[(i/2)]
  sub.data <- droplevels(subset(themes, variety==v)) # drop unused levels
  sub.data.d <- droplevels(subset(sub.data, pref.occur=="ditransitive"))
  sub.data.p <- droplevels(subset(sub.data, pref.occur=="prepositional"))
  data.list2[[i - 1]] <- sub.data.d; names(data.list2)[i - 1] <- paste(v, "d.data", sep = '.')
  data.list2[[i]] <- sub.data.p; names(data.list2)[i] <- paste(v, "p.data", sep = '.')
  rm(list = c("sub.data.d", "sub.data.p")) # remove temporary stuff
}

collstr.mean.th <- summarySE(themes, measurevar = "coll.strength", groupvars = c("variety", "pref.occur"))

# change order of factors to be displayed for variety:
collstr.mean.th$variety <- factor(collstr.mean.th$variety, levels=c("CAN", "GB", "IRE", "NZ","JA", "SIN", "HK", "IND", "PHI"))


# barplot of mean coll strength (df.list2)
cols <- c(ditransitive="dimgray",prepositional="gainsboro")

(mean.collstr.th <- ggplot(collstr.mean.th, aes(variety,coll.strength, fill=pref.occur)) +
  geom_bar(stat="identity", position=position_dodge(), color="black", width=.8) +
    geom_text(aes(y=0.03, label=round(coll.strength, 2)), position=position_dodge(0.8), size=5) +
  geom_errorbar(aes(ymin=coll.strength-ci, ymax=coll.strength+ci), colour="black", width=.1, position=position_dodge(width=0.8), linetype=1) +
    scale_fill_manual(values=cols, name = "Variant")+
    labs(title="theme collostructions", y="mean collstr. strength", x="") +
    theme(axis.title.y = element_text(size=20)) +
    theme(axis.text.x = element_text(size=20), plot.title = element_text(size=20, face="bold", hjust = 0.5)) +
    guides(fill=guide_legend(ncol=2,title.position = "top", title.hjust = 0.5))+
    theme(legend.position = "bottom",legend.title = element_text(size=20,face = "bold"), legend.text = element_text(size=20)))

pdf("Rplot_lexical_collostruction_mean_theme.pdf", width = 12, height = 7)
mean.collstr.th
dev.off()

# -- ...t-test to compare means ------------------------------------------
head(themes)

length(as.vector(themes[themes$variety=="CAN" & themes$pref.occur=="ditransitive", "coll.strength"])) # CanE, ditrans, number of tokens in coll.strength

# is the data normally distributed?
plot(density(themes[themes$variety=="CAN" & themes$pref.occur=="ditransitive", "coll.strength"]))
plot(density(themes[themes$variety=="CAN" & themes$pref.occur=="prepositional", "coll.strength"]))
plot(density(themes[themes$variety=="GB" & themes$pref.occur=="ditransitive", "coll.strength"]))
plot(density(themes[themes$variety=="GB" & themes$pref.occur=="prepositional", "coll.strength"]))
plot(density(themes[themes$variety=="HK" & themes$pref.occur=="ditransitive", "coll.strength"]))
plot(density(themes[themes$variety=="HK" & themes$pref.occur=="prepositional", "coll.strength"]))
plot(density(themes[themes$variety=="IND" & themes$pref.occur=="ditransitive", "coll.strength"]))
plot(density(themes[themes$variety=="IND" & themes$pref.occur=="prepositional", "coll.strength"]))
plot(density(themes[themes$variety=="IRE" & themes$pref.occur=="ditransitive", "coll.strength"]))
plot(density(themes[themes$variety=="IRE" & themes$pref.occur=="prepositional", "coll.strength"]))
plot(density(themes[themes$variety=="JA" & themes$pref.occur=="ditransitive", "coll.strength"]))
plot(density(themes[themes$variety=="JA" & themes$pref.occur=="prepositional", "coll.strength"]))
plot(density(themes[themes$variety=="NZ" & themes$pref.occur=="ditransitive", "coll.strength"]))
plot(density(themes[themes$variety=="NZ" & themes$pref.occur=="prepositional", "coll.strength"]))
plot(density(themes[themes$variety=="PHI" & themes$pref.occur=="ditransitive", "coll.strength"]))
plot(density(themes[themes$variety=="PHI" & themes$pref.occur=="prepositional", "coll.strength"]))
plot(density(themes[themes$variety=="SIN" & themes$pref.occur=="ditransitive", "coll.strength"]))
plot(density(themes[themes$variety=="SIN" & themes$pref.occur=="prepositional", "coll.strength"]))

# very long tails, binomial indications, but lots of data points

qqnorm(themes[themes$variety=="SIN" & themes$pref.occur=="prepositional", "coll.strength"])
qqline(themes[themes$variety=="SIN" & themes$pref.occur=="prepositional", "coll.strength"])

shapiro.test(collexeme.list$CAN.d[,3])
shapiro.test(collexeme.list$CAN.p[,3])
shapiro.test(collexeme.list$GB.d[,3])
shapiro.test(collexeme.list$GB.p[,3])


# run t-tests on all separately by variety, correct for Bonferroni
test.can <- t.test(themes[themes$variety=="CAN" & themes$pref.occur=="ditransitive", "coll.strength"], themes[themes$variety=="CAN" & themes$pref.occur=="prepositional", "coll.strength"])
pCAN <- test.can$p.value
test.gb <- t.test(themes[themes$variety=="GB" & themes$pref.occur=="ditransitive", "coll.strength"], themes[themes$variety=="GB" & themes$pref.occur=="prepositional", "coll.strength"])
pGB <- test.gb$p.value
test.hk <- t.test(themes[themes$variety=="HK" & themes$pref.occur=="ditransitive", "coll.strength"], themes[themes$variety=="HK" & themes$pref.occur=="prepositional", "coll.strength"])
pHK <- test.hk$p.value
test.ind <- t.test(themes[themes$variety=="IND" & themes$pref.occur=="ditransitive", "coll.strength"], themes[themes$variety=="IND" & themes$pref.occur=="prepositional", "coll.strength"])
pIND <- test.ind$p.value
test.ire <- t.test(themes[themes$variety=="IRE" & themes$pref.occur=="ditransitive", "coll.strength"], themes[themes$variety=="IRE" & themes$pref.occur=="prepositional", "coll.strength"])
pIRE <- test.ire$p.value
test.ja <- t.test(themes[themes$variety=="JA" & themes$pref.occur=="ditransitive", "coll.strength"], themes[themes$variety=="JA" & themes$pref.occur=="prepositional", "coll.strength"])
pJA <- test.ja$p.value
test.nz <- t.test(themes[themes$variety=="NZ" & themes$pref.occur=="ditransitive", "coll.strength"], themes[themes$variety=="NZ" & themes$pref.occur=="prepositional", "coll.strength"])
pNZ <- test.nz$p.value
test.phi <- t.test(themes[themes$variety=="PHI" & themes$pref.occur=="ditransitive", "coll.strength"], themes[themes$variety=="PHI" & themes$pref.occur=="prepositional", "coll.strength"])
pPHI <- test.phi$p.value
test.sin <- t.test(themes[themes$variety=="SIN" & themes$pref.occur=="ditransitive", "coll.strength"], themes[themes$variety=="SIN" & themes$pref.occur=="prepositional", "coll.strength"])
pSIN <- test.sin$p.value

p.adjust(c(pCAN, pGB, pHK, pIND, pIRE, pJA, pNZ, pPHI, pSIN), method = "bonferroni")


# compare mean ditrans IND vs global average
# get global subset -IND
th.not.ind <- droplevels(subset(themes, themes$variety != "IND" & themes$pref.occur=="ditransitive"))
head(th.not.ind)
t.test(themes[themes$variety=="IND" & themes$pref.occur=="ditransitive", "coll.strength"], th.not.ind$coll.strength)



#-- ..Top 3 themes per variety and variant----
collstr.th.do <- as.data.frame(data.list2$CAN.d.data[1:3,]) %>% droplevels
collstr.th.do <- rbind(collstr.th.do, as.data.frame(data.list2$GB.d.data[1:3,]), as.data.frame(data.list2$HK.d.data[1:3,]), as.data.frame(data.list2$IND.d.data[1:3,]), as.data.frame(data.list2$IRE.d.data[1:3,]), as.data.frame(data.list2$JA.d.data[1:3,]), as.data.frame(data.list2$NZ.d.data[1:3,]), as.data.frame(data.list2$PHI.d.data[1:3,]), as.data.frame(data.list2$SIN.d.data[1:3,])) %>% droplevels


row.names(collstr.th.do) <- NULL
collstr.th.do$variety <- factor(collstr.th.do$variety, levels=c("CAN", "GB", "IRE", "NZ","JA", "SIN", "HK", "IND", "PHI"))

# reorder to plot from highest to lowest collstrength
collstr.th.do$words <- reorder(collstr.th.do$words, -collstr.th.do$coll.strength)

unique(collstr.th.do$words)

(ex.collstr.th.do <- ggplot(collstr.th.do, aes(variety, coll.strength, fill=words)) +
  geom_bar(stat="identity", position=position_dodge(), color="black", width=.8) +
  geom_text(aes(label=round(coll.strength,1)), vjust=-0.4, position=position_dodge(0.8), size=4) +
    scale_fill_brewer(palette = "BrBG", direction=-1, name="theme") +
    labs(title="ditransitive dative", y="collostructional strength", x="") +
    theme(axis.title.x = element_text(size=12), axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=16), plot.title = element_text(size=16, face="bold", hjust = 0.5)) +
    guides(fill=guide_legend(ncol=3,title.position = "top", title.hjust = 0.5))+
    theme(legend.position = "bottom",legend.title = element_text(size=16,face = "bold"), legend.text = element_text(size=16)))

pdf("Rplot_lexical_collostruction_examples_th_do.pdf", width = 12, height = 7)
ex.collstr.th.do
dev.off()


# Summary of top 10 themes and in which varieties they occur
addmargins(xtabs(~ words + variety, data=collstr.th.do))  
nrow(xtabs(~words+variety, data=collstr.th.do))


# PREPOSITIONAL DATIVE
as.data.frame(data.list2$CAN.p.data[1:7,])
as.data.frame(data.list2$GB.p.data[1:7,])
as.data.frame(data.list2$IRE.p.data[1:7,])
as.data.frame(data.list2$NZ.p.data[1:7,])
as.data.frame(data.list2$JA.p.data[1:7,])
as.data.frame(data.list2$SIN.p.data[1:7,])
as.data.frame(data.list2$HK.p.data[1:7,])
as.data.frame(data.list2$IND.p.data[1:7,])
as.data.frame(data.list2$PHI.p.data[1:7,])

collstr.th.pd <- as.data.frame(data.list2$CAN.p.data[1:3,]) %>% droplevels
collstr.th.pd <- rbind(collstr.th.pd, as.data.frame(data.list2$GB.p.data[1:3,]), as.data.frame(data.list2$HK.p.data[1:3,]), as.data.frame(data.list2$IND.p.data[1:3,]), as.data.frame(data.list2$IRE.p.data[1:3,]), as.data.frame(data.list2$JA.p.data[1:3,]), as.data.frame(data.list2$NZ.p.data[1:3,]), as.data.frame(data.list2$PHI.p.data[1:3,]), as.data.frame(data.list2$SIN.p.data[1:3,])) %>% droplevels


row.names(collstr.th.pd) <- NULL
collstr.th.pd$variety <- factor(collstr.th.pd$variety, levels=c("CAN", "GB", "IRE", "NZ","JA", "SIN", "HK", "IND", "PHI"))

# reorder to plot from highest to lowest collstrength
collstr.th.pd$words <- reorder(collstr.th.pd$words, -collstr.th.pd$coll.strength)

unique(collstr.th.pd$words)

(ex.collstr.th.pd <- ggplot(collstr.th.pd, aes(variety, coll.strength, fill=words)) +
  geom_bar(stat="identity", position=position_dodge(), color="black", width=.8) +
  geom_text(aes(label=round(coll.strength,1)), vjust=-0.4, position=position_dodge(0.8), size=4) +
  scale_fill_brewer(palette = "BrBG", direction=-1, name="theme") +
    labs(title="prepositional dative", y="collostructional strength", x="") +
    theme(axis.title.x = element_text(size=12), axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=16), plot.title = element_text(size=16, face="bold", hjust = 0.5)) +
    guides(fill=guide_legend(ncol=3,title.position = "top", title.hjust = 0.5))+
    theme(legend.position = "bottom",legend.title = element_text(size=16,face = "bold"), legend.text = element_text(size=16)))


pdf("Rplot_lexical_collostruction_examples_th_pd.pdf", width = 12, height = 7)
ex.collstr.th.pd
dev.off()


# Summary of top 10 themes and in which varieties they occur
addmargins(xtabs(~ words + variety, data=collstr.th.do))  
nrow(xtabs(~words+variety, data=collstr.th.do))


pdf("Rplot_lexical_collostruction_examples_th_allvariants.pdf", width = 12, height = 15)
grid.arrange(arrangeGrob(ex.collstr.th.do, ex.collstr.th.pd))
dev.off()



# -- type-token ratio by variant ---------------------------
# H: TTR is lower for recipients and themes in ditransitive datives in L2 compared to L1 

length(unique(data[data$Resp=="do" & data$Nativity=="L2", 'RecHeadLemma']))/length(data[data$Resp=="do" & data$Nativity=="L2", 'RecHeadLemma'])

length(unique(data[data$Resp=="do" & data$Nativity=="L1", 'RecHeadLemma']))/length(data[data$Resp=="do" & data$Nativity=="L2", 'RecHeadLemma'])

length(unique(data[data$Resp=="do" & data$Nativity=="L2", 'ThemeHeadLemma']))/length(data[data$Resp=="do" & data$Nativity=="L2", 'ThemeHeadLemma'])

length(unique(data[data$Resp=="do" & data$Nativity=="L1", 'ThemeHeadLemma']))/length(data[data$Resp=="do" & data$Nativity=="L1", 'ThemeHeadLemma'])


save.image("datives_lexical.RData")