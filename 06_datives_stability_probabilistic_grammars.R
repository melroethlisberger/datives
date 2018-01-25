#=======================================================================
#   Script: 06_datives_stability_probabilistic_grammars.R
#   Date last changed: 11 January 2018
#
#   Rversion: 3.3.2
#   Author: MR
#   This script provides the code to calculate the probabilistic distance
#   between varieties along the three lines of evidence proposed by 
#   Tagliamonte and colleagues (Comparative Sociolinguistic methods).
#=======================================================================


#----SETUP: LIBRARIES-----------------------------------------
library(reshape2);library(plyr);library(dplyr);library(magrittr);library(lme4);library(ggplot2);library(gridExtra); library(JGmisc); library(beepr); library(effects); library(Hmisc); library(data.table); library(MASS); library(vegan); library(cluster)

setwd("/Users/melanierothlisberger/Dropbox/Universitaet/04_Statistics/statistical analyses/20170209_THESIS_full_model/")
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


# -- compute by-variety models with 5 predictors ------------------------
f_rand <- Resp ~ (1|Verb) + (1|ThemeHeadFilter) + (1|RecHeadFilter) + (1|FileID)
f_five <- update.formula(f_rand, .~. +
                           z.logWeightRatio +
                           c.RecPron +
                           c.ThemeBinComplexity +
                           c.ThemePron +
                           z.ThemeHeadFreq)

vars <- levels(data$Variety)

model_list5 <- list()
model_params5 <- list()


t1 <- proc.time()
for (i in 1:9){
  print(paste("Working on", vars[i], "model..."))
  d <- droplevels(subset(data, data$Variety == vars[i]))
  
  cat("\tscaling and pruning factors...\n ")
  
  # Filter out infrequent verbs 
  d$Verb <- filter.infrequent(d$Verb, threshold = 5)
  d$ThemeHeadFilter <- filter.infrequent(d$ThemeHeadLemma, threshold = 5)
  d$RecHeadFilter <- filter.infrequent(d$RecHeadLemma, threshold = 5)
  
  # set reference levels if not done so to combine DO-friendly features in ref
  d$RecPron <- relevel(d$RecPron, ref="pron")
  d$ThemePron <- relevel(d$ThemePron, ref="non-pron")
  d$ThemeBinComplexity <- relevel(d$ThemeBinComplexity, ref="complex")

  
  # rescale predictors
  d$z.logWeightRatio = z.(d$logWeightRatio, factor=2)
  d$c.RecPron = c.(as.numeric(d$RecPron))
  d$c.ThemeBinComplexity = c.(as.numeric(d$ThemeBinComplexity))
  d$c.ThemePron = c.(as.numeric(d$ThemePron))
  d$z.ThemeHeadFreq = z.(d$ThemeHeadFreq, factor=2)

  
  cat("fitting model...\n")
  # model
  m <- glmer(f_five, data = d, family = binomial,
             control = glmerControl(optimizer = 'bobyqa',
                                    optCtrl = list(maxfun = 1e6)))
  
  model_list5[[i]] <- m
  names(model_list5)[i] <- paste(vars[i], "glmm", sep = '.')
  
  cat("add predictions ...\n")
  # add predictions
  d$preds <- predict(m)
  name <- paste(vars[i], "df", sep = '.') # make name for variety dataframe
  assign(name, d)
  
  parameters <- c(as.vector(summary(m)$coefficients[, 1]), # fixed effects
                  as.data.frame(summary(m)$varcor)[, "sdcor"] # random effects SDs
  )
  r2 <- MuMIn::r.squaredGLMM(m)
  model_params5[[i]] <- summary(m)$coefficients
  
  if(i == 1){
    # on first run create a dataframe of model parameters
    pars5 <- data.frame(row.names = c(rownames(summary(m)$coefficients), 
                                     names(summary(m)$varcor)))
    # add first model results to the dataframe
    pars5 <- cbind(pars5, parameters)
    
    # create matrix for r-squared
    squared5 <- data.frame(row.names = rownames(as.matrix(r2)))
    
    # add results of first r-squared to df
    squared5 <- cbind(squared5, as.matrix(r2))
  }
  else {
    pars5 <- cbind(pars5, parameters)
    squared5 <- cbind(squared5, as.matrix(r2))
  }
  rm("d") # remove var specific df (optional)
}
t2 <- proc.time() - t1

names(pars5) <- vars
round(pars5, 3)

for (i in 1:9){
  names(squared5)[i] <- paste(vars[i], "glmm", sep=".")
}
for (i in 1:9){
  names(model_params5)[i] <- paste(vars[i])
}


# -- .. evaluation of by-variety models MDS-------------
summary(model_list5$GB.glmm, correlation=FALSE) # ok
summary(model_list5$CAN.glmm, correlation=FALSE) # ok
summary(model_list5$HK.glmm, correlation=FALSE) # ok
summary(model_list5$IND.glmm, correlation=FALSE) # ok
summary(model_list5$IRE.glmm, correlation=FALSE) # ok, failed to converge with SpeakerID
summary(model_list5$JA.glmm, correlation=FALSE) #ok
summary(model_list5$NZ.glmm, correlation=FALSE) # ok
summary(model_list5$PHI.glmm, correlation=FALSE) # ok
summary(model_list5$SIN.glmm, correlation=FALSE) # ok


ggLogit.plot(model_list5$GB.glmm, GB.df) # good
ggLogit.plot(model_list5$HK.glmm, HK.df) # good
ggLogit.plot(model_list5$CAN.glmm, CAN.df) # bad
ggLogit.plot(model_list5$IND.glmm, IND.df) # good
ggLogit.plot(model_list5$IRE.glmm, IRE.df) # bad
ggLogit.plot(model_list5$JA.glmm, JA.df) # s, bad
ggLogit.plot(model_list5$NZ.glmm, NZ.df) # s
ggLogit.plot(model_list5$PHI.glmm, PHI.df) # bad
ggLogit.plot(model_list5$SIN.glmm, SIN.df) # good

GB.fits5 <- fitted(model_list5$GB.glmm)
sum(ifelse(GB.fits5 > .5, "pd", "do") == GB.df$Resp)/nrow(GB.df) # 90.7 %

CAN.fits5 <- fitted(model_list5$CAN.glmm)
sum(ifelse(CAN.fits5 > .5, "pd", "do") == CAN.df$Resp)/nrow(CAN.df) # 95.1 %

HK.fits5 <- fitted(model_list5$HK.glmm)
sum(ifelse(HK.fits5 > .5, "pd", "do") == HK.df$Resp)/nrow(HK.df) # 92.6 %

IND.fits5 <- fitted(model_list5$IND.glmm)
sum(ifelse(IND.fits5 > .5, "pd", "do") == IND.df$Resp)/nrow(IND.df) # 93.4 %

IRE.fits5 <- fitted(model_list5$IRE.glmm)
sum(ifelse(IRE.fits5 > .5, "pd", "do") == IRE.df$Resp)/nrow(IRE.df) # 94.5 %

JA.fits5 <- fitted(model_list5$JA.glmm)
sum(ifelse(JA.fits5 > .5, "pd", "do") == JA.df$Resp)/nrow(JA.df) # 94.8 %

NZ.fits5 <- fitted(model_list5$NZ.glmm)
sum(ifelse(NZ.fits5 > .5, "pd", "do") == NZ.df$Resp)/nrow(NZ.df) # 92.5 %

PHI.fits5 <- fitted(model_list5$PHI.glmm)
sum(ifelse(PHI.fits5 > .5, "pd", "do") == PHI.df$Resp)/nrow(PHI.df) # 94.9 %

SIN.fits5 <- fitted(model_list5$SIN.glmm)
sum(ifelse(SIN.fits5 > .5, "pd", "do") == SIN.df$Resp)/nrow(SIN.df) # 90.7 %

somers.C5 <- list()
for(i in 1:9) {
  somers.C5 <- cbind(somers.C5, somers.mer(model_list5[[i]])[[1]])
}
somers.C5 <- as.data.frame(somers.C5)
names(somers.C5) <- vars

for(i in 1:9){
  print(paste(overdisp.mer(model_list5[[i]])))
}





# -- 1st line of evidence: statistical significance --------
(significance_matrix <- matrix(nrow = length(rownames(model_params5[[1]])),
                               ncol = length(vars),
                               dimnames = list(rownames(model_params5[[1]]),
                                               vars),
                               data = sapply(model_params5,
                                             function(x) {
                                               ifelse(x[, 4] < .05,
                                                      TRUE, FALSE)
                                             }
                               )))

significance_matrix <- significance_matrix[-1, ] # Ignoring the intercept

# Euclidean distances
(significance_dist <- dist(t(significance_matrix), method = "euclidean"))
(significance_dist <- significance_dist^2) # squared Euclidean distances

# creating similarity matrix
(significance_sim <- 5-significance_dist)



# Putting distances in data table (long format)
dt_dist_sig <- data.table(t(combn(rownames(t(significance_matrix)), 2)), 
                          as.numeric(significance_dist))
names(dt_dist_sig) <- c("var1", "var2", "distance")

# Doubling dt_dist_sig in order to compute mean distances
dt_dist_sig2 <- rbind(dt_dist_sig, dt_dist_sig[, .(var1 = var2, var2 = var1, 
                                                   distance)])

# Computing mean distances
dt_dist_sig2_means <- dt_dist_sig2[, .(Mean_Dist = round(mean(distance), 2)), 
                                   by = (Variety = var1)][
                                     order(Mean_Dist, decreasing = T)]

# Computing scaled version (between 0 and 1), 1 being the most similar
n_comparisons <- nrow(significance_matrix)
dt_dist_sig2[, similarity := n_comparisons - distance][]
dt_dist_sig2[, sim_scaled := similarity / n_comparisons][]

# Computing mean similarities
dt_sim_sig2_means <- dt_dist_sig2[, .(Mean_Similarity = round(mean(sim_scaled), 3)), 
                                   by = (Variety = var1)][
                                     order(Mean_Similarity, decreasing = T)]
dt_sim_sig2_means

# Saving mean similarities
write.table(dt_sim_sig2_means, "1st_line_mean_similarities5.txt", quote = F)

# Finally, computing single coefficient
round(dt_dist_sig2[, mean(sim_scaled)], 3)


# -- .. plot MDS --------------------------------------
dist1 <- dist(t(significance_matrix), method="manhattan")
# replace 0 values with very small positive numbers (Levshina 345)
dist1b <- dist1
dist1b[dist1b == 0] <- 0.0001
require(vegan)
fit <- isoMDS(dist1b, k=2)

# check out reduce in stress vs. dimensions:
str <- vector()
for (i in 1:5){
  f <- isoMDS(dist1b, k=i)
  s <- f$stress
  str <- cbind(str, s)
}


# Shepard's plot:
require(MASS)
sign.sh <- Shepard(dist1b, fit$points)
plot(sign.sh, main="Shepard plot", pch=".")
lines(sign.sh$x, sign.sh$yf, type="S")

# use smacof due to ties in distances:
require(smacof)
smac <- smacofSym(dist1b, type="interval", ties="secondary") # it is interval data
plot(smac) # the plot is exactly the same as when I use cmdscale!
summary(smac)
smac$conf # x and y axis to plot

# 2 Dimensions
fit2d <- as.data.frame(fit$points) # the same picture for all
names(fit2d) <- c("x", "y")
# relabel to fit other data
rownames(fit2d) <- c("CanE", "BrE", "HKE", "IndE", "IrE", "JamE", "NZE", "PhiE", "SinE")
require(ggrepel)
(plot.mds1st <- ggplot(fit2d, aes(x=x, y=y, label = rownames(fit2d))) + 
    #geom_point(colour=brewerblue, size=5) +
    geom_text_repel(size=6.5) +
    theme(axis.title = element_text(size=17))
  + labs(x = "Dimension 1", y='Dimension 2'))

pdf("Rplot_1stline_mds.pdf", width=6, height=4)
plot.mds1st
dev.off()



# -- 2nd line of evidence: coefficient estimates -----------
# Effect sizes from model
(param_matrix <- matrix(nrow = length(rownames(model_params5[[1]])), 
                        ncol = length(vars), 
                        dimnames = list(rownames(model_params5[[1]]), 
                                        vars),
                        data = sapply(model_params5, function(x) { 
                          x[, 1] 
                        } )))

param_matrix <- param_matrix[-1, ] # Ignoring the intercept
round(param_matrix,3)
# Euclidean distances
(param_strength <- dist(t(param_matrix), method = "euclidean"))
#(param_strength <- param_strength^2) # squared Euclidean distance
round(param_strength,3)


# Putting distances in data table
dt_dist_strength <- data.table(t(combn(rownames(t(param_matrix)),2)), 
                               as.numeric(param_strength))
names(dt_dist_strength) <- c("var1", "var2", "distance")

# Doubling dt_dist_strength in order to compute mean distances
dt_dist_strength2 <- rbind(dt_dist_strength, dt_dist_strength[, .(var1 = var2, 
                                                                  var2 = var1, 
                                                                  distance)])

# Computing mean distances
dt_dist_strength2_means <- dt_dist_strength2[, .(Mean_Dist = round(mean(distance), 3)), 
                                             by = (Variety = var1)][
                                               order(Mean_Dist, decreasing = T)]
dt_dist_strength2_means

# Saving mean distances
write.table(dt_dist_strength2_means, "2nd_line_mean_dist.txt", quote = F)

# Distance --> similarity

# The question is: What is maximally different? Let's say it is if there are 
# no effects (the problem with Euclidean distance is that there is no upper threshold and no indication how distant two groups are wrt to their maximal distance, we only know that 0=most similar.)

# How much does a model in which nothing has an effect differ on average from 
# the others?

param_matrix_null <- cbind(param_matrix, rep(0, nrow(param_matrix)))
colnames(param_matrix_null)[length(vars) + 1] <- "NULL"
(param_strength_null <- dist(t(param_matrix_null), method = "euclidean"))
#(param_strength_null <- param_strength_null^2) # squared Euclidean distance

# create data table for ease of calculations
dt_dist_strength_null <- data.table(t(combn(rownames(t(param_matrix_null)), 2)), 
                                    as.numeric(param_strength_null))
names(dt_dist_strength_null) <- c("var1", "var2", "distance")

# Doubling param_strength_null in order to compute mean distances
dt_dist_strength_null2 <- rbind(dt_dist_strength_null, dt_dist_strength_null[, .(var1 = var2, var2 = var1, distance)])
dt_dist_strength_null_means <- dt_dist_strength_null2[, .(Mean_Dist = round(mean(distance), 3)),
                                                     by = (Variety = var1)][
                                                       order(Mean_Dist, decreasing = T)]
dt_dist_strength_null_means

# The transformation to similarities can now be done by normalizing by 
# distances to the null variety (for the original distance scores without NULL)
dt_dist_strength2_means[, Similarity := 1 - (Mean_Dist / max(dt_dist_strength_null_means$Mean_Dist))]
dt_dist_strength2_means[order(Similarity)]
dt_dist_strength2_means$Variety <- factor(dt_dist_strength2_means$Variety,
                                          levels = dt_dist_strength2_means[
                                            order(Similarity, decreasing = T), 
                                            Variety])


# Saving mean similarities
dt_dist_strength2_means[order(Similarity, decreasing = T), .(Variety, round(Similarity, 3))]
write.table(dt_dist_strength2_means[order(Similarity, decreasing = T), .(Variety, round(Similarity, 3))], 
            "2nd_line_mean_similarity_withNULL.txt", quote = F)

# Finally, computing single coefficient
dt_dist_strength2_means[, round(mean(Similarity), 4)]

dt_dist_strength2_means[,c(1,3)]



#--..plot MDS ---------------------------------
dist2 <- dist(t(param_matrix), method = "euclidean")
fit <- cmdscale(dist2, eig = TRUE, k=2)
eig_df <- data.frame(Eigenvalue = fit$eig, Dim=factor(1:length(fit$eig)))
ggplot(eig_df, aes(Dim, Eigenvalue)) +
  geom_bar(stat="identity", fill = "steelblue", color="black", width = .7) +
  geom_line(aes(group = 1)) + geom_point()

eig <- fit$eig
var.dim1 <- fit$eig[1]/sum(fit$eig) # X
var.dim2 <- fit$eig[2]/sum(fit$eig) # Y
var.dim3 <- fit$eig[3]/sum(fit$eig) # Z
eig/sum(eig)
cumsum(eig/sum(eig)) # the first two dimensions accounts for 89.0% of the variance

# GOF = 0.8901736
# stress:
sqrt(sum((dist2 - dist(fit$points))^2)/sum(dist2^2))




# 2 Dimensions needed
fit2d <- as.data.frame(fit[[1]])
names(fit2d) <- c("x", "y")
# relabel to fit other data
rownames(fit2d) <- c("CanE", "BrE", "HKE", "IndE", "IrE", "JamE", "NZE", "PhiE", "SinE")
require(ggrepel)
(plot.mds2nd <- ggplot(fit2d, aes(x=x, y=y, label = rownames(fit2d))) + 
    #geom_point(colour='black', size=5) +
    geom_text(size=6.5) +
    theme(axis.title = element_text(size=17)) +
    labs(x = "53.3% variance", y='35.8% variance'))

pdf("Rplot_2ndline_mds.pdf", width=9.5, height=7)
plot.mds2nd
dev.off()



# -- 3rd line of evidence: ranking of predictors -----------
f <- formula(Resp ~      # using same predictors as in glmer                         
               z.logWeightRatio +
               RecPron +
               ThemeBinComplexity +
               ThemePron +
               z.ThemeHeadFreq)

# Control settings
library(party);set.seed(12578)
forest.controls = cforest_unbiased(ntree=2000, mtry=3)

rf <- cforest(f, data = data, controls = forest.controls)
varimpAUC.rf <- party::varimpAUC(rf)


# ..Building Per-Variety Random Forest Models -----------------------------
forest.controls = cforest_unbiased(ntree=2000, mtry=3)
set.seed(230573)

varimpAUC.rf <- as.data.frame(varimpAUC.rf)

# Setting up a data table that captures variable importance values
varimp <- data.table(expand.grid(pred = rownames(varimpAUC.rf), 
                                 var = levels(data$Variety)))
varimp[, MeanDecreaseVarimp := 0][] # Adding column for varimp values


t1 <- proc.time()
for(i in 1:length(levels(data$Variety))) {
  v  <- levels(data$Variety)[i]
  cat("Fitting model on", v,  "\n")
  d  <- subset(data, data$Variety == v) %>% droplevels
  
  # rescale predictors
  d$z.logWeightRatio = z.(d[,"logWeightRatio"], factor=2)
  d$z.ThemeHeadFreq = z.(d[, "ThemeHeadFreq"], factor =2)

  # run model
  rf <- cforest(f, data = d, controls = forest.controls)
  assign(paste(v, "rf", sep = '.'), rf)
  
  # add predictions
  d$preds <- predict(rf)
  name <- paste(v, "df", sep= '.')
  assign(name,d)
  
  # get varimp
  varimp <- party::varimpAUC(rf)
  assign(paste(v, "varimp", sep = '.'), varimp)
  
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
names(pars.crf.rk) <- vars
# transform pars.crf to long data format
crf_long <- melt(t(pars.crf.chisq))
names(crf_long) <- c("Variety", "Factor", "Varimp")

crf_long.dt <- data.table(crf_long)
rk_long <- melt(t(pars.crf.rk)); names(rk_long) <- c("Variety", "Factor", "Rank")
rk_long.dt <- data.table(rk_long)


# ..Calculating Pairwise Correlations -------------------------------------

# Getting global order of predictors
(pred_order <- crf_long.dt[, .(imp = mean(Varimp)), by=Factor][order(imp, decreasing=T), Factor])

# Releveling varimp accordingly
#varimp[, pred := factor(pred, levels = pred_order)]
pairs <- combn(levels(crf_long$Variety), 2)

# Setting up a data frame that holds pairwise correlation coefficients
pairwise_correlations <- data.frame(
  matrix(nrow = ncol(pairs), ncol = 4, 
         dimnames = list(NULL, c("VAR1", "VAR2", "COR", "P"))))

# Filling the data frame
for(i in 1:ncol(pairs)) { # For all variety combinations ...
  var1 <- pairs[, i][1]
  var2 <- pairs[, i][2]
  
  # Compare the rankings from rk_long.dt by sorting the values according to
  # the global order to make sure that both vectors have the same ordering
  var1_values <- crf_long.dt[Variety==var1][order(match(Factor, pred_order)), 
                                   Varimp]
  var2_values <- crf_long.dt[Variety==var2][order(match(Factor, pred_order)), 
                                   Varimp]
  
  # Determine correlation
  correlation <- cor.test(rank(var1_values), rank(var2_values), method = "spearman")
  
  # Saving values
  pairwise_correlations[i, ] <- c(var1, var2, correlation$estimate, 
                                  correlation$p.value)
}

# to test whether the cor.test from the rk.dt would have given the same:
#cor.test(pars.crf.rk[,1], pars.crf.rk[,2], method="spearman")

# ..Evaluation ------------------------------------------------------------

# Overall mean similarity of predictor rankings
pairwise_correlations <- data.table(pairwise_correlations)
pairwise_correlations[, mean(as.numeric(COR))] # 0.839 (BH's data: 0.958)

# Doubling data table to calculate mean correlations per variety
pairwise_correlations2 <- rbind(pairwise_correlations,
                                pairwise_correlations[, .(VAR1=VAR2, 
                                                          VAR2=VAR1, COR, P)])
mean_sim <- pairwise_correlations2[, .(Mean_Cor = mean(as.numeric(COR))), by = VAR1][
  order(Mean_Cor, decreasing = T)]

mean_sim

# Pairwise correlations to matrix
pairwise_correlations_wide <- dcast(pairwise_correlations2, VAR1 ~ VAR2, value.var = "COR")
rownames <- pairwise_correlations_wide$VAR1 # Saving rownames
# Ignore first column
pairwise_correlations_matrix <- as.matrix(pairwise_correlations_wide[, -1, with = F])
# Set row names
rownames(pairwise_correlations_matrix) <- rownames
# Get rid of NAs
pairwise_correlations_matrix <- ifelse(is.na(pairwise_correlations_matrix), 1, 
                                       pairwise_correlations_matrix)
# Convert matrix to numeric
class(pairwise_correlations_matrix) <- "numeric"
# Display
round(as.dist(pairwise_correlations_matrix), 3)
round(pairwise_correlations_matrix[,-1],3)

# Saving 
write.table(round(pairwise_correlations_matrix, 3), "3rd_line_pc_matrix.txt", quote = F)

# ..display ranks compared to global ranking --------------------------
ranking <- pars.crf.rk
ranking$global <- order(match(pred_order, rownames(ranking)))
# order according to global predictor ranking:
ranking[order(ranking$global, decreasing = F),]
# see differences from global average:
ranking[order(ranking$global, decreasing = F),]-(ranking[order(ranking$global, decreasing = F),]$global)


# ..plot MDS -------------------------------------------------------------------

# Distance matrix = inverse of pairwise correlations
dist_matrix <- 1 - pairwise_correlations_matrix

# As with dist(), only focus on lower triangle without diagonal in MDS
dist_matrix <- as.dist(dist_matrix)
# add very low values to replace 0:
dist_matrix[dist_matrix == 0] <- 0.0001

fit <- isoMDS(dist_matrix, k = 3)

# stress = 0.007
fit$points

# 3 Dimensions needed
fit.df3d <- as.data.frame(fit$points)
names(fit.df3d) <- c("x","y", "z")


# change labels for variety to make more readable in the output
rownames(fit.df3d) <- c("CanE", "BrE", "HKE", "IndE", "IrE", "JamE", "NZE", "PhiE", "SinE")

# this creates a static 3D plot with vertical lines from the data points
require(scatterplot3d)
source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')
#par(pty='m', new=TRUE) # to have axes with the same length 
with(fit.df3d, {
  s3d <- scatterplot3d(x,   # x axis
                       y,     # y axis
                       z,    # z axis
                       main="", 
                       xlab = "",
                       ylab = '',
                       zlab = '',
                       pch="", # filled black circles
                       #type="h", # horizontal lines 
                       asp=NA,
                       angle=55,
                       grid=F,
                       box=F,
                       type="h",lty.hplot=1,
                       x.ticklabs = '',
                       y.ticklabs='',
                       z.ticklabs='')
  addgrids3d(fit.df3d[, 1:3], grid = c("xy", "xz", "yz"), angle=55, col.grid = "gainsboro")
  s3d$points3d(fit.df3d[, 1:3], pch = 21, col="black", bg="dimgray", type="h")
  s3d.coords <- s3d$xyz.convert(x, y, z) # convert 3D coords to 2D projection
  text(s3d.coords$x, s3d.coords$y,             # x and y coordinates
       labels=row.names(fit.df3d), # text to plot
       cex=1.2, pos=4)           # shrink text 50% and place to right of points
})


pdf("Rplot_3rdline_mds_3D.pdf", width=7, height=5)
plot.mds3rd.3d
dev.off()



# -- Correlation between the three distance matrices --------
require(vegan)
# dist1 - significance, dist2 - effect size, dist_matrix - ranking
sign <- as.matrix(dist1); row.names(sign) <- NULL; colnames(sign) <- NULL
eft <- as.matrix(dist2); row.names(eft) <- NULL; colnames(eft) <- NULL
ranks <- as.matrix(dist_matrix); row.names(ranks) <- NULL; colnames(ranks) <- NULL
mantel(sign,eft, method="pearson")
cor(sign, eft)
mantel(eft,ranks)
mantel(sign, ranks)
mantel.partial(sign, eft, ranks)

# -- Correlation with HVE -------------------
hve <- read.csv("Features_BS/table_binary.csv", header=T, sep=";", row.names=1)
names(hve)
rownames(hve)
colnames(hve)
head(hve)
str(hve)

# -- Correlation with ewave ---------------
#ewave <- read.delim("Features_BS/ewave.txt", header=T)
ewave <- read.csv("Features_BS/ewave_short.csv", header=T, sep="\t")
names(ewave)
rownames(ewave)
colnames(ewave)
head(ewave)
str(ewave)


#--..plot MDS--------------------------------------------

MDSdata <- ewave
MDSdata <- MDSdata[-(1:4), ] # delete rows with variety info (type, region) (for hve = rows 1-2)
rownames(MDSdata) <- NULL
rownames(MDSdata) <- MDSdata$Feature
head(MDSdata)
MDSdata <- MDSdata[,-1] # get rid of features
# combine A+B as attested (1) and C,D,X and ? as not attested (0)

# Ordinal predictor (Levshina, p. 343)
MDSdata_ord <- MDSdata
MDSdata_ord[MDSdata_ord=="?"] <- NA
MDSdata_ord[MDSdata_ord=="X"] <- NA
MDSdata_ord <- lapply(MDSdata_ord[,1:79], function(x) ordered(x, levels=c("D", "C", "B", "A")))
MDSdata_ord <- data.frame(MDSdata_ord)
str(MDSdata_ord)


# use ordinal scale
MDSdata_transposed <- data.frame(t(MDSdata_ord))
str(MDSdata_transposed)

# restrict attention to nine varieties, incl. Southwest for BrE
MDSrestr <- subset(MDSdata_transposed, rownames(MDSdata_transposed) %in% c("IrE", "SE", "JamE", "IndE", "NZE", "HKE", "CollSgE", "PhilE"))

# Gower distance to take ordinal scales into account
require(cluster)
ling.dist <- daisy(MDSrestr, metric = "gower")

summary(ling.dist)
ling.mds <- isoMDS(ling.dist, k=2)

fit.gower <- data.frame(ling.mds$points)
rownames(fit.gower) <- c("IrE", "BrE", "JamE", "IndE", "HKE", "SinE", "PhiE", "NZE")
(ewave.plot.gower <- ggplot(fit.gower, aes(x=X1, y=X2, label = rownames(fit.gower))) + 
    #geom_point(colour=brewerblue, size=5) +
    geom_text(size=6.5) +
    theme(axis.title = element_text(size=17)) +
    labs(x='Dimension 1', y='Dimension 2'))


pdf("Rplot_distances_mds_ewave_gower.pdf", width=8, height=6)
ewave.plot.gower
dev.off()


#-- ..correlation between distance matrices----------
require(vegan)

# create distance matrices without CanE
# 1st line
cor.sign <- significance_matrix[,-c(1)]
d1 <- dist(t(cor.sign), method="manhattan")
m.d1 <- as.matrix(d1); row.names(m.d1) = NULL; colnames(m.d1) = NULL

# 2nd line
cor.eff <- param_matrix[,-c(1)]
d2 <- dist(t(cor.eff), method = "euclidean")
m.d2 <- as.matrix(d2); row.names(m.d2) = NULL; colnames(m.d2) = NULL

# 3rd line
dist_matrix <- 1 - pairwise_correlations_matrix
cor.rk <- dist_matrix[-c(1),-c(1)]
d3 <- as.dist(cor.rk)
m.d3 <- as.matrix(d3); row.names(m.d3) = NULL; colnames(m.d3) = NULL

# order ewave matrix the same way as other distance matrices re Varieties
rownames(MDSrestr) # should be GB-HK-IND-IRE-JA-NZ-PHI-SIN
MDSreordered <- as.data.frame(rbind(MDSrestr[2,], MDSrestr[5,], MDSrestr[4,], MDSrestr[1,], MDSrestr[3,], MDSrestr[8,], MDSrestr[7,], MDSrestr[6,])); rownames(MDSreordered) <- rownames(cor.rk)

d4 <- daisy(MDSreordered)

wave.dist <- as.matrix(d4); row.names(wave.dist) <- NULL; colnames(wave.dist) <- NULL


mantel(m.d1,wave.dist)
mantel(m.d2,wave.dist)
mantel(m.d3, wave.dist)



# -- DISCUSSION: random-intercept plots by speaker ID-----------
# ==> use results from speaker ID to plot random intercept adjustments

intercept_df <- model_list5 #(with SpeakerIDs, exclude IRE for plot since model didn't converge for SpeakerID

# get standard deviations per random intercept for Speaker
sds.CAN <- round(attributes(VarCorr(intercept_df$CAN.glmm)$SpeakerID)$stddev[[1]],2)
sds.GB <- round(attributes(VarCorr(intercept_df$GB.glmm)$SpeakerID)$stddev[[1]],2)
sds.HK <- round(attributes(VarCorr(intercept_df$HK.glmm)$SpeakerID)$stddev[[1]],2)
sds.IND <- round(attributes(VarCorr(intercept_df$IND.glmm)$SpeakerID)$stddev[[1]],2)
sds.JA <- round(attributes(VarCorr(intercept_df$JA.glmm)$SpeakerID)$stddev[[1]],2)
sds.NZ <- round(attributes(VarCorr(intercept_df$NZ.glmm)$SpeakerID)$stddev[[1]],2)
sds.PHI <- round(attributes(VarCorr(intercept_df$PHI.glmm)$SpeakerID)$stddev[[1]],2)
sds.SIN <- round(attributes(VarCorr(intercept_df$SIN.glmm)$SpeakerID)$stddev[[1]],2)

sds.all <- c(sds.CAN, sds.GB, sds.HK, sds.IND, sds.JA, sds.NZ, sds.PHI, sds.SIN)


r.CAN <- ranef(intercept_df$CAN.glmm)[[1]] # get by-speaker adjustements for CAN
r.CAN <- as.data.frame(cbind(Variety = rep("CAN", nrow(r.CAN)), r.CAN))
r.GB <- ranef(intercept_df$GB.glmm)[[1]] # get by-speaker adjustements for GB
r.GB <- as.data.frame(cbind(Variety = rep("GB", nrow(r.GB)), r.GB))
r.HK <- ranef(intercept_df$HK.glmm)[[1]] # get by-speaker adjustements for HK
r.HK <- as.data.frame(cbind(Variety = rep("HK", nrow(r.HK)), r.HK))
r.IND <- ranef(intercept_df$IND.glmm)[[1]] # get by-speaker adjustements for IND
r.IND <- as.data.frame(cbind(Variety = rep("IND", nrow(r.IND)), r.IND))
r.JA <- ranef(intercept_df$JA.glmm)[[1]] # get by-speaker adjustements for JA
r.JA <- as.data.frame(cbind(Variety = rep("JA", nrow(r.JA)), r.JA))
r.NZ <- ranef(intercept_df$NZ.glmm)[[1]] # get by-speaker adjustements for NZ
r.NZ <- as.data.frame(cbind(Variety = rep("NZ", nrow(r.NZ)), r.NZ))
r.PHI <- ranef(intercept_df$PHI.glmm)[[1]] # get by-speaker adjustements for PHI
r.PHI <- as.data.frame(cbind(Variety = rep("PHI", nrow(r.PHI)), r.PHI))
r.SIN <- ranef(intercept_df$SIN.glmm)[[1]] # get by-speaker adjustements for SIN
r.SIN <- as.data.frame(cbind(Variety = rep("SIN", nrow(r.SIN)), r.SIN))

r.all <- rbind(r.CAN, r.GB, r.HK, r.IND, r.JA, r.NZ, r.PHI, r.SIN)
names(r.all)
names(r.all)[2] <- "Intercept"
str(r.all)
#r.all$Variety <- factor(td02.df$Morph, levels = c("semi-weak", "monomorphemes", "regular past"))

varieties <- levels(r.all$Variety)

#--.. Plot the BLUBs ----------

(speaker.plot <- ggplot(r.all, aes(x = Variety, y = Intercept)) +
   geom_point(position= position_jitter(width = .15, height = 0), color="gray50", alpha=.5) +
   geom_boxplot(width = .5, alpha = .7) +
   theme(legend.position = "none") +
   annotate("text", x = varieties, y = rep(2.2, 8),
            label = paste("SD =", sds.all),
            col = 'black', 
            size=5) +
   labs(x = "", y = "by speaker adjustment") +
   theme(axis.text.x = element_text(size=18), axis.title.y = element_text(size=15)))

pdf("Rplot_by_speaker_intercept_boxplots.pdf", width = 10, height = 6.5)
speaker.plot
dev.off()

# -- DISCUSSION: random-intercept plots by recipient-----------
# ==> use results from speaker ID to plot random intercept adjustments

intercept_rec <- model_list5 # with FileID to include IrE as well

# get standard deviations per random intercept for Speaker
sds.CAN <- round(attributes(VarCorr(intercept_rec$CAN.glmm)$RecHeadFilter)$stddev[[1]],2)
sds.GB <- round(attributes(VarCorr(intercept_rec$GB.glmm)$RecHeadFilter)$stddev[[1]],2)
sds.HK <- round(attributes(VarCorr(intercept_rec$HK.glmm)$RecHeadFilter)$stddev[[1]],2)
sds.IND <- round(attributes(VarCorr(intercept_rec$IND.glmm)$RecHeadFilter)$stddev[[1]],2)
sds.IRE <- round(attributes(VarCorr(intercept_rec$IRE.glmm)$RecHeadFilter)$stddev[[1]],2)
sds.JA <- round(attributes(VarCorr(intercept_rec$JA.glmm)$RecHeadFilter)$stddev[[1]],2)
sds.NZ <- round(attributes(VarCorr(intercept_rec$NZ.glmm)$RecHeadFilter)$stddev[[1]],2)
sds.PHI <- round(attributes(VarCorr(intercept_rec$PHI.glmm)$RecHeadFilter)$stddev[[1]],2)
sds.SIN <- round(attributes(VarCorr(intercept_rec$SIN.glmm)$RecHeadFilter)$stddev[[1]],2)

sds.all.rec <- c(sds.CAN, sds.GB, sds.HK, sds.IND, sds.IRE, sds.JA, sds.NZ, sds.PHI, sds.SIN)


r.CAN <- ranef(intercept_rec$CAN.glmm)[["RecHeadFilter"]] # get by-rec adjustements for CAN
r.CAN <- as.data.frame(cbind(Variety = rep("CAN", nrow(r.CAN)), r.CAN))
r.GB <- ranef(intercept_rec$GB.glmm)[["RecHeadFilter"]] # get by-rec adjustements for GB
r.GB <- as.data.frame(cbind(Variety = rep("GB", nrow(r.GB)), r.GB))
r.HK <- ranef(intercept_rec$HK.glmm)[["RecHeadFilter"]] # get by-rec adjustements for HK
r.HK <- as.data.frame(cbind(Variety = rep("HK", nrow(r.HK)), r.HK))
r.IND <- ranef(intercept_rec$IND.glmm)[["RecHeadFilter"]] # get by-rec adjustements for IND
r.IND <- as.data.frame(cbind(Variety = rep("IND", nrow(r.IND)), r.IND))
r.IRE <- ranef(intercept_rec$IRE.glmm)[["RecHeadFilter"]] # get by-rec adjustements for IRE
r.IRE <- as.data.frame(cbind(Variety = rep("IRE", nrow(r.IRE)), r.IRE))
r.JA <- ranef(intercept_rec$JA.glmm)[["RecHeadFilter"]] # get by-rec adjustements for JA
r.JA <- as.data.frame(cbind(Variety = rep("JA", nrow(r.JA)), r.JA))
r.NZ <- ranef(intercept_rec$NZ.glmm)[["RecHeadFilter"]] # get by-rec adjustements for NZ
r.NZ <- as.data.frame(cbind(Variety = rep("NZ", nrow(r.NZ)), r.NZ))
r.PHI <- ranef(intercept_rec$PHI.glmm)[["RecHeadFilter"]] # get by-rec adjustements for PHI
r.PHI <- as.data.frame(cbind(Variety = rep("PHI", nrow(r.PHI)), r.PHI))
r.SIN <- ranef(intercept_rec$SIN.glmm)[["RecHeadFilter"]] # get by-rec adjustements for SIN
r.SIN <- as.data.frame(cbind(Variety = rep("SIN", nrow(r.SIN)), r.SIN))

r.all.rec <- rbind(r.CAN, r.GB, r.HK, r.IND, r.IRE, r.JA, r.NZ, r.PHI, r.SIN)
names(r.all.rec)
names(r.all.rec)[2] <- "Intercept"
str(r.all.rec)

varieties <- levels(r.all.rec$Variety)

#--.. Plot the BLUBs ----------

(recipient.plot <- ggplot(r.all.rec, aes(x = Variety, y = Intercept)) +
   geom_point(position= position_jitter(width = .15, height = 0), color="gray50", alpha=.5) +
   geom_boxplot(width = .5, alpha = .7) +
   theme(legend.position = "none") +
   annotate("text", x = varieties, y = rep(0.8, 9),
            label = paste("SD =", sds.all.rec),
            col = 'black', 
            size=5) +
   labs(x = "", y = "by-recipient adjustment", title="Random intercept by recipient") +
   theme(axis.text.x = element_text(size=18), axis.title.y = element_text(size=15), plot.title = element_text(size=16, face="bold")))

pdf("Rplot_by_recipient_intercept_boxplots.pdf", width = 10, height = 6.5)
recipient.plot
dev.off()

# -- DISCUSSION: random-intercept plots by theme-----------
# ==> use results from speaker ID to plot random intercept adjustments

intercept_df <- model_list5 # with FileID to include IrE as well

# get standard deviations per random intercept for Speaker
sds.CAN <- round(attributes(VarCorr(intercept_df$CAN.glmm)$ThemeHeadFilter)$stddev[[1]],2)
sds.GB <- round(attributes(VarCorr(intercept_df$GB.glmm)$ThemeHeadFilter)$stddev[[1]],2)
sds.HK <- round(attributes(VarCorr(intercept_df$HK.glmm)$ThemeHeadFilter)$stddev[[1]],2)
sds.IND <- round(attributes(VarCorr(intercept_df$IND.glmm)$ThemeHeadFilter)$stddev[[1]],2)
sds.IRE <- round(attributes(VarCorr(intercept_df$IRE.glmm)$ThemeHeadFilter)$stddev[[1]],2)
sds.JA <- round(attributes(VarCorr(intercept_df$JA.glmm)$ThemeHeadFilter)$stddev[[1]],2)
sds.NZ <- round(attributes(VarCorr(intercept_df$NZ.glmm)$ThemeHeadFilter)$stddev[[1]],2)
sds.PHI <- round(attributes(VarCorr(intercept_df$PHI.glmm)$ThemeHeadFilter)$stddev[[1]],2)
sds.SIN <- round(attributes(VarCorr(intercept_df$SIN.glmm)$ThemeHeadFilter)$stddev[[1]],2)

sds.all.th <- c(sds.CAN, sds.GB, sds.HK, sds.IND, sds.IRE, sds.JA, sds.NZ, sds.PHI, sds.SIN)


r.CAN <- ranef(intercept_df$CAN.glmm)[["ThemeHeadFilter"]] # get by-th adjustements for CAN
r.CAN <- as.data.frame(cbind(Variety = rep("CAN", nrow(r.CAN)), r.CAN))
r.GB <- ranef(intercept_df$GB.glmm)[["ThemeHeadFilter"]] # get by-th adjustements for GB
r.GB <- as.data.frame(cbind(Variety = rep("GB", nrow(r.GB)), r.GB))
r.HK <- ranef(intercept_df$HK.glmm)[["ThemeHeadFilter"]] # get by-th adjustements for HK
r.HK <- as.data.frame(cbind(Variety = rep("HK", nrow(r.HK)), r.HK))
r.IND <- ranef(intercept_df$IND.glmm)[["ThemeHeadFilter"]] # get by-th adjustements for IND
r.IND <- as.data.frame(cbind(Variety = rep("IND", nrow(r.IND)), r.IND))
r.IRE <- ranef(intercept_df$IRE.glmm)[["ThemeHeadFilter"]] # get by-th adjustements for IRE
r.IRE <- as.data.frame(cbind(Variety = rep("IRE", nrow(r.IRE)), r.IRE))
r.JA <- ranef(intercept_df$JA.glmm)[["ThemeHeadFilter"]] # get by-th adjustements for JA
r.JA <- as.data.frame(cbind(Variety = rep("JA", nrow(r.JA)), r.JA))
r.NZ <- ranef(intercept_df$NZ.glmm)[["ThemeHeadFilter"]] # get by-th adjustements for NZ
r.NZ <- as.data.frame(cbind(Variety = rep("NZ", nrow(r.NZ)), r.NZ))
r.PHI <- ranef(intercept_df$PHI.glmm)[["ThemeHeadFilter"]] # get by-th adjustements for PHI
r.PHI <- as.data.frame(cbind(Variety = rep("PHI", nrow(r.PHI)), r.PHI))
r.SIN <- ranef(intercept_df$SIN.glmm)[["ThemeHeadFilter"]] # get by-th adjustements for SIN
r.SIN <- as.data.frame(cbind(Variety = rep("SIN", nrow(r.SIN)), r.SIN))

r.all.th <- rbind(r.CAN, r.GB, r.HK, r.IND, r.IRE, r.JA, r.NZ, r.PHI, r.SIN)
names(r.all.th)
names(r.all.th)[2] <- "Intercept"
str(r.all.th)

varieties <- levels(r.all.th$Variety)

#--.. Plot the BLUBs ----------

(theme.plot <- ggplot(r.all.th, aes(x = Variety, y = Intercept)) +
   geom_point(position= position_jitter(width = .15, height = 0), color="gray50", alpha=.5) +
   geom_boxplot(width = .5, alpha = .7) +
   theme(legend.position = "none") +
   annotate("text", x = varieties, y = rep(8, 9),
            label = paste("SD =", sds.all.th),
            col = 'black', 
            size=5) +
   labs(x = "", y = "by-theme adjustment", title="Random intercept by theme") +
   theme(axis.text.x = element_text(size=18), axis.title.y = element_text(size=15), plot.title = element_text(size=16, face="bold")))

pdf("Rplot_by_theme_intercept_boxplots.pdf", width = 10, height = 6.5)
theme.plot
dev.off()


# -- DISCUSSION: random-intercept plots by verb-----------
# ==> use results from verb to plot random intercept adjustments

intercept_df <- model_list5 # with FileID to include IrE as well

# get standard deviations per random intercept for Speaker
sds.CAN <- round(attributes(VarCorr(intercept_df$CAN.glmm)$Verb)$stddev[[1]],2)
sds.GB <- round(attributes(VarCorr(intercept_df$GB.glmm)$Verb)$stddev[[1]],2)
sds.HK <- round(attributes(VarCorr(intercept_df$HK.glmm)$Verb)$stddev[[1]],2)
sds.IND <- round(attributes(VarCorr(intercept_df$IND.glmm)$Verb)$stddev[[1]],2)
sds.IRE <- round(attributes(VarCorr(intercept_df$IRE.glmm)$Verb)$stddev[[1]],2)
sds.JA <- round(attributes(VarCorr(intercept_df$JA.glmm)$Verb)$stddev[[1]],2)
sds.NZ <- round(attributes(VarCorr(intercept_df$NZ.glmm)$Verb)$stddev[[1]],2)
sds.PHI <- round(attributes(VarCorr(intercept_df$PHI.glmm)$Verb)$stddev[[1]],2)
sds.SIN <- round(attributes(VarCorr(intercept_df$SIN.glmm)$Verb)$stddev[[1]],2)

sds.all.verb <- c(sds.CAN, sds.GB, sds.HK, sds.IND, sds.IRE, sds.JA, sds.NZ, sds.PHI, sds.SIN)


r.CAN <- ranef(intercept_df$CAN.glmm)[["Verb"]] # get by-th adjustements for CAN
r.CAN <- as.data.frame(cbind(Variety = rep("CAN", nrow(r.CAN)), r.CAN))
r.GB <- ranef(intercept_df$GB.glmm)[["Verb"]] # get by-th adjustements for GB
r.GB <- as.data.frame(cbind(Variety = rep("GB", nrow(r.GB)), r.GB))
r.HK <- ranef(intercept_df$HK.glmm)[["Verb"]] # get by-th adjustements for HK
r.HK <- as.data.frame(cbind(Variety = rep("HK", nrow(r.HK)), r.HK))
r.IND <- ranef(intercept_df$IND.glmm)[["Verb"]] # get by-th adjustements for IND
r.IND <- as.data.frame(cbind(Variety = rep("IND", nrow(r.IND)), r.IND))
r.IRE <- ranef(intercept_df$IRE.glmm)[["Verb"]] # get by-th adjustements for IRE
r.IRE <- as.data.frame(cbind(Variety = rep("IRE", nrow(r.IRE)), r.IRE))
r.JA <- ranef(intercept_df$JA.glmm)[["Verb"]] # get by-th adjustements for JA
r.JA <- as.data.frame(cbind(Variety = rep("JA", nrow(r.JA)), r.JA))
r.NZ <- ranef(intercept_df$NZ.glmm)[["Verb"]] # get by-th adjustements for NZ
r.NZ <- as.data.frame(cbind(Variety = rep("NZ", nrow(r.NZ)), r.NZ))
r.PHI <- ranef(intercept_df$PHI.glmm)[["Verb"]] # get by-th adjustements for PHI
r.PHI <- as.data.frame(cbind(Variety = rep("PHI", nrow(r.PHI)), r.PHI))
r.SIN <- ranef(intercept_df$SIN.glmm)[["Verb"]] # get by-th adjustements for SIN
r.SIN <- as.data.frame(cbind(Variety = rep("SIN", nrow(r.SIN)), r.SIN))

r.all.verb <- rbind(r.CAN, r.GB, r.HK, r.IND, r.IRE, r.JA, r.NZ, r.PHI, r.SIN)
names(r.all.verb)
names(r.all.verb)[2] <- "Intercept"
str(r.all.verb)

varieties <- levels(r.all.verb$Variety)

#--.. Plot the BLUBs ----------

(verb.plot <- ggplot(r.all.verb, aes(x = Variety, y = Intercept)) +
   geom_point(position= position_jitter(width = .15, height = 0), color="gray50", alpha=.5) +
   geom_boxplot(width = .5, alpha = .7) +
   theme(legend.position = "none") +
   annotate("text", x = varieties, y = rep(8, 9),
            label = paste("SD =", sds.all.th),
            col = 'black', 
            size=5) +
   labs(x = "", y = "by-verb adjustment", title="Random intercept by verb") +
   theme(axis.text.x = element_text(size=18), axis.title.y = element_text(size=15), plot.title = element_text(size=16, face="bold")))

pdf("Rplot_by_verb_intercept_boxplots.pdf", width = 10, height = 6.5)
verb.plot
dev.off()


library(gridExtra)

pdf("Rplot_by_lexical_all_intercept_boxplots.pdf", width = 10, height = 15)
grid.arrange(arrangeGrob(recipient.plot, theme.plot, verb.plot, ncol=1))
dev.off()




save.image("datives_probabilistic_distance.RData")