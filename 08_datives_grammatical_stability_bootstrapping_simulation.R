#=======================================================================
#   Script: datives_grammatical_stability_bootstrapping_simulation.R
#   Date last changed: 26/10/17
#
#   Rversion: 3.3.2
#   Author: MR
#   This script intends to assess the concept reliability of the methods
#   in comparative sociolinguistics (three lines of evidence) by using
#   bootstrapping on a randomly selected subsample of the data.
#=======================================================================


#----SETUP: LIBRARIES-----------------------------------------
library(reshape2);library(plyr);library(dplyr);library(magrittr);library(lme4);library(ggplot2);library(gridExtra); library(JGmisc); library(beepr); library(effects); library(Hmisc); library(data.table); library(MASS); library(vegan); library(cluster)

setwd("/Users/melanierothlisberger/Dropbox/Universitaet/04_Statistics/statistical analyses/20170209_THESIS_full_model")
source("/Users/melanierothlisberger/Dropbox/Universitaet/03_Programming/R/helper.R")
load("datives_probabilistic_distance.RData") # image file from 06_datives_stability_probabilistic_grammars.R



#---SETUP: PLOTTING----------------
theme_mr = theme_set(theme_light())
theme_mr = theme_update(axis.text = element_text(size = rel(1.0), color="black"),
                        axis.ticks = element_line(colour = "grey90", size = 0.25),
                        axis.title = element_text(size=rel(0.9)),
                        panel.border = element_rect(color = "black"),
                        strip.background=element_rect(fill="grey95", color="black"), 
                        strip.text.x=element_text(color="black"))


# -- Settings bootstrapping --------------
f_rand <- Resp ~ (1|Verb) + (1|ThemeHeadFilter) + (1|RecHeadFilter) + (1|FileID)
f_five <- update.formula(f_rand, .~. +
                           z.logWeightRatio +
                           c.RecPron +
                           c.ThemeBinComplexity +
                           c.ThemePron +
                           z.ThemeHeadFreq)
f_glm <- Resp ~ z.logWeightRatio +  c.RecPron + c.ThemeBinComplexity + c.ThemePron + z.ThemeHeadFreq
f_crf <- Resp ~ z.logWeightRatio +  RecPron + ThemeBinComplexity + ThemePron + z.ThemeHeadFreq

# Simulation settings
n.runs <- 1000

# --- compute the models -----------------------
# take only 50% of the data when running the model and run the model 1000 times

vars <- levels(data$Variety)

# create empty matrixes in which to save the runs per variety model
runs.sign <- matrix(nrow=length(vars),
                ncol=6,
                dimnames = list(paste0(vars),
                                rownames(summary(model_list5$CAN.glmm)$coefficients)))

runs.eff <- matrix(nrow=length(vars),
                    ncol=6,
                    dimnames = list(paste0(vars),
                                    rownames(summary(model_list5$CAN.glmm)$coefficients)))

runs.rank <- matrix(nrow=length(vars),
                    ncol=5,
                    dimnames = list(paste0(vars), c("z.logWeightRatio", "RecPron", "ThemeBinComplexity", "ThemePron", "z.ThemeHeadFreq")))



model_vars <- list()
model_list_boot <- list()
boots_sign <- list() # to save sign values of models by variety
boots_eff <- list() # to save effects of models by variety

t1 <- proc.time()
for (j in 1:n.runs){
  print(paste("Run", j))
    
  for (i in 1:9){
    print(paste("Working on", vars[i], "model..."))
    d <- droplevels(subset(data, data$Variety == vars[i]))
    half <- d[sample(nrow(d), replace=F, nrow(d)/2), ] %>% droplevels # take rand subsample
    
    # rescale predictors
    half$z.logWeightRatio = z.(half$logWeightRatio, factor=2)
    half$c.RecPron = c.(as.numeric(half$RecPron))
    half$c.ThemeBinComplexity = c.(as.numeric(half$ThemeBinComplexity))
    half$c.ThemePron = c.(as.numeric(half$ThemePron))
    half$z.ThemeHeadFreq = z.(half$ThemeHeadFreq, factor=2)
    
    m <- glm(f_glm, data=half, family="binomial")
    model_list_boot[[i]] <- m
    
    runs.sign[i,] <- as.vector(summary(m)$coefficients[,4])
    runs.eff[i,] <- as.vector(summary(m)$coefficients[,1])
    
    
  }
  
  model_vars[[j]] <- model_list_boot
  boots_sign[[j]] <- runs.sign
  boots_eff[[j]] <- runs.eff
  
  rm("d") # remove var specific df (optional)
}
t2 <- proc.time() - t1



# -- 1st line of evidence: stat sign -----------------------

boots_sign[[1]] # check
bootruns_sign <- data.frame()

for (i in 1:n.runs){
  print(paste("Run", i))
significance_matrix <- lapply(boots_sign[i], function(x) { 
  matrix(nrow = length(rownames(x)),
         ncol = length(colnames(x)),
         dimnames = list(rownames(x),
                         colnames(x)),
         data = ifelse(x < .05, TRUE, FALSE)
  )})

significance_dist <- dist(significance_matrix[[1]], method = "euclidean")
significance_dist <- significance_dist^2 # squared Euclidean distances

# Putting distances in data table (long format)
dt_dist_sig <- data.table(t(combn(rownames(significance_matrix[[1]]), 2)), 
                          as.numeric(significance_dist))
names(dt_dist_sig) <- c("var1", "var2", "distance")

# Doubling dt_dist_sig in order to compute mean distances
dt_dist_sig2 <- rbind(dt_dist_sig, dt_dist_sig[, .(var1 = var2, var2 = var1, 
                                                   distance)])


# Computing scaled version (between 0 and 1), 1 being the most similar
n_comparisons <- nrow(significance_matrix[[1]])
dt_dist_sig2[, similarity := n_comparisons - distance][]
dt_dist_sig2[, sim_scaled := similarity / n_comparisons][]

# Computing mean similarities
dt_sim_sig2_means <- dt_dist_sig2[, .(Mean_Similarity = round(mean(sim_scaled), 3)), 
                                  by = (Variety = var1)][
                                    order(Mean_Similarity, decreasing = T)]

# Finally, computing single coefficient and adding to vector

bootruns_sign <- rbind(bootruns_sign, round(dt_dist_sig2[, mean(sim_scaled)], 3))

}

bootruns_sign$Run <- rownames(bootruns_sign)
names(bootruns_sign)[1] <- "value"

# -- 2nd line of evidence: effect sizes ------------------------------

boots_eff[1]
bootruns_eff <- data.frame()

for (i in 1:n.runs){
  print(paste("Run", i))
  param_matrix <- boots_eff[[i]][, -1] # Ignoring the intercept

  # Add Null model and compute Euclidean distances
  param_matrix_null <- rbind(param_matrix, rep(0, ncol(param_matrix)))
  rownames(param_matrix_null)[length(vars) + 1] <- "NULL"
  (param_strength_null <- dist(param_matrix_null, method = "euclidean"))

  # create data table for ease of calculations
  dt_dist_strength_null <- data.table(t(combn(rownames(param_matrix_null), 2)), 
                                    as.numeric(param_strength_null))
  names(dt_dist_strength_null) <- c("var1", "var2", "distance")

  # Doubling param_strength_null in order to compute mean distances
  dt_dist_strength_null2 <- rbind(dt_dist_strength_null, dt_dist_strength_null[, .(var1 = var2, var2 = var1, distance)])
  dt_dist_strength_null_means <- dt_dist_strength_null2[, .(Mean_Dist = round(mean(distance), 3)), 
                                                      by = (Variety = var1)][
                                                        order(Mean_Dist, decreasing = T)]

  # The transformation to similarities can now be done by normalizing by 
  # distances to the null variety
  dt_dist_strength2_means[, Similarity := 1 - Mean_Dist / max(dt_dist_strength_null_means$Mean_Dist)]
  dt_dist_strength2_means[order(Similarity)]
  dt_dist_strength2_means$Variety <- factor(dt_dist_strength2_means$Variety,
                                          levels = dt_dist_strength2_means[
                                            order(Similarity, decreasing = T), 
                                            Variety])


  # Finally, computing single coefficient
  bootruns_eff <- rbind(bootruns_eff, dt_dist_strength2_means[, round(mean(Similarity), 3)])
}

bootruns_eff$Run <- rownames(bootruns_eff)
names(bootruns_eff)[1] <- "value"


# -- 3rd line of evidence (AWS server) ----------------------------
# This takes ages with the party package but is pretty fast with the 
# RandomForest package. Run on external supercluster.
# Fitting the random forests
set.seed(1245)
forest.controls = cforest_unbiased(ntree=500, mtry=3)


forest_list_boot <- list()
varimp_list_boot <- list()
forest_vars <- list()
varimp_vars <- list()
boots_rank <- list()

t1 <- proc.time()
for (j in 1:n.runs){
  print(paste("Run", j))
  
  for (i in 1:9){
    print(paste("Working on", vars[i], "model..."))
    d <- droplevels(subset(data, data$Variety == vars[i]))
    half <- droplevels(d[sample(nrow(d), replace=F, nrow(d)/2), ]) # take rand subsample
    
    # rescale predictors
    half$z.logWeightRatio = z.(half$logWeightRatio, factor=2)
    half$z.ThemeHeadFreq = z.(half$ThemeHeadFreq, factor=2)
    
    rf <- cforest(f_crf, data = half, controls = forest.controls)
    varimp <- party::varimpAUC(rf)
    varimp_list_boot[[i]] <- varimp
    forest_list_boot[[i]] <- rf
    
    runs.rank[i,] <- as.vector(rank(-varimp, ties.method="first"))
  }
  names(varimp_list_boot) <- vars
  names(forest_list_boot) <- vars
  forest_vars[[j]] <- forest_list_boot
  varimp_vars[[j]] <- varimp_list_boot
  boots_rank[[j]] <- runs.rank
  
  rm("d") # remove var specific df (optional)
}
t2 <- proc.time() - t1


# -- .. plotting ----------------------------
bootruns_rank1 <- data.frame()

for (j in 1:n.runs){
  print(paste("Run", j))
  d <- as.data.frame(boots_rank[j])
  d$Variety <- as.factor(rownames(d))
  pairs <- combn(as.factor(levels(d$Variety)), 2)
  
  # Setting up a data frame that holds pairwise correlation coefficients
  pairwise_correlations <- data.frame(
    matrix(nrow = ncol(pairs), ncol = 4, 
           dimnames = list(NULL, c("VAR1", "VAR2", "COR", "P"))))
  
  # Filling the data frame
  for(i in 1:ncol(pairs)) { # For all variety combinations ...
    var1 <- pairs[, i][1]
    var2 <- pairs[, i][2]
    
    # convert string of ranks into vector for each variety comb
    var1_values <- as.vector(as.numeric(d[d$Variety==var1,][1:5]))
    var2_values <- as.vector(as.numeric(d[d$Variety==var2,][1:5]))
    
    # Determine correlation
    correlation <- cor.test(var1_values, var2_values, method = "spearman")
    
    # Saving values
    pairwise_correlations[i, ] <- c(var1, var2, correlation$estimate, 
                                    correlation$p.value)
  }
  pairwise_correlations <- data.table(pairwise_correlations)
  
  bootruns_rank1 <- rbind(bootruns_rank1, pairwise_correlations[, mean(as.numeric(COR))])
  
}

bootruns_rank1$Run <- rownames(bootruns_rank1)
names(bootruns_rank1)[1] <- "value"
bootruns_rank1$variable <- rep("ranking", nrow(bootruns_rank1))

min.mean.sd.max <- function(x) {
  r <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

(ggplot(bootruns_rank1, aes(x="", y=value))
  + stat_summary(fun.data = min.mean.sd.max, geom="boxplot") 
  + labs(x = "constraint hierarchy")
)


# -- Plotting combined ---------------------------

bootruns_eff$variable <- "effects_value"
bootruns_sign$variable <- "sign_value"
combined <- rbind(bootruns_eff, bootruns_sign, bootruns_rank)

levels(combined$variable)
combined$variable <- as.factor(combined$variable)
combined$variable <- factor(combined$variable, levels=c("effects_value", "sign_value", "ranking"))



(boots.plot <- ggplot(combined, aes(x=factor(variable), y=value))
  + stat_summary(fun.data = min.mean.sd.max, geom="boxplot")
  + ggtitle("Boxplots with mean, 95% CI and min, max values")
  + theme(plot.title = element_text(hjust=.5))
  + labs(x="\n 3 lines of evidence", y="stability score")
  + scale_x_discrete(labels = c("coefficient estimates", "significance", "constraint ranking"))
  + theme(axis.text.x = element_text(size=13), axis.title = element_text(size=13), plot.title = element_text(face="bold"))
)

pdf("Rplot_bootstrapping_boxplots.pdf", width = 7, height = 5)
boots.plot
dev.off()

save.image("datives_simulation_bootstrapping.RData")
