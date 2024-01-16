###############################################################################
# Load libraries
###############################################################################

library(depmixS4)
library(dplyr)
library(tidyr)
library(fitdistrplus)
library(MASS)
library(caret)
library(pROC)

###############################################################################
# Load the input data
###############################################################################

setwd("C:\\Users\\liber\\Desktop\\Study\\LMU\\Thesis Project - MariaCT's Lab\\Data\\HMM inputs\\Alleloscope_batch")

theta_vals <- read.csv("high_content_cells.csv", row.names = 1)
theta_vals$cnv <- NULL

cnv_labels <- read.csv("cnv_labels_comparison.csv", row.names = 2)
cnv_labels$X <- NULL

counts <- read.csv("read_counts_filtered.csv", row.names = 1)

#reorder theta_vals by regions
theta_vals <- theta_vals[rownames(cnv_labels),]

###############################################################################
# Prepare the data for analysis
###############################################################################

#find the best cell for a test run
foo <- c()
for (i in 1:ncol(theta_vals)){
  foo <- append(foo, sum(!is.na(theta_vals[,i])))
}

bar <- c()
for (i in 1:ncol(counts)){
  bar <- append(bar, sum(!is.na(counts[,i])))
}

best_cell <- colnames(theta_vals)[which(foo == max(foo))]

#combine data in one dataframe
dataset <- as.data.frame(matrix(NA, nrow(cnv_labels), ncol = 0))
rownames(dataset) <- rownames(cnv_labels)
dataset$counts <- counts[, best_cell]
dataset$theta_hat <- theta_vals[, best_cell]
dataset$wgs_cnv_score <- cnv_labels$wgs_score

#get a dataset with multiple cells
#transform counts data
counts_long <- counts %>% mutate(region = rownames(.))
counts_long <- counts_long %>% 
  pivot_longer(cols = -region, names_to = "cell", values_to = "counts")

#add regions as a column for further merging
cnv_labels$region <- rownames(cnv_labels)

#transform theta_hat data
theta_long <- theta_vals %>% mutate(region = rownames(.))
theta_long <- theta_long %>% 
  pivot_longer(cols = -region, names_to = "cell", values_to = "theta_hat")

#combine the data into one dataframe
all_cells_data <- full_join(counts_long, theta_long, by = c("region", "cell"))
all_cells_data <- full_join(all_cells_data, cnv_labels, by = "region")

###############################################################################
# Explorative data analysis
###############################################################################
#poisson fit test
fit_pois <- fitdist(dataset$counts, "pois")
plotdist(dataset$counts, "pois", para = list(lambda = fit_pois$estimate))
gofstat(fit_pois)

#multinomial distribution
fit_negbin <- fitdist(dataset$counts, "nbinom", 
                      start = list(size = 1, mu = mean(dataset$counts)))
plotdist(dataset$counts, "nbinom", 
         para = list(size = fit_negbin$estimate[1], mu = fit_negbin$estimate[2]))
gofstat(fit_negbin)

#log gaussian
log_data <- log(dataset$counts)
fit_normal <- fitdist(log_data, "norm")
plotdist(log_data, "norm", para = list(mean = fit_normal$estimate["mean"],
                                       sd = fit_normal$estimate["sd"]))
foo <- gofstat(fit_normal)
ks_result <- ks.test(log_data, "pnorm", mean = fit_normal$estimate["mean"], 
                     sd = fit_normal$estimate["sd"])

#class-by-class tests
dataset_base <- dataset[which(dataset$wgs_cnv_score<2.5),]
dataset_gain <- dataset[which(dataset$wgs_cnv_score>=2.5),]

#poisson
fit_pois <- fitdist(dataset_base$counts, "pois")
plotdist(dataset_base$counts, "pois", para = list(lambda = fit_pois$estimate))
gofstat(fit_pois)

fit_pois <- fitdist(dataset_gain$counts, "pois")
plotdist(dataset_gain$counts, "pois", para = list(lambda = fit_pois$estimate))
gofstat(fit_pois)

#neg.binom.
fit_negbin <- fitdist(dataset_base$counts, "nbinom", 
                      start = list(size = 1, mu = mean(dataset_base$counts)))
plotdist(dataset_base$counts, "nbinom", 
         para = list(size = fit_negbin$estimate[1], mu = fit_negbin$estimate[2]))
gofstat(fit_negbin)

fit_negbin <- fitdist(dataset_gain$counts, "nbinom", 
                      start = list(size = 1, mu = mean(dataset_gain$counts)))
plotdist(dataset_gain$counts, "nbinom", 
         para = list(size = fit_negbin$estimate[1], mu = fit_negbin$estimate[2]))
gofstat(fit_negbin)
###############################################################################
# Construct an HMM
###############################################################################

set.seed(2024)

#build a simple univariate model
# train_mod <- depmix(response = counts ~ 1, data = dataset, nstates = 2,
#                      family = gaussian("log"))
# fm <- fit(train_mod, verbose = TRUE, emc = em.control(rand = FALSE))
# print(fm)
# summary(fm)

#predict states for new data and estimate accuracy
# cells <- unique(all_cells_data$cell)
# cells <- cells[-which(cells == best_cell)]
# set.seed(2024)
# test_set <- all_cells_data[which(all_cells_data$cell == sample(cells, 1)),]
# set.seed(2024)
# test_mod <- depmix(response = counts ~ 1, data = test_set, nstates = 2,
#                    family = poisson())
# test_mod <- setpars(test_mod, getpars(f_train))

#create a dataframe for metrics
metrics <- data.frame(matrix(ncol = 9, nrow = length(unique(all_cells_data$cell))))
colnames(metrics) <- c("Accuracy", "Accuracy P-value", "Kappa", "Sensitivity", 
                       "Specificity", "Precision", "AUC", "region_count", "NA_theta")
rownames(metrics) <- unique(all_cells_data$cell)

#build a multivariate model for each cell
for (cell in unique(all_cells_data$cell)){
  set.seed(2024)
  data_sub <- all_cells_data[all_cells_data$cell == cell,]
  
  #record library size and missing theta-hat count
  metrics[cell, "region_count"] <- nrow(data_sub)
  metrics[cell, "NA_theta"] <- sum(is.na(data_sub$theta_hat))
  
  mod <- depmix(response = list(counts ~ 1, theta_hat ~ 1), data = data_sub, 
                       nstates = 2, family = list(gaussian("log"), gaussian("identity")))
  fm <- fit(mod, verbose = TRUE)
  # print(fm)
  # summary(fm)
  
  #determine which class is which CNV
  response_mtx <- summary(fm, which = "response")
  base_class <- ifelse(response_mtx[1, "Re1.(Intercept)"] < response_mtx[2, "Re1.(Intercept)"],
                       1,2)
  gain_class <- ifelse(base_class == 1, 2, 1)
  
  #create prediction and reference vectors
  true_classes <- ifelse(data_sub$wgs_score<2.5, base_class, gain_class)
  pred <- posterior(fm)[,1]
  
  #calculate performance metrics
  #create a confusion matrix and extract the metrics
  confmat <- confusionMatrix(data = as.factor(pred), reference <- as.factor(true_classes))
  
  metrics[cell, "Accuracy"] <- confmat$overall["Accuracy"]
  metrics[cell, "Accuracy P-value"] <- confmat$overall["AccuracyPValue"]
  metrics[cell, "Kappa"] <- confmat$overall["Kappa"]
  metrics[cell, "Sensitivity"] <- confmat$byClass["Sensitivity"]
  metrics[cell, "Specificity"] <- confmat$byClass["Specificity"]
  metrics[cell, "Precision"] <- confmat$byClass["Precision"]
  
  #build a ROC curve
  roc_curve <- roc(true_classes, pred)
  metrics[cell, "AUC"] <- auc(roc_curve)[1]
}

#get metric summaries
for (metric in colnames(metrics)){
  print(metric)
  print(summary(metrics[,metric]))
}

#get metric distributions
for (metric in colnames(metrics)){
  hist(metrics[,metric], main = metric, xlab = metric)
}

#evaluate epiAneufinder performance for comparison
epiAneufinder_preds <- ifelse(cnv_labels$cnv == 2, base_class, gain_class)
confmat_epiA <- confusionMatrix(data = as.factor(epiAneufinder_preds), reference <- as.factor(true_classes))
print(confmat_epiA$overall["Accuracy"])
print(confmat_epiA$overall["AccuracyPValue"])
print(confmat_epiA$overall["Kappa"])
print(confmat_epiA$byClass["Sensitivity"])
print(confmat_epiA$byClass["Specificity"])
print(confmat_epiA$byClass["Precision"])

roc_curve <- roc(true_classes, epiAneufinder_preds)
print("AUC")
print(auc(roc_curve)[1])

#check for correlation between metrics and the amount of missing theta hat values
correlations <- cor(metrics)[, "NA_theta"]
print(correlations)

###############################################################################
# Test accuracy
###############################################################################

true_classes <- ifelse(dataset$wgs_cnv_score<2.5, 1, 2)
pred <- posterior(fm)[,1]
roc_curve <- roc(true_classes, pred)
plot(roc_curve)
auc(roc_curve)

confmat <- confusionMatrix(data = as.factor(pred), reference <- as.factor(true_classes))
confmat