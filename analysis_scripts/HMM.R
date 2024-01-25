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
library(ggpubr)

###############################################################################
# Define settings values
###############################################################################

#TODO: start substituting paths with variables for further conversion into a function
stat_flag <- FALSE
models_to_test <- "filtered" #options: "full", "filtered", c("full", "filtered")

###############################################################################
# Load the input data
###############################################################################

setwd("C:\\Users\\liber\\Desktop\\Study\\LMU\\Thesis Project - MariaCT's Lab\\Data\\HMM inputs\\Alleloscope_batch_fixed")

#load theta vals (full and filtered)
theta_vals_full <- read.csv("HMM_theta_full.csv", row.names = 1)
theta_vals_full$cnv <- NULL

theta_vals_filtered <- read.csv("HMM_theta_high_content.csv", row.names = 1)
theta_vals_filtered$cnv <- NULL

#load CNV labels (pseudobulk and per cell)
cnv_labels_pb <- read.csv("cnv_labels_comparison.csv", row.names = 2)
cnv_labels_pb$X <- NULL

cnv_labels_per_cell <- read.csv("epiAneufinder_per_cell_cnv.csv", row.names = 1)

#load read counts (full and filtered)
counts_full <- read.csv("HMM_read_counts_full.csv", row.names = 1)

counts_filtered <- read.csv("HMM_read_counts_filtered.csv", row.names = 1)

#reorder theta_vals by regions
theta_vals_full <- theta_vals_full[rownames(cnv_labels_pb),]
theta_vals_filtered <- theta_vals_filtered[rownames(cnv_labels_pb),]

###############################################################################
# Prepare the data for analysis
###############################################################################

#find the best cell for a test run
non_na_theta <- c()
for (i in 1:ncol(theta_vals_filtered)){
  non_na_theta <- append(non_na_theta, sum(!is.na(theta_vals_filtered[,i])))
}

best_cell <- colnames(theta_vals_filtered)[which(non_na_theta == max(non_na_theta))]

#combine data in one dataframe
dataset <- as.data.frame(matrix(NA, nrow(cnv_labels_pb), ncol = 0))
rownames(dataset) <- rownames(cnv_labels_pb)
dataset$counts <- counts_filtered[, best_cell]
dataset$theta_hat <- theta_vals_filtered[, best_cell]
dataset$wgs_cnv_score <- cnv_labels_pb$wgs_score

#get a dataset with multiple cells
#transform counts data
counts_long <- counts_filtered %>% mutate(region = rownames(.))
counts_long <- counts_long %>% 
  pivot_longer(cols = -region, names_to = "cell", values_to = "counts")

#add regions as a column for further merging
cnv_labels_pb$region <- rownames(cnv_labels_pb)

#transform theta_hat data
theta_long <- theta_vals_filtered %>% mutate(region = rownames(.))
theta_long <- theta_long %>% 
  pivot_longer(cols = -region, names_to = "cell", values_to = "theta_hat")

#transform the per-cell predictions
pred_per_cell_long <- cnv_labels_per_cell %>% mutate(region = rownames(.))
pred_per_cell_long <- pred_per_cell_long %>%
  pivot_longer(cols = -region, names_to = "cell", values_to = "cnv_per_cell")

#combine the data into one dataframe
all_cells_data_filtered <- full_join(counts_long, theta_long, by = c("region", "cell"))
all_cells_data_filtered <- full_join(all_cells_data_filtered, cnv_labels_pb, by = "region")
all_cells_data_filtered <- left_join(all_cells_data_filtered, pred_per_cell_long, by = c("region", "cell"))

#get a dataset for unfiltered data
#transform counts data
counts_trimmed <- counts_full[, colnames(counts_full) %in% colnames(theta_vals_full)]

counts_long <- counts_trimmed %>% mutate(region = rownames(.))
counts_long <- counts_long %>% 
  pivot_longer(cols = -region, names_to = "cell", values_to = "counts")

#transform theta_hat data
#TODO: figure out why we have more cells for theta vals than counts and preds
#for now, remove extra cells
theta_trimmed <- theta_vals_full[, colnames(theta_vals_full) %in% colnames(counts_trimmed)]

theta_long <- theta_trimmed %>% mutate(region = rownames(.))
theta_long <- theta_long %>% 
  pivot_longer(cols = -region, names_to = "cell", values_to = "theta_hat")

#combine the data into one dataframe
all_cells_data_full <- full_join(counts_long, theta_long, by = c("region", "cell"))
all_cells_data_full <- full_join(all_cells_data_full, cnv_labels_pb, by = "region", "cell")
all_cells_data_full <- left_join(all_cells_data_full, pred_per_cell_long, by = c("region", "cell"))
###############################################################################
# Explorative data analysis
###############################################################################

if (stat_flag == TRUE){
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
  
  #log gaussian
  log_data <- log(dataset_base$counts)
  fit_normal <- fitdist(log_data, "norm")
  plotdist(log_data, "norm", para = list(mean = fit_normal$estimate["mean"],
                                         sd = fit_normal$estimate["sd"]))
  foo <- gofstat(fit_normal)
  ks_result <- ks.test(log_data, "pnorm", mean = fit_normal$estimate["mean"], 
                       sd = fit_normal$estimate["sd"])
  
  log_data <- log(dataset_gain$counts)
  fit_normal <- fitdist(log_data, "norm")
  plotdist(log_data, "norm", para = list(mean = fit_normal$estimate["mean"],
                                         sd = fit_normal$estimate["sd"]))
  foo <- gofstat(fit_normal)
  ks_result <- ks.test(log_data, "pnorm", mean = fit_normal$estimate["mean"], 
                       sd = fit_normal$estimate["sd"])
}

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

#create a list for storage of models
models <- replicate(length(models_to_test), list(), simplify = FALSE)
names(models) <- models_to_test
if ("filtered" %in% models_to_test){
  models$filtered[["data"]] <- all_cells_data_filtered
}
if ("full" %in% models_to_test){
  models$full[["data"]] <- all_cells_data_full
}

#train HMM models and calculate performance metrics
for (model in names(models)){
  #create a list for models
  models[[model]][["fitted_models"]] <- list()
  #create a dataframe for metrics
  metrics <- data.frame(matrix(ncol = 9, nrow = length(unique(models[[model]][["data"]]$cell))))
  colnames(metrics) <- c("Accuracy", "Accuracy P-value", "Kappa", "Sensitivity",
                         "Specificity", "Precision", "AUC", "region_count", "NA_theta")
  rownames(metrics) <- unique(models[[model]][["data"]]$cell)
  models[[model]][["metrics"]] <- metrics
  
  #create a dataframe for predictions
  HMM_preds <- data.frame(fragment = character(0), cell = character(0), HMM_pred = numeric(0), WGS_score = numeric(0))
  
  #build a multivariate model for each cell
  for (cell in unique(models[[model]][["data"]]$cell)){
    set.seed(2024)
    data_sub <- models[[model]][["data"]][models[[model]][["data"]]$cell == cell,]
    
    #skip the cell if there are less than 5 non-NA theta-hat values present for it
    #TODO: figure out why different models break at different non_NA theta-hat numbers
    if (sum(!is.na(data_sub$theta_hat))>=5){
      #record library size and missing theta-hat count
      metrics[cell, "region_count"] <- nrow(data_sub)
      metrics[cell, "NA_theta"] <- sum(is.na(data_sub$theta_hat))
      
      mod <- depmix(response = list(counts ~ 1, theta_hat ~ 1), data = data_sub, 
                    nstates = 2, family = list(gaussian("log"), gaussian("identity")))
      print(paste("cell:", cell, "model", model))
      fm <- fit(mod, verbose = FALSE)
      models[[model]][["fitted_models"]][[cell]] <- fm
      # print(fm)
      # summary(fm)
      
      #determine which class is which CNV
      response_mtx <- summary(fm, which = "response")
      #TODO: when we have 3 classes, change the label assigment to include losses and do it through min/max functions
      #loss_class <- which(response_mtx[,"Re1.(Intercept)"] == min(response_mtx[,"Re1.(Intercept)"]))
      base_class <- ifelse(response_mtx[1, "Re1.(Intercept)"] < response_mtx[2, "Re1.(Intercept)"],
                           1,2)
      gain_class <- ifelse(base_class == 1, 2, 1)
      
      #create prediction and reference vectors
      #TODO: add class assignment for 3 classes when we have losses in the data
      true_classes <- ifelse(data_sub$wgs_score<2.5, base_class, gain_class)
      pred <- posterior(fm)[,1]
      
      #map predictions to the HMM_preds
      insert_df <- data_sub[, c("region", "cell")]
      insert_df$HMM_pred <- ifelse(pred == base_class, 2,3)
      insert_df$WGS_score <- ifelse(data_sub$wgs_score<2.5, 2, 3)
      HMM_preds <- rbind(HMM_preds, insert_df)
      
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
  }
  models[[model]][["predictions"]] <- HMM_preds
  models[[model]][["metrics"]] <- metrics
  
  #get metric summaries
  models[[model]][["metrics_summary"]] <- summary(metrics)
}

# #build a multivariate model for each cell
# for (cell in unique(all_cells_data$cell)){
#   set.seed(2024)
#   data_sub <- all_cells_data[all_cells_data$cell == cell,]
#   
#   #record library size and missing theta-hat count
#   metrics[cell, "region_count"] <- nrow(data_sub)
#   metrics[cell, "NA_theta"] <- sum(is.na(data_sub$theta_hat))
#   
#   mod <- depmix(response = list(counts ~ 1, theta_hat ~ 1), data = data_sub, 
#                        nstates = 2, family = list(gaussian("log"), gaussian("identity")))
#   fm <- fit(mod, verbose = TRUE)
#   # print(fm)
#   # summary(fm)
#   
#   #determine which class is which CNV
#   response_mtx <- summary(fm, which = "response")
#   base_class <- ifelse(response_mtx[1, "Re1.(Intercept)"] < response_mtx[2, "Re1.(Intercept)"],
#                        1,2)
#   gain_class <- ifelse(base_class == 1, 2, 1)
#   
#   #create prediction and reference vectors
#   true_classes <- ifelse(data_sub$wgs_score<2.5, base_class, gain_class)
#   pred <- posterior(fm)[,1]
#   
#   #map predictions to the HMM_preds
#   insert_df <- data_sub[, c("region", "cell")]
#   insert_df$HMM_pred <- ifelse(pred == base_class, 2,3)
#   insert_df$WGS_score <- ifelse(data_sub$wgs_score<2.5, 2, 3)
#   HMM_preds <- rbind(HMM_preds, insert_df)
#   
#   #calculate performance metrics
#   #create a confusion matrix and extract the metrics
#   confmat <- confusionMatrix(data = as.factor(pred), reference <- as.factor(true_classes))
#   
#   metrics[cell, "Accuracy"] <- confmat$overall["Accuracy"]
#   metrics[cell, "Accuracy P-value"] <- confmat$overall["AccuracyPValue"]
#   metrics[cell, "Kappa"] <- confmat$overall["Kappa"]
#   metrics[cell, "Sensitivity"] <- confmat$byClass["Sensitivity"]
#   metrics[cell, "Specificity"] <- confmat$byClass["Specificity"]
#   metrics[cell, "Precision"] <- confmat$byClass["Precision"]
#   
#   #build a ROC curve
#   roc_curve <- roc(true_classes, pred)
#   metrics[cell, "AUC"] <- auc(roc_curve)[1]
# }
# write.csv(metrics, "..//..//HMM outputs//Alleloscope_old_epiA_batch//metrics.csv")

#get metric summaries
# summary(metrics)
# write.csv(summary(metrics), "..//..//HMM outputs//Alleloscope_old_epiA_batch//metrics_summary.csv")

###############################################################################
# Test accuracy
###############################################################################

#get metric distributions
for (metric in names(metrics)){
  plot_data <- data.frame()
  for (model in names(models)){
    plot_data[model] <- models[[model]][["metrics"]]$metric
  }
}

#get metric distributions
pdf(file <- "..//..//HMM outputs//Alleloscope_old_epiA_batch//metrics_distributions.pdf")
for (metric in colnames(metrics)){
  hist(metrics[,metric], main = metric, xlab = metric)
}
dev.off()

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
# Plot the results on a karyogram
###############################################################################

#plot CNV predictions
#extract chromosome number and starting position
HMM_preds$chr <- gsub("chr([0-9]+):.*", "chr\\1", HMM_preds$region)
HMM_preds$fragment <- as.numeric(gsub("chr[0-9]+:([0-9]+)", "\\1", HMM_preds$region))

#make chr column a factor
HMM_preds$chr <- factor(HMM_preds$chr, levels = paste0("chr", 1:22))

#reorder the dataframe
HMM_preds <- HMM_preds[order(HMM_preds$chr, HMM_preds$fragment), ]

#add row number and record positions for labels
HMM_preds <- HMM_preds %>%
  mutate(X_pos = as.numeric(factor(region, levels = unique(region))))

chromosome_label_positions <- HMM_preds %>%
  filter(!duplicated(chr)) %>%
  pull(X_pos)

#create a plot for the predictions
p <- ggplot(HMM_preds, aes(x = X_pos, y = cell, fill = factor(HMM_pred))) +
  geom_tile() +
  geom_vline(xintercept = chromosome_label_positions - 0.5, color = "black", linetype = "dotted", size = 0.5) +
  labs(title = "CNV predictions by the HMM",
       x = "Region",
       y = "Cell",
       fill = "CNV class") +
  scale_x_continuous(breaks = chromosome_label_positions, labels = unique(HMM_preds$chr),
                     expand = c(0.05, 0)) +
  scale_fill_manual(values = c("#1fb424", "#b4421f"),
                    labels = c("base", "gain")) +
  theme_minimal() +
  theme(legend.position = "top", 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_blank()) 

#create a plot for the ground truth
p_wgs <- ggplot(HMM_preds, aes(x = X_pos, y = 1, fill = factor(WGS_score))) +
  geom_tile() +
  geom_vline(xintercept = chromosome_label_positions - 0.5, color = "black", linetype = "dotted", size = 0.5) +
  labs(title = "CNV predictions by the HMM",
       x = "Region",
       y = "WGS",
       fill = "CNV class") +
  scale_x_continuous(breaks = chromosome_label_positions, labels = unique(HMM_preds$chr),
                     expand = c(0.05, 0)) +
  scale_fill_manual(values = c("#1fb424", "#b4421f"),
                    labels = c("base", "gain")) +
  theme_minimal() +
  theme(legend.position = "top", 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_blank())

#combine the plots
combined_plot <- ggarrange(p, p_wgs + labs(title = NULL) + theme(legend.position="none", ), 
                           ncol = 1, heights = c(88, 12))

ggsave("..//..//HMM outputs//Alleloscope_old_epiA_batch//HMM_karyogram.pdf", 
       plot = combined_plot, width = 40, height = 20, units = "cm")
