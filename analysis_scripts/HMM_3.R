###############################################################################
# Load libraries
###############################################################################

library(depmixS4)
library(dplyr)
library(tidyr)
library(grid)
library(gridExtra)
library(patchwork)
library(RColorBrewer)
library(fitdistrplus)
library(MASS)
library(caret)
library(pROC)
library(ggpubr)

###############################################################################
# Load functions
###############################################################################

source("C:\\Users\\liber\\Documents\\GitHub\\lmu-thesis-project\\analysis_scripts\\by_cell_plotter.R")

###############################################################################
# Define settings values
###############################################################################

#TODO: start substituting paths with variables for further conversion into a function
seed_value <- 2024
stat_flag <- FALSE
impute_theta_flag <- TRUE
#options: "counts", "counts_and_theta", "counts_and_theta_imp" - for multiple use a vector
models_to_test <- c("counts", "counts_and_theta", "counts_and_theta_imp")

#adapt for the two-state model
#TODO: implement it in the code
cnv_state_vector <- c("loss", "base", "gain")

input_dir <- "C:\\Users\\liber\\Desktop\\Study\\LMU\\Thesis Project - MariaCT's Lab\\Data\\HMM inputs\\Alleloscope_1mb_100nSNP"
output_dir <- paste("..//..//HMM outputs//Alleloscope_1mb_100nSNP//", "3_state_model//", sep = "")

dir.create(output_dir)

###############################################################################
# Load the input data
###############################################################################

setwd(input_dir)

#load theta vals (full and filtered)
theta_vals_full <- read.csv("HMM_theta_full.csv", row.names = 1)
theta_vals_full$cnv <- NULL

theta_vals_imputed <- read.csv("HMM_theta_full_imputed.csv", row.names = 1)
theta_vals_imputed$cnv <- NULL

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
theta_vals_imputed <- theta_vals_imputed[rownames(cnv_labels_pb),]
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

#get a dataset for imputed_theta_data
#change theta_vals_full to be identical to the imputed version if the impute_theta flag is enabled
if (impute_theta_flag){
  theta_vals_full <- theta_vals_imputed
  theta_vals_full <- theta_vals_full[,colnames(theta_vals_filtered)]
}
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
  # #poisson fit test
  # fit_pois <- fitdist(dataset$counts, "pois")
  # plotdist(dataset$counts, "pois", para = list(lambda = fit_pois$estimate))
  # gofstat(fit_pois)
  # 
  # #multinomial distribution
  # fit_negbin <- fitdist(dataset$counts, "nbinom", 
  #                       start = list(size = 1, mu = mean(dataset$counts)))
  # plotdist(dataset$counts, "nbinom", 
  #          para = list(size = fit_negbin$estimate[1], mu = fit_negbin$estimate[2]))
  # gofstat(fit_negbin)
  # 
  # #log gaussian
  # log_data <- log(dataset$counts)
  # fit_normal <- fitdist(log_data, "norm")
  # plotdist(log_data, "norm", para = list(mean = fit_normal$estimate["mean"],
  #                                        sd = fit_normal$estimate["sd"]))
  # foo <- gofstat(fit_normal)
  # ks_result <- ks.test(log_data, "pnorm", mean = fit_normal$estimate["mean"], 
  #                      sd = fit_normal$estimate["sd"])
  
  #class-by-class tests
  dataset_loss <- dataset[which(dataset$wgs_cnv_score<=1.5),]
  dataset_base <- dataset[which(dataset$wgs_cnv_score<2.5) && which(dataset$wgs_cnv_score>1.5),]
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

set.seed(seed_value)

#build a simple univariate model
# train_mod <- depmix(response = counts ~ 1, data = dataset, nstates = 2,
#                      family = gaussian("log"))
# fm <- fit(train_mod, verbose = TRUE, emc = em.control(rand = FALSE))
# print(fm)
# summary(fm)

#predict states for new data and estimate accuracy
# cells <- unique(all_cells_data$cell)
# cells <- cells[-which(cells == best_cell)]
# set.seed(seed_value)
# test_set <- all_cells_data[which(all_cells_data$cell == sample(cells, 1)),]
# set.seed(seed_value)
# test_mod <- depmix(response = counts ~ 1, data = test_set, nstates = 2,
#                    family = poisson())
# test_mod <- setpars(test_mod, getpars(f_train))

#create a list for storage of models
models <- replicate(length(models_to_test), list(), simplify = FALSE)
names(models) <- models_to_test
for (model in names(models)){
  if (model == "counts_and_theta_imp"){
    models[[model]][["data"]] <- all_cells_data_full
  } else {
    models[[model]][["data"]] <- all_cells_data_filtered
  }
}

#train HMM models and calculate performance metrics
for (model in names(models)){
  #create a list for models
  models[[model]][["fitted_models"]] <- list()

  #create dataframes for metrics storage
  overall_metrics_df <- data.frame(matrix(ncol = 3, nrow = length(unique(models[[model]][["data"]]$cell))))
  colnames(overall_metrics_df) <- c("overall_Acc", "NA_theta", "region_count")
  rownames(overall_metrics_df) <- unique(models[[model]][["data"]]$cell)
  
  per_class_metrics_df <- data.frame(matrix(ncol = 12, nrow = length(unique(models[[model]][["data"]]$cell))))
  colnames(per_class_metrics_df) <- c("Sensitivity", "Specificity", "Pos Pred Value",
                                      "Neg Pred Value", "Precision", "Recall", "F1",
                                      "Prevalence", "Detection Rate", "Detection Prevalence",
                                      "Balanced Accuracy", "AUC")
  rownames(per_class_metrics_df) <- unique(models[[model]][["data"]]$cell)
  
  #create a list for metrics storage
  metrics_list <- list("loss" = per_class_metrics_df, "base" = per_class_metrics_df, 
                       "gain" = per_class_metrics_df, "overall" = overall_metrics_df)
  models[[model]][["metrics"]] <- metrics_list
  
  #create a dataframe for predictions
  HMM_preds <- data.frame(fragment = character(0), cell = character(0), HMM_pred = numeric(0), WGS_score = numeric(0))
  
  #build a multivariate model for each cell
  for (cell in unique(models[[model]][["data"]]$cell)){
    set.seed(seed_value)
    data_sub <- models[[model]][["data"]][models[[model]][["data"]]$cell == cell,]
    
    #skip the cell if there are less than 5 non-NA theta-hat values present for it
    #TODO: figure out why different models break at different non_NA theta-hat numbers
    # if (sum(!is.na(data_sub$theta_hat))>=5){}
    
    #record library size and missing theta-hat count
    metrics_list[["overall"]][cell, "region_count"] <- nrow(data_sub)
    metrics_list[["overall"]][cell, "NA_theta"] <- sum(is.na(data_sub$theta_hat))
    
    #choose a model for training
    if (model == "counts"){
      mod <- depmix(response = counts ~ 1, data = data_sub, 
                    nstates = 3, family = gaussian("log"))
    } else {
      mod <- depmix(response = list(counts ~ 1, theta_hat ~ 1), data = data_sub, 
                    nstates = 3, family = list(gaussian("log"), gaussian("identity")))
    }
    print(paste("cell:", cell, "model:", model))
    fm <- fit(mod, verbose = FALSE)
    models[[model]][["fitted_models"]][[cell]] <- fm
    # print(fm)
    # summary(fm)
    
    #determine which class is which CNV
    response_mtx <- summary(fm, which = "response")
    
    loss_class <- which(response_mtx[,"Re1.(Intercept)"] == min(response_mtx[,"Re1.(Intercept)"]))
    gain_class <- which(response_mtx[,"Re1.(Intercept)"] == max(response_mtx[,"Re1.(Intercept)"]))
    base_class <- which(!(c(1,2,3) %in% c(loss_class, gain_class)))
    
    #'@two_state_relabeling
    # base_class <- ifelse(response_mtx[1, "Re1.(Intercept)"] < response_mtx[2, "Re1.(Intercept)"],
    #                      1,2)
    # gain_class <- ifelse(base_class == 1, 2, 1)
    
    #create prediction and reference vectors
    true_classes <- ifelse(data_sub$wgs_score >= 2.5, gain_class, ifelse(data_sub$wgs_score <= 1.5, loss_class, base_class))
    pred <- posterior(fm)[,1]
    posterior_df <- posterior(fm)[,c(2:4)]
    posterior_df <- posterior_df[, c(loss_class, base_class, gain_class)]
    
    #map predictions to the HMM_preds
    insert_df <- data_sub[, c("region", "cell")]
    insert_df$HMM_pred <- ifelse(pred == loss_class, 1, ifelse(pred == base_class, 2, 3))
    insert_df$WGS_score <- ifelse(data_sub$wgs_score<=1.5, 1, ifelse(data_sub$wgs_score>=2.5, 3, 2))
    HMM_preds <- rbind(HMM_preds, insert_df)
    
    #create a confusion matrix and extract the metrics
    confmat <- confusionMatrix(data = as.factor(insert_df$HMM_pred), reference <- as.factor(insert_df$WGS_score))
    for (cnv_type in c(1:3)){
      #extract the metrics into a single-row dataframe
      extracted_metrics <- confmat$byClass[cnv_type,]
      
      #build a ROC curve
      binary_ref <- ifelse(insert_df$WGS_score == cnv_type, 1, 0)
      #OLD roc curve with binary preds
      #binary_pred <- ifelse(pred == cnv_type, 1, 0)
      #roc_curve <- roc(binary_ref, binary_pred)
      roc_curve <- roc(binary_ref, posterior_df[,cnv_type])
      
      #add AUC as a column to the extracted_metrics df
      extracted_metrics$AUC <- auc(roc_curve)[1]
      
      #add the metrics as a row to the according dataframe in the model storage list
      metrics_list[[cnv_type]][cell, ] <- extracted_metrics
    }
    metrics_list[["overall"]][cell, "overall_Acc"] <- confmat$overall[["Accuracy"]]
  }
  models[[model]][["predictions"]] <- HMM_preds
  models[[model]][["metrics"]] <- metrics_list
  
  #get metric summaries
  models[[model]][["metrics_summary"]] <- list("loss" = NA, "base" = NA, "gain" = NA)
  for (cnv_type in names(models[[model]][["metrics_summary"]])){
    models[[model]][["metrics_summary"]][[cnv_type]] <- summary(models[[model]][["metrics"]][[cnv_type]])
  }
}

#leave only the relevant metrics
for (model in names(models)){
  for (cnv_class in c("loss", "base", "gain")){
    for (metric in colnames(models[[model]][["metrics"]][[cnv_class]])){
      if (all(is.nan(models[[model]][["metrics"]][[cnv_class]][,metric]))){
        models[[model]][["metrics"]][[cnv_class]][,metric] <- 0
      }
    }
  }
}

#save the models list as an RDS
saveRDS(models, paste(output_dir, "models.RDS", sep = ""))

# #build a multivariate model for each cell
# for (cell in unique(all_cells_data$cell)){
#   set.seed(seed_value)
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
# Evaluate model metrics and save their distributions as plots
###############################################################################

#check for correlation between metrics and the amount of missing theta hat values
for (model in names(models)){
  #create an output dataframe
  corr_df <- data.frame(matrix(ncol = 15, nrow = 3))
  rownames(corr_df) <- c("loss", "base", "gain")    #TODO: change if more than one model will be implemented
  
  for (cnv_type in names(models[[model]][["metrics"]])[-4]){
    print(paste("Correlation between metrics and missing theta-values for the", cnv_type, "class (", model, "model)"))
    metrics_merged <- models[[model]][["metrics"]][[cnv_type]]
    metrics_merged <- cbind(models[[model]][["metrics"]][[cnv_type]], models[[model]][["metrics"]][["overall"]])
    correlations <- cor(metrics_merged)[, "NA_theta"]
    corr_df[cnv_type,] <- correlations
    print(correlations)
  }
  #save the output dataframe
  colnames(corr_df) <- names(correlations)
  write.csv(corr_df, paste(output_dir, model, "_", "metrics_correlation_with_NA_theta.csv", sep = ""))
}


#get metric distributions
#for per_class metrics
pdf(paste(output_dir, "metrics_distributions_per_class.pdf", sep = ""))
for (metric in colnames(per_class_metrics_df)){
  for (model in names(models)){
    for (cnv_class in c("loss", "base", "gain")){
      print(paste(metric, model))
      plot_data <- models[[model]][["metrics"]][[cnv_class]][,metric]
      
      par(mar = c(5.1, 4.1, 4.1, 8.1))
      hist(plot_data, breaks = 10, main = paste(metric, "distribution for", model, 
                                                "model", paste("(", cnv_class, ")", sep = "")), 
           xlab = metric)
      summarized_data <- summary(plot_data)
      summarized_data <- paste(names(summarized_data), round(summarized_data, 3), sep = ":")
      legend("topright", legend = summarized_data, xpd = TRUE, inset = c(-0.3, 0))
    }
  }
}
dev.off()

#for overall metrics
pdf(paste(output_dir, "metrics_distributions_overall.pdf", sep = ""))
for (metric in colnames(overall_metrics_df)){
  for (model in names(models)){
    plot_data <- models[[model]][["metrics"]][["overall"]][, metric]
    
    par(mar = c(5.1, 4.1, 4.1, 8.1))
    hist(plot_data, breaks = 10, main = paste(metric, "distribution for", model, 
                                              "model"), 
         xlab = metric)
    summarized_data <- summary(plot_data)
    summarized_data <- paste(names(summarized_data), round(summarized_data, 3), sep = ":")
    legend("topright", legend = summarized_data, xpd = TRUE, inset = c(-0.3, 0))
  }
}
dev.off()

#plot model metrics against each other
#per-class metrics
pdf(paste(output_dir, "different_HMM_per_class_metrics_comparison.pdf", sep = ""))
for (metric in colnames(per_class_metrics_df)){
  for (cnv_type in c("loss", "base", "gain")){
    extracted_metrics <- list()
    for (model in names(models)){
      extracted_metrics[[model]] <- models[[model]][["metrics"]][[cnv_type]][, metric] 
    }
    value_column <- c()
    for (model in names(extracted_metrics)){
      value_column <- c(value_column, extracted_metrics[[model]])
    }
    plot_data <- data.frame(value = value_column, 
                            source_model = rep(names(extracted_metrics),
                                               each = length(extracted_metrics[[1]])))
    
    p <- ggplot(plot_data, aes(value, fill = source_model)) + 
      geom_histogram(position = "dodge", alpha = 0.5, colour = "black", bins = 15) +
      scale_fill_brewer(palette = "Set1") +
      theme_minimal() +
      labs(title = paste(metric, "distribution for the ", names(extracted_metrics)[1], 
                         "and", names(extracted_metrics)[2],"HMM models \n",
                         paste("(", cnv_type, ")", sep = "")), 
           x = metric, y = "Count")
    print(p)
  }
}
dev.off()

#overall metrics
pdf(paste(output_dir, "different_HMM_overall_metrics_comparison.pdf", sep = ""))
for (metric in colnames(overall_metrics_df)){
  extracted_metrics <- list()
  for (model in names(models)){
    extracted_metrics[[model]] <- models[[model]][["metrics"]][["overall"]][, metric] 
  }
  value_column <- c()
  for (model in names(extracted_metrics)){
    value_column <- c(value_column, extracted_metrics[[model]])
  }
  plot_data <- data.frame(value = value_column, 
                          source_model = rep(names(extracted_metrics),
                                             each = length(extracted_metrics[[1]])))
  
  p <- ggplot(plot_data, aes(value, fill = source_model)) + 
    geom_histogram(position = "dodge", alpha = 0.5, colour = "black", bins = 15) +
    scale_fill_brewer(palette = "Set1") +
    theme_minimal() +
    labs(title = paste(metric, "distribution for the ", names(extracted_metrics)[1], 
                       "and", names(extracted_metrics)[2],"HMM models"), 
         x = metric, y = "Count")
  print(p)
}
dev.off()
###############################################################################
# Test epiAneufinder accuracy
###############################################################################

#evaluate epiAneufinder performance for comparison
#create a discretized ground truth vector
epiA_truth <- ifelse(cnv_labels_pb$wgs_score <= 1.5, 1, ifelse(cnv_labels_pb$wgs_score >= 2.5, 3, 2))
confmat_epiA <- confusionMatrix(data = as.factor(cnv_labels_pb$cnv), 
                                reference = as.factor(epiA_truth))

print(confmat_epiA$overall["Accuracy"])
print(confmat_epiA$byClass)

#save results in an output dataframe
epiA_eval_df <- as.data.frame(confmat_epiA$byClass)
rownames(epiA_eval_df) <- c("loss", "base", "gain")
epiA_eval_df$overall_Acc <- confmat_epiA$overall["Accuracy"]

#build ROC curves
epiA_eval_df$AUC <- NA
for (cnv_type in c(1:3)){
  binary_ref <- ifelse(epiA_truth == cnv_type, 1, 0)
  binary_pred <- ifelse(cnv_labels_pb$cnv == cnv_type, 1, 0)
  roc_curve <- roc(binary_ref, binary_pred)
  
  #add AUC as a column to the output df
  print(paste("AUC for", c("loss", "base", "gain")[cnv_type], "class"))
  print(auc(roc_curve)[1])
  
  epiA_eval_df$AUC[cnv_type] <- auc(roc_curve)[1]
}

write.csv(epiA_eval_df, paste(output_dir, "epiAneufinder_pseudobulk_metrics.csv", sep = ""))

#calculate metrics for per-cell predictions
#subset the per-cell predictions to only include regions, present in Alleloscope output and cells with enough theta hat values
cnv_labels_per_cell_filtered <- cnv_labels_per_cell[rownames(cnv_labels_per_cell) %in% rownames(cnv_labels_pb), 
                                                    colnames(cnv_labels_per_cell) %in% all_cells_data_filtered$cell]

#create an output dataframe
epiA_eval <- data.frame(matrix(ncol = ncol(per_class_metrics_df) + 3,
                               nrow = 0))
colnames(epiA_eval) <- c(colnames(per_class_metrics_df), "overall_Acc", "cell", "cnv_class")

#create a confusion matrix for each cell and save metrics as a dataframe
for (cell in colnames(cnv_labels_per_cell_filtered)){
  epiA_pred <- cnv_labels_per_cell_filtered[, cell]
  confmat_epiA <- confusionMatrix(data = as.factor(epiA_pred),
                                  reference = as.factor(epiA_truth))
  for (cnv_class in (c("loss", "base", "gain"))){
    #extract confusion matrix metrics
    extracted_metrics <- as.data.frame(t(confmat_epiA$byClass[which(c("loss", "base", "gain") == cnv_class),]))
    extracted_metrics$overall_Acc <- confmat_epiA$overall[["Accuracy"]]
    
    #build ROC curves
    binary_ref <- ifelse(epiA_truth == which(c("loss", "base", "gain") == cnv_class), 1, 0)
    roc_curve <- roc(binary_ref, epiA_pred)
    
    #add AUC as a column to the output df
    extracted_metrics$AUC <- auc(roc_curve)[1]
    
    #save cell and cnv_type
    extracted_metrics$cell <- cell
    extracted_metrics$cnv_class <- cnv_class
    
    #save metrics to the output dataframe
    epiA_eval <- rbind(epiA_eval, extracted_metrics) 
  }
}
#save the output dataframe
write.csv(epiA_eval, paste(output_dir, "epiAneufinder_per_cell_metrics.csv", sep = ""))

#plot epiAneufinder distributions with summaries
pdf(paste(output_dir, "epiA_metrics_dists_with_summary.pdf", sep = ""))
par(mar = c(5.1, 4.1, 4.1, 8.1))
for (cnv_class in c("loss", "base", "gain")){
  plot_data <- epiA_eval[epiA_eval$cnv_class == cnv_class,]
  for (metric in c("overall_Acc", "Balanced Accuracy", "Sensitivity", 
                   "Specificity", "Precision", "Prevalence", "AUC")){
    print(paste(cnv_type, metric))
    hist(plot_data[,metric], breaks = 10, 
         main = paste(metric, "distribution for", model,
                      "model", paste("(", cnv_class, ")", sep = "")), 
         xlab = metric)
    summarized_data <- summary(plot_data[,metric])
    summarized_data <- paste(names(summarized_data), round(summarized_data, 3), sep = ":")
    legend("topright", legend = summarized_data, xpd = TRUE, inset = c(-0.3, 0))
  }
}
dev.off()

#plot epiAneufinder per-cell metrics against HMM metrics
#one model against epiAneufinder
for (model in names(models)){
  pdf(paste(output_dir, model, "_epiAneufinder_and_HMM_metrics_comparison.pdf", sep = ""))
  for (cnv_type in c("loss", "base", "gain")){
    #subset the data by cnv type and model name
    data_sub <- models[[model]][["metrics"]][[cnv_type]]
    epiA_sub <- epiA_eval[(epiA_eval$cnv_class == cnv_type),]
    #plot each metric
    for (metric in colnames(data_sub)){
      #create a dataframe for plotting
      plot_data <- data.frame(value = c(data_sub[, metric], epiA_sub[, metric]),
                              source_model = rep(c("HMM", "epiAneufinder"), each = nrow(data_sub)))
      # create a ggplot object
      p <- ggplot(plot_data, aes(value, fill = source_model)) +
        geom_histogram(position = "identity", alpha = 0.5, colour = "black", bins = 30) +
        scale_fill_brewer(palette = "Set1") +
        theme_minimal() +
        labs(title = paste(metric, "distribution for the ", model, "HMM model and epiAneufinder \n",
                           paste("(", cnv_type, ")", sep = "")),
             x = metric, y = "Count")
      print(p)
    }
  }
  dev.off()
}

#epiAneufinder against all models
#set the y axis limits
y_limits <- c(0, nrow(models[[1]][["metrics"]][[cnv_type]]))
#create a list for plot_storage
plot_list <- list()
for (cnv_type in c("loss", "base", "gain")){
  #subset the data by cnv type and model name
  extracted_data <- list()
  for (model in names(models)){
    extracted_data[[model]] <- models[[model]][["metrics"]][[cnv_type]]
    extracted_data[[model]]$overall_Acc <- models[[model]][["metrics"]][["overall"]]$overall_Acc
    extracted_data[[model]]$source_model <- model
    rownames(extracted_data[[model]]) <- c(1:nrow(extracted_data[[model]]))
  }
  extracted_data[["epiAneufinder"]] <- epiA_eval[(epiA_eval$cnv_class == cnv_type), 
                                                 !(colnames(epiA_eval) %in% c("cell", "cnv_class"))]
  extracted_data[["epiAneufinder"]]$source_model <- "epiAneufinder"
  rownames(extracted_data[["epiAneufinder"]]) <- c(1:nrow(extracted_data[["epiAneufinder"]]))
  
  #join data into a dataframe
  plot_data <- data.frame(matrix(nrow = 0, ncol = ncol(extracted_data[[1]])))
  colnames(plot_data) <- colnames(extracted_data[[1]])
  for (model in names(extracted_data)){
    plot_data <- rbind(plot_data, extracted_data[[model]])
  }
  
  #plot each metric (drop overall_Acc to avoid duplicate plots)
  for (metric in colnames(plot_data)[colnames(plot_data) != c("overall_Acc", "source_model")]){
    #if there is no sublist for the metric in the plot_list, create one
    if (!(metric %in% names(plot_list))){
      plot_list[[metric]] <- list()
    }
    
    # create a ggplot object
    p <- ggplot(plot_data, aes(.data[[metric]], fill = source_model)) + 
      geom_histogram(position = "dodge", alpha = 0.5, colour = "black", bins = 15) +
      scale_fill_brewer(palette = "Set1",
                        labels = c("Counts", "Counts and AF", "Counts and imputed AF", "epiAneufinder")) +
      ylim(y_limits) +
      theme_minimal() +
      theme(legend.position = "bottom") +
      labs(title = paste(metric, "distribution for\nHMM models and epiAneufinder",
                         paste("(", cnv_type, ")", sep = "")), 
           x = metric, y = "Number of cells", fill = "Model")
    print(p)
    
    #save the plot into the plot_list
    plot_list[[metric]][[cnv_type]] <- p
  }
}

#save plots, grouped by metric
pdf(paste(output_dir, "all_HMM_and_epiAneufinder_metrics_comparison.pdf", sep = ""),
    width = 15, height = 6)

for (metric in names(plot_list)){
  #join plots for the same metric on the same page
  do.call(grid.arrange, c(plot_list[[metric]], ncol = 3))
}
#plot overall_Acc
p <- ggplot(plot_data, aes(overall_Acc, fill = source_model)) + 
  geom_histogram(position = "dodge", alpha = 0.5, colour = "black", bins = 15) +
  scale_fill_brewer(palette = "Set1",
                    labels = c("Counts", "Counts and AF", "Counts and imputed AF", "epiAneufinder")) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = paste("Overall Accuracy distribution for HMM models and epiAneufinder\n",
                     paste("(", cnv_type, ")", sep = "")), 
       x = metric, y = "Number of cells", fill = "Model")
print(p)

dev.off()
###############################################################################
# Plot the results on a karyogram
###############################################################################

#plot CNV predictions
#extract chromosome number and starting position
for (model in names(models)){
  HMM_preds <- models[[model]][["predictions"]]
  HMM_preds$chr <- gsub("chr([0-9]+):.*", "chr\\1", HMM_preds$region)
  HMM_preds$fragment <- as.numeric(gsub("chr[0-9]+:([0-9]+)", "\\1", HMM_preds$region))
  
  #add cell hierarchical clustering
  #transform the predictions table into a cell by region form
  wide_HMM_preds <- as.data.frame(pivot_wider(HMM_preds[, colnames(HMM_preds) %in% c("region", "cell", "HMM_pred")],
                                              names_from = region, values_from = HMM_pred))
  rownames(wide_HMM_preds) <- as.vector(wide_HMM_preds$cell)
  wide_HMM_preds <- wide_HMM_preds[,-1]
  #hierarchical clustering itself
  set.seed(seed_value)
  dist_matrix <- dist(wide_HMM_preds)
  dist_matrix[is.na(dist_matrix)] <- 0
  set.seed(seed_value)
  hc_counts <- hclust(dist_matrix, method = "ward.D")
  
  #extract the order of cells
  cell_order <- rownames(wide_HMM_preds)[hc_counts$order]
  cell_order <- data.frame(Y_pos = c(1:length(cell_order)), cell = cell_order)
  
  HMM_preds <- left_join(HMM_preds, cell_order, by = "cell")
  
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
  p <- ggplot(HMM_preds, aes(x = X_pos, y = Y_pos, fill = factor(HMM_pred))) +
    geom_tile() +
    geom_vline(xintercept = chromosome_label_positions - 0.5, color = "black", linetype = "dotted", size = 0.5) +
    labs(title = paste("CNV predictions by the", model, "HMM"),
         x = "Region",
         y = "Cell",
         fill = "CNV class") +
    scale_x_continuous(breaks = chromosome_label_positions, labels = unique(HMM_preds$chr),
                       expand = c(0.05, 0)) +
    scale_fill_manual(values = c("#2171B5", "#1fb424", "#b4421f"),
                      labels = c("loss", "base", "gain")) +
    theme_minimal() +
    theme(legend.position = "top", 
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank()) 
  
  #create a plot for the ground truth
  p_wgs <- ggplot(HMM_preds, aes(x = X_pos, y = 1, fill = factor(WGS_score))) +
    geom_tile() +
    geom_vline(xintercept = chromosome_label_positions - 0.5, color = "black", linetype = "dotted", size = 0.5) +
    labs(title = paste("CNV predictions by the", model, "HMM"),
         x = "Region",
         y = "WGS",
         fill = "CNV class") +
    scale_x_continuous(breaks = chromosome_label_positions, labels = unique(HMM_preds$chr),
                       expand = c(0.05, 0)) +
    scale_fill_manual(values = c("#2171B5", "#1fb424", "#b4421f"),
                      labels = c("loss", "base", "gain")) +
    theme_minimal() +
    theme(legend.position = "top", 
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank())
  
  #combine the plots
  combined_plot <- ggarrange(p, p_wgs + labs(title = NULL) + theme(legend.position="none", ), 
                             ncol = 1, heights = c(88, 12))
  
  ggsave(paste(output_dir, model, "_HMM_karyogram.pdf", sep = ""), 
         plot = combined_plot, width = 40, height = 20, units = "cm")
}


#create a caryogram for epiAneufinder results for comparison
#construct an input table
epiA_karyo_data <- cnv_labels_per_cell_filtered
epiA_karyo_data$WGS_score <- epiA_truth
epiA_karyo_data$region <- rownames(epiA_karyo_data)

#add cell hierarchical clustering
#transpose the data
t_epiA_karyo <- as.data.frame(t(epiA_karyo_data[, !(colnames(epiA_karyo_data) %in% c("WGS_score", "region"))]))

#hierarchical clustering itself
dist_matrix <- dist(t_epiA_karyo)
dist_matrix[is.na(dist_matrix)] <- 0
hc_counts <- hclust(dist_matrix, method = "ward.D")

#extract the order of cells
cell_order <- rownames(t_epiA_karyo)[hc_counts$order]
cell_order <- data.frame(Y_pos = c(1:length(cell_order)), cell = cell_order)

#create a long format dataframe
epiA_karyo_data <- pivot_longer(epiA_karyo_data, cols = -c(region, WGS_score), 
                                names_to = "cell", values_to = "pred")

#add Y_pos
epiA_karyo_data <- left_join(epiA_karyo_data, cell_order, by = "cell")

#add segment and chromosome columns
epiA_karyo_data$chr <- gsub("chr([0-9]+):.*", "chr\\1", epiA_karyo_data$region)
epiA_karyo_data$segment <- as.numeric(gsub("chr[0-9]+:([0-9]+)", "\\1", epiA_karyo_data$region))

#make chr column a factor
epiA_karyo_data$chr <- factor(epiA_karyo_data$chr, levels = paste0("chr", 1:22))

#reorder the dataframe
epiA_karyo_data <- epiA_karyo_data[order(epiA_karyo_data$chr, 
                                         epiA_karyo_data$segment), ]

#add row number and record positions for labels
epiA_karyo_data <- epiA_karyo_data %>%
  mutate(X_pos = as.numeric(factor(region, levels = unique(region))))

chromosome_label_positions <- epiA_karyo_data %>%
  filter(!duplicated(chr)) %>%
  pull(X_pos)

#create a plot for the predictions
p <- ggplot(epiA_karyo_data, aes(x = X_pos, y = Y_pos, fill = factor(pred))) +
  geom_tile() +
  geom_vline(xintercept = chromosome_label_positions - 0.5, color = "black", linetype = "dotted", size = 0.5) +
  labs(title = paste("CNV predictions by epiAneufinder"),
       x = "Region",
       y = "Cell",
       fill = "CNV class") +
  scale_x_continuous(breaks = chromosome_label_positions, labels = unique(epiA_karyo_data$chr),
                     expand = c(0.05, 0)) +
  scale_fill_manual(values = c("#2171B5", "#1fb424", "#b4421f"),
                    labels = c("loss", "base", "gain")) +
  theme_minimal() +
  theme(legend.position = "top", 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_blank()) 

#create a plot for the ground truth
p_wgs <- ggplot(epiA_karyo_data, aes(x = X_pos, y = 1, fill = factor(WGS_score))) +
  geom_tile() +
  geom_vline(xintercept = chromosome_label_positions - 0.5, color = "black", linetype = "dotted", size = 0.5) +
  labs(title = paste("CNV predictions by epiAneufinder"),
       x = "Region",
       y = "WGS",
       fill = "CNV class") +
  scale_x_continuous(breaks = chromosome_label_positions, labels = unique(epiA_karyo_data$chr),
                     expand = c(0.05, 0)) +
  scale_fill_manual(values = c("#2171B5", "#1fb424", "#b4421f"),
                    labels = c("loss", "base", "gain")) +
  theme_minimal() +
  theme(legend.position = "top", 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_blank())

#combine the plots
combined_plot <- ggarrange(p, p_wgs + labs(title = NULL) + theme(legend.position="none", ), 
                           ncol = 1, heights = c(88, 12))

ggsave(paste(output_dir, "epiAneufinder_karyogram.pdf", sep = ""), 
       plot = combined_plot, width = 40, height = 20, units = "cm")

################################################################################
# Make a pdf with per-cell plots
################################################################################

per_cell_summary(models, output_dir)
