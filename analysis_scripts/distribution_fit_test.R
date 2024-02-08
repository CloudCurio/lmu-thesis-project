#this script is for testing the response fit to different distributions 

###############################################################################
# Load libraries
###############################################################################

library(fitdistrplus)
library(MASS)

###############################################################################
# Define settings values
###############################################################################

input_dir <- "C:\\Users\\liber\\Desktop\\Study\\LMU\\Thesis Project - MariaCT's Lab\\Data\\HMM inputs\\Alleloscope_1mb_batch\\"
output_dir <- paste("..//..//HMM outputs//Alleloscope_1mb_batch//", "fitness_tests//", sep = "")

#change to 2 for the two-state model
nstates <- 2
#options: counts, theta_hat
responses <- c("counts", "theta_hat")

#do we want imputation?
impute_theta_flag <- TRUE
test_imputed <- FALSE

setwd(input_dir)
dir.create(output_dir)

###############################################################################
# Load the input data
###############################################################################

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

if (test_imputed){
  dataset <- all_cells_data_full
} else {
  dataset <- all_cells_data_filtered
}

#subset data by state
state_subdata_list <- list()
if (nstates == 2){
  state_subdata_list[["base"]] <- dataset[which(dataset$wgs_score<2.5),]
  state_subdata_list[["gain"]] <- dataset[which(dataset$wgs_score>=2.5),]
} else if (nstates == 3){
  state_subdata_list[["loss"]] <- dataset[which(dataset$wgs_score<=1.5),]
  state_subdata_list[["base"]] <- dataset[which(dataset$wgs_score<2.5) && which(dataset$wgs_score>1.5),]
  state_subdata_list[["gain"]] <- dataset[which(dataset$wgs_score>=2.5),]
}

for (cnv_state in names(state_subdata_list)){
  data_sub <- state_subdata_list[[cnv_state]]
  
  for (resp in responses[1]){
    print(paste(cnv_state, resp))
    pdf(paste(output_dir, nstates, resp, cnv_state, "fits.pdf", sep = "_"))
    
    data_vector <- unlist(data_sub[,resp])
    if (resp == "theta_hat"){
      data_vector <- data_vector[!is.na(data_vector)]
      warning("Please consider that for AF fit NA values are dropped")
    }
    
    #poisson
    fit_pois <- fitdist(data_vector, "pois")
    plotdist(data_vector, "pois", para = list(lambda = fit_pois$estimate))
    gofstat(fit_pois)
    
    #neg.binom.
    fit_negbin <- fitdist(data_vector, "nbinom", 
                          start = list(size = 1, mu = mean(data_vector)))
    plotdist(data_vector, "nbinom", 
             para = list(size = fit_negbin$estimate[1], mu = fit_negbin$estimate[2]))
    gofstat(fit_negbin)
    
    #log gaussian
    log_data <- log(data_vector)
    fit_normal <- fitdist(log_data, "norm")
    plotdist(log_data, "norm", para = list(mean = fit_normal$estimate["mean"],
                                           sd = fit_normal$estimate["sd"]))
    foo <- gofstat(fit_normal)
    ks_result <- ks.test(log_data, "pnorm", mean = fit_normal$estimate["mean"], 
                         sd = fit_normal$estimate["sd"])
    
    dev.off()
  }
  #temporary theta hat solution
  if ("theta_hat" %in% responses){
    warning("Please consider that for AF fit NA values are dropped")
    pdf(paste(output_dir, nstates, "AF", cnv_state, "fits.pdf", sep = "_"))
    data_vector <- data_sub$theta_hat[!is.na(data_sub$theta_hat)]
    fit_normal <- fitdist(data_vector, "norm")
    plotdist(data_vector, "norm", para = list(mean = fit_normal$estimate["mean"],
                                               sd = fit_normal$estimate["sd"]))
    foo <- gofstat(fit_normal)
    ks_result <- ks.test(data_vector, "pnorm", mean = fit_normal$estimate["mean"], 
                         sd = fit_normal$estimate["sd"])
    
    dev.off()
  }
}

#figure out why theta hat estimation fails

