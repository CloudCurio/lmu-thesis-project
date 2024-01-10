###############################################################################
# Load libraries
###############################################################################

library(depmixS4)
library(dplyr)
library(tidyr)

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
# Construct an HMM
###############################################################################

set.seed(2024)

#build a simple univariate model
train_mod <- depmix(response = counts ~ 1, data = dataset, nstates = 2,
                     family = poisson())
f_train <- fit(train_mod, verbose = TRUE, emc = em.control(rand = FALSE))
print(f_train)
summary(f_train)

#predict states for new data and estimate accuracy
cells <- unique(all_cells_data$cell)
cells <- cells[-which(cells == best_cell)]
set.seed(2024)
test_set <- all_cells_data[which(all_cells_data$cell == sample(cells, 1)),]
set.seed(2024)
test_mod <- depmix(response = counts ~ 1, data = test_set, nstates = 2,
                   family = poisson())
test_mod <- setpars(test_mod, getpars(f_train))

#try a multivariate model with gaussian distributions
multiV_mod <- depmix(response = list(counts ~ 1, theta_hat ~ 1), data = all_cells_data, 
                     nstates = 3, family = list(gaussian("identity"), multinomial()))
f_mV <- fit(multiV_mod, verbose = TRUE)
print(f_mV)
summary(f_mV)

#try a multivariate model with tailored distributions
mV_tuned_mod <- depmix(response = list(counts ~ 1, theta_hat ~ 1), data = dataset, 
                     nstates = 3, family = list(poisson(), gaussian()))
f_mVt <- fit(mV_tuned_mod, verbose = TRUE)
print(f_mVt)
summary(f_mVt)