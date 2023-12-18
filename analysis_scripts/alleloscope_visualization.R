#loading libraries
library(ggplot2)

#set the working directory
setwd("C:\\Users\\liber\\Desktop\\Study\\LMU\\Thesis Project - MariaCT's Lab\\Data\\Alleloscope_batch_run")

#load input files
theta_values <- read.csv("bin_by_cell_theta.csv")
seg_table <- readRDS("..//seg_table_500k_epiAneuFinder_SNU601.rds")

#obtain cnv values from the segmentation table
seg_table$ID <- NA
for (i in 1:nrow(seg_table)){
  seg_table$ID[i] <- paste("chr", seg_table$chr[i],":", seg_table$start[i], sep = "")
}

#add cnv class information to the theta_values
cnv_data <- data.frame(X = seg_table$ID, cnv = seg_table$cnv)
theta_values <- merge(theta_values, cnv_data, all.x = T, by = "X", sort = F)
rownames(theta_values) <- theta_values$X
theta_values <- theta_values[, !names(theta_values) %in% "X"]

#check NA proportion for each cell
na_proportion <- c()
for (i in colnames(theta_values)){
  if (i == "cnv"){
    next
  }
  na_proportion <- append(na_proportion, sum(is.na(theta_values[,i]))/length(theta_values[,i]))
}
names(na_proportion) <- colnames(theta_values)[-ncol(theta_values)]
plot_data <- data.frame("val" = as.numeric(na_proportion))

ggplot(plot_data, aes(x = val)) +
  geom_histogram(binwidth = 0.02, fill = "#008080", color = "black") +
  labs(
    title = "NA proportion per cell",
    x = "NA proportion",
    y = "Frequency"
  ) +
  theme_minimal()

#plot mean values per region, colored by CNV type
plot_data <- data.frame("val" = rowMeans(theta_values[, !names(theta_values) %in% "cnv"], na.rm = T), "cnv" = theta_values$cnv)

ggplot(plot_data, aes(x = val, fill = cnv)) +
  geom_histogram(binwidth = 0.01, color = "black", position = "dodge") +
  labs(
    title = "Mean theta-hat values",
    x = "Theta-hat",
    y = "Frequency"
  ) +
  theme_minimal()

ggplot(data = plot_data, aes(x = val, fill = cnv)) +
  geom_histogram(binwidth = 0.01, position = "dodge", na.rm = TRUE) +
  labs(
    title = "Mean theta hat",
    x = "Values",
    y = "Frequency"
  ) +
  theme_minimal()
