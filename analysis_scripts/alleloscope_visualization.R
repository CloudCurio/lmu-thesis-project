#loading libraries
library(ggplot2)
library(dplyr)

#set the working directory
setwd("C:\\Users\\liber\\Desktop\\Study\\LMU\\Thesis Project - MariaCT's Lab\\Data\\Alleloscope_batch_run")

#load input files
theta_values <- read.csv("bin_by_cell_theta.csv")
colnames(theta_values)[which(colnames(theta_values) == "X")] <- "region"
seg_table <- readRDS("..//seg_table_500k_epiAneuFinder_SNU601.rds")
SNP_counts <- read.csv("SNP_counts.csv")

SNP_counts <- SNP_counts[,-1]
colnames(SNP_counts) <- c("region", "nSNP")

#obtain cnv values from the segmentation table
seg_table$region <- NA
for (i in 1:nrow(seg_table)){
  seg_table$region[i] <- paste("chr", seg_table$chr[i],":", seg_table$start[i], sep = "")
}

#add cnv class information to the theta_values
cnv_data <- data.frame(region = seg_table$region, cnv = seg_table$cnv)
theta_values <- merge(theta_values, cnv_data, all.x = T, by = "region", sort = F)
rownames(theta_values) <- theta_values$region
theta_values <- theta_values[, !names(theta_values) %in% "region"]

#add cnv class information to the SNP_counts
SNP_counts <- merge(SNP_counts, cnv_data, all.x = T, by = "region", sort = F)

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

#save a subset of cells with at least 50% of the regions
good_cells <- theta_values[, which(na_proportion <= 0.5)]
write.csv(good_cells, "high_content_cells.csv")

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

#plot SNP counts per region
#extract chromosome number and starting position
SNP_counts$chr <- gsub("chr([0-9]+):.*", "chr\\1", SNP_counts$region)
SNP_counts$start <- as.numeric(gsub("chr[0-9]+:([0-9]+)", "\\1", SNP_counts$region))

#make chr column a factor
SNP_counts$chr <- factor(SNP_counts$chr, levels = paste0("chr", 1:22))

#reorder the dataframe
SNP_counts <- SNP_counts[order(SNP_counts$chr, SNP_counts$start), ]

#add row number and record positions for labels
SNP_counts <- SNP_counts %>%
  mutate(X_pos = row_number())
chromosome_label_positions <- SNP_counts %>%
  filter(!duplicated(chr)) %>%
  pull(X_pos)

#create the plot
p <- ggplot(SNP_counts, aes(x = X_pos, y = nSNP, fill = cnv)) +
  geom_bar(stat = "identity") +
  geom_vline(xintercept = chromosome_label_positions - 0.5, color = "black", linetype = "dotted", size = 0.5) +
  labs(title = "nSNP per region",
       x = "Region",
       y = "nSNP") +
  scale_x_continuous(breaks = chromosome_label_positions, labels = SNP_counts$chr[chromosome_label_positions],
                     expand = c(0.05, 0)) +
  theme_minimal() +
  theme(legend.position = "top", axis.text.x = element_text(angle = 90, hjust = 1)) 

ggsave("SNP_by_region.pdf", plot = p, width = 40, height = 20, units = "cm")
