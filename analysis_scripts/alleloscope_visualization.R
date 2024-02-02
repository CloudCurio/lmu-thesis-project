#loading libraries
library(ggplot2)
library(dplyr)
library(Seurat)

#set the working directory
setwd("C:\\Users\\liber\\Desktop\\Study\\LMU\\Thesis Project - MariaCT's Lab\\Data\\Alleloscope_runs\\Alleloscope_batch_fixed")

#load input files
theta_values <- read.csv("bin_by_cell_theta.csv")
colnames(theta_values)[which(colnames(theta_values) == "X")] <- "region"

seg_table <- readRDS("..//..//seg_tables//seg_table_500k_epiAneuFinder_SNU601.rds")
SNP_counts <- read.csv("SNP_counts.csv")
count_matrix <- readRDS("..//..//epiAneufinder runs//500kb bins//counts_gc_corrected.rds")
colnames(count_matrix) <- sub("^cell-", "", colnames(count_matrix))

count_summary <- readRDS("..//..//epiAneufinder runs//500kb bins//count_summary.rds")
rowinfo <- rowRanges(count_summary)

wgs_cnv_class <- read.csv("..//..//wgs_results_formated.csv")
wgs_cnv_class$X <- NULL
wgs_cnv_class$region <- apply(wgs_cnv_class, 1, function(row) {
  paste("chr", row['chr'], ":", row['start'], sep = "")
})
wgs_cnv_class$wgs_score <-  1*wgs_cnv_class$loss_wgs+2*wgs_cnv_class$base_wgs+3*wgs_cnv_class$gain_wgs

SNP_counts <- SNP_counts[,-1]
colnames(SNP_counts) <- c("region", "nSNP")

#reform the per-cell CNV calls
per_cell_CNVs <- as.data.frame(readRDS("..//..//epiAneufinder runs//500kb bins//cnv_calls.rds"))
colnames(per_cell_CNVs) <- sub("^cell.", "", colnames(per_cell_CNVs))
colnames(per_cell_CNVs) <- gsub("\\.", "-", colnames(per_cell_CNVs))
per_cell_CNVs <- per_cell_CNVs+1

#obtain cnv values from the segmentation table
seg_table$region <- NA
for (i in 1:nrow(seg_table)){
  seg_table$region[i] <- paste("chr", seg_table$chr[i],":", seg_table$start[i], sep = "")
}
rownames(per_cell_CNVs) <- seg_table$region
write.csv(per_cell_CNVs, "..//..//HMM inputs//epiAneufinder_per_cell_cnv.csv")

#add Alleloscope cnv class information to the theta_values
cnv_data <- data.frame(region = seg_table$region, cnv = seg_table$cnv)
theta_values <- merge(theta_values, cnv_data, all.x = T, by = "region", sort = F)
rownames(theta_values) <- theta_values$region
theta_values <- theta_values[, !names(theta_values) %in% "region"]
write.csv(theta_values, "HMM_theta_full.csv")

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
print(summary(na_proportion))
plot_data <- data.frame("val" = as.numeric(na_proportion))

p <- ggplot(plot_data, aes(x = val)) +
  geom_histogram(binwidth = 0.02, fill = "#008080", color = "black") +
  labs(
    title = "NA proportion per cell",
    x = "NA proportion",
    y = "Frequency"
  ) +
  theme_minimal()
ggsave("NA_proportion.pdf", p)

#save a subset of cells with at least 50% of the regions
good_cells <- theta_values[, which(na_proportion <= 0.5)]
good_cells$cnv <- theta_values$cnv
colnames(good_cells) <- sub("\\.", "-", colnames(good_cells))
write.csv(good_cells, "HMM_theta_high_content.csv")

#subset the count matrix by selected cells
#introduce rownames from the summary GRanges
for (fragment in 1:nrow(seg_table)){
  rownames(count_matrix)[fragment] <- paste(rowinfo$wSeq[fragment], rowinfo$wStart[fragment], sep = ":")
}

#keep only the regions present in Alleloscope results
good_counts <- as.data.frame(count_matrix[rownames(count_matrix) %in% rownames(good_cells),])
rownames(good_counts) <- rownames(count_matrix)[rownames(count_matrix) %in% rownames(good_cells)]
write.csv(good_counts, "HMM_read_counts_full.csv")

#leave only the previously selected cells
good_counts <- good_counts[, colnames(good_counts) %in% colnames(good_cells)]
foo <- good_counts[, which(colnames(good_counts) %in% colnames(good_cells)), drop = FALSE]
write.csv(good_counts, file = "HMM_read_counts_filtered.csv")

#plot mean values per region, colored by CNV type
plot_data <- data.frame("val" = rowMeans(theta_values[, !names(theta_values) %in% "cnv"], na.rm = T), "cnv" = theta_values$cnv)

p <- ggplot(plot_data, aes(x = val, fill = cnv)) +
  geom_histogram(binwidth = 0.01, color = "black", position = "dodge") +
  labs(
    title = "Mean theta-hat values",
    x = "Theta-hat",
    y = "Frequency"
  ) +
  theme_minimal()
ggsave("mean_theta_per_region.pdf", p)

#TODO:plot reads per region

#create a file with CNV label comparison (WGS and epiAneufinder predictions)
cnv_data <- merge(cnv_data, wgs_cnv_class, by.x = "region", sort = F)
cnv_data <- cnv_data[which(cnv_data$region %in% rownames(good_cells)),colnames(cnv_data) %in% c("region", "cnv", "wgs_score")]
cnv_data$cnv <- factor(cnv_data$cnv, levels = c("loss", "base", "gain"))
cnv_data$cnv <- as.integer(cnv_data$cnv)
write.csv(cnv_data, file = "cnv_labels_comparison.csv")

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
