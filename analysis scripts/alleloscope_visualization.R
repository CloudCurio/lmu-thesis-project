library(ggplot2)

setwd("C:\\Users\\liber\\Desktop\\Study\\LMU\\Thesis Project - MariaCT's Lab\\Data")

#read Alleloscope output files
files_to_read <- list.files(path = "Alleloscope - EpiAneufinder bins merged\\rds\\EMresults",pattern = "\\.rds$",full.names = T)
all_frags <- list()
all_frags <- lapply(files_to_read,function(x) {
  readRDS(file = x)
})

#read the segmentation table (for CNV labels)
seg_table <- readRDS("seg_table_epiAneuFinder_SNU601.rds")

#add seq names to the Alleloscope output objects
names <- sapply(files_to_read, basename)
names <- gsub(".rds", "", names)
names(all_frags) <- names

#create a dataframe for the reformated outputs
frag_summary <- as.data.frame(matrix(nrow = length(all_frags), ncol = 3))
colnames(frag_summary) <- c("chr", "start", "theta_hat_avg")

#split names and starting positions
split_names <- strsplit(names(all_frags), "_")
chr <- c()
start <- c()

for (i in 1:length(split_names)){
  chr <- append(chr, split_names[[i]][1])
  start <- append(start, split_names[[i]][2])
}

#remove "chr" from chromosome identifiers
chr <- gsub("chr", "", chr)

frag_summary$chr <- as.numeric(chr)
frag_summary$start <- as.numeric(start)

#get average theta hat values for each fragment
for (frag in 1:length(all_frags)){
  frag_summary$theta_hat_avg[frag] <- mean(mean(unlist(all_frags[[frag]]["theta_hat"])))  
}

#merge Alleloscope output with the segmentation table to obtain EpiAneufinder-predicted cnv states (as well as more complete size metrics)
frag_summary <- merge(frag_summary, seg_table)
frag_summary <- frag_summary[order(frag_summary$chr, frag_summary$start),]

#visualize the distribution of theta hat values
ggplot(frag_summary, aes(x = theta_hat_avg, y = length, fill = cnv, color = cnv)) +
  geom_col(alpha = 0.5) +
  theme_classic()
ggsave("plots//alleloscope_merged//100kb_all_regions.pdf")

ggplot(frag_summary, aes(x = theta_hat_avg, fill = cnv, color = cnv)) +
  geom_histogram(bins = length(unique(frag_summary$theta_hat_avg)), alpha = 0.5) +
  theme_classic()
ggsave("plots//alleloscope_merged//100kb_all_regions_hist.pdf")

#select a region for per cell analysis
#region <- "chr14_93700001"
for (region in names(all_frags)){
  reg_theta <- unlist(all_frags[[region]]['theta_hat'])
  foo <- data.frame("region" = region, "theta" = reg_theta)
  ggplot(foo, aes(x = reg_theta, fill = region))+
    geom_histogram(bins = 50, alpha = 0.5)+
    theme_classic()
  ggsave(paste("plots//alleloscope_merged//100kb_", region, ".pdf"))
}
