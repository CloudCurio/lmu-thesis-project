################################################################################
#'This function grabs results .rds files created by Alleloscope and creates
#'a number of plots using the data from them. It can be used with batch runs
#'from Alleloscope_batch_run.R thanks to the chr and seg_table parameters 
#'(that script iterates through seg_table by chromosome and packages outputs
#'in folders for each of them)
#'@param wd path to where the output folder of Alleloscope_batch_run.R. chr
#'folders are there.
#'@param chrom chromosome number for batch runs. Allows to perform analysis by
#'chromosome.
#'@param seg_table segmentation table used by Alleloscope. When used in batch
#'analysis, only a subset for a particular chromosome is taken.
#'@param out_dir path to where the output folder should be created. By default,
#'results are saved in the current chr folder, in a subdirectory "summary_plots"
################################################################################

#wd value for my pc: "C:\\Users\\liber\\Desktop\\Study\\LMU\\Thesis Project - MariaCT's Lab\\Data"

plot_alleloscope <- function(wd = "/work/project/ladcol_014/thesis_cnvcalling/output/Alleloscope_batch/",
                             chrom = NULL,
                             seg_table = NULL,
                             out_dir = paste(chrom, "//summary_plots", sep = "")){
  #load libraries
  library(ggplot2)
  
  #set the working environment
  setwd(wd)
  
  #check the seg_table parameter for correctness
  if(is.null(seg_table)){
    stop("Segmentation table not selected")
  } else if (!is.data.frame(seg_table)){
    stop("seg_table must be a data frame")
  } else if (colnames(seg_table) != c("chr", "start", "end", "length")){
    if (all(c("chr", "start", "end", "length")) %in% colnames(seg_table)){
      seg_table <- seg_table[, c("chr", "start", "end", "length")]
      warning("More columns found than required, extra columns were removed")
    } else {
      stop('seg_table must have columns "chr", "start", "end" and "length')
    }
  }
  
  #read Alleloscope output files
  files_to_read <- list.files(path = paste("chr", chrom, "\\rds\\EMresults", sep = ""), 
                              pattern = "\\.rds$", full.names = T)
  all_frags <- list()
  all_frags <- lapply(files_to_read,function(x) {
    readRDS(file = x)
  })
  
  #read the segmentation table (for CNV labels)
  seg_table <- seg_table
  
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
  
  #visualize the distribution of mean theta hat values
  dir.create(out_dir)
  ggplot(frag_summary, aes(x = theta_hat_avg, y = length, fill = cnv, color = cnv)) +
    geom_col(alpha = 0.5) +
    theme_classic()
  ggsave(paste(out_dir, "//mean_theta.pdf", sep = ""))
  
  ggplot(frag_summary, aes(x = theta_hat_avg, fill = cnv, color = cnv)) +
    geom_histogram(bins = length(unique(frag_summary$theta_hat_avg)), alpha = 0.5) +
    theme_classic()
  ggsave(paste(out_dir, "//mean_theta_hist.pdf", sep = ""))
  
  #select a region for per cell analysis
  #region <- "chr14_93700001"
  for (region in names(all_frags)){
    reg_theta <- unlist(all_frags[[region]]['theta_hat'])
    foo <- data.frame("region" = region, "theta" = reg_theta)
    ggplot(foo, aes(x = reg_theta, fill = region))+
      geom_histogram(bins = 50, alpha = 0.5)+
      theme_classic()
    ggsave(paste(out_dir, "//", region, ".pdf", sep = ""))
  }
}