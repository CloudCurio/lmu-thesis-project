#set the working directory
setwd("C:\\Users\\liber\\Desktop\\Study\\LMU\\Thesis Project - MariaCT's Lab\\Data\\Alleloscope_batch_fixed")

#create a list to contain sub-dataframes
output_list <- list()

#create a dataframe for SNP counts per region
SNP_counts <- data.frame("region" = character(0), "nSNP" = numeric(0))

#create a dataframe of region by cell theta hat values
for (chr in c(1:22)){
  #read the summary RDS file
  #file existence check
  file_check <- list.files(paste("chr", chr, "//summary_files//", sep = ""), 
                           pattern = "\\.rds$", full.names = T)
  if (length(file_check) == 0){
    print(paste("no segments for chr", chr))
    next
  }
  
  summary <- readRDS(paste("chr", chr, "//summary_files//chr", chr, "all_frags.rds", sep = ""))
  
  #construct a vector of all present barcodes
  barcodes <- c()
  for (segment in summary){
    barcodes <- append(barcodes, segment[['barcodes']])
  }
  barcodes <- unique(barcodes)
  barcodes <- barcodes[order(barcodes)]
  
  #create a dataframe for theta_hat values
  theta_summary <- data.frame(matrix(nrow = 0, ncol = length(barcodes))) 
  colnames(theta_summary) = barcodes
  
  #add data per segment
  for (segment in names(summary)){
    theta_slice <- summary[[segment]]$theta_hat
    theta_slice <- as.data.frame(t(theta_slice))
    theta_summary <- merge(theta_summary, theta_slice, all.x = T, all.y = T, sort = F)
    
    #get SNP numbers per region
    SNP_counts <- rbind(SNP_counts, c(segment, length(summary[[segment]]$SNPs)))
  }
  rownames(theta_summary) <- names(summary)
  
  #add data for chr to the output list
  output_list[[chr]] <- theta_summary
}

#collect all barcodes and segment IDs
total_barcodes <- c()
total_segments <- c()

for (chr_df in output_list){
  total_barcodes <- append(total_barcodes, colnames(chr_df))
  total_segments <- append(total_segments, rownames(chr_df))
}
total_barcodes <- unique(total_barcodes)
total_barcodes <- total_barcodes[order(total_barcodes)]

#construct a final df
output_df <- data.frame(matrix(nrow = 0, ncol = length(total_barcodes)))
colnames(output_df) <- total_barcodes
for (chr_df in output_list){
  print(paste("out:", nrow(output_df)))
  print(paste("chr:", nrow(chr_df)))
  print(paste("sum", (nrow(output_df)+nrow(chr_df))))
  if (!is.null(chr_df)){
    output_df <- merge(output_df, chr_df, all.x = T, all.y = T, sort = F)
  }
}
rownames(output_df) <- total_segments

write.csv(output_df, file = "bin_by_cell_theta.csv")
write.csv(SNP_counts, file = "SNP_counts.csv")
