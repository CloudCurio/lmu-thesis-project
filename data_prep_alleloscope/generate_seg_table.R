################################################################################
#'This function gets a segmentation table from the epiAneufinder output folder
#'(results_table.tsv) and calculates a pseudobulk value, with CNV label cutoffs
#'at >=1.5 for gain and <=0.5 for loss. The values are then transformed into
#'character labels and the table is returned. Optionally, the table can be
#'saved at desired location
#'@param epiA_results_path - the path to the results_table.tsv
#'@param out_path - the path where the table should be saved for future use
#'(optional)
################################################################################

generate_seg_table <- function(epiA_results_path = "epiAneufinder runs//500kb bins//results_table.tsv",
                               out_path = NA){
  #load the input file
  results_table <- read.table(epiA_results_path, header = TRUE)
  
  #calculate mean cnv values
  results_table$cnv <- rowMeans(results_table[, !(colnames(results_table) %in% c("seq", "start", "end"))])
  
  #maintenance changes, such as renaming columns and adding length
  colnames(results_table)[1] <- "chr"
  results_table$length <- results_table$end - results_table$start + 1
  
  #change cnv values to character labels
  results_table$cnv <- ifelse(results_table$cnv >= 1.5, "gain", 
                ifelse(results_table$cnv <= 0.5, "loss", "base"))
  
  #save the resulting file
  if (!is.na(out_path)){
    saveRDS(results_table, out_path)
  }
  
  return(results_table)
} 