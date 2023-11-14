#This script reads the result file from epiAneufinder and converts
#it into a format suitable for Alleloscope. It is possible to also
#merge neighbouring bins of the same identity to obtain regions.

#wd_path <- "//work//project//ladcol_014//thesis_cnvcalling//data//SNU601_scATACseq"
wd_path <- "C:\\Users\\liber\\Desktop\\Study\\LMU\\Thesis Project - MariaCT's Lab\\Data"

to_merge <- F

loss_thresh <- 0.5
gain_thresh <- 1.5
  
#set the path and load the input file
setwd(wd_path)

epiAneufinder_res <- read.table("epiAneufinder runs\\500kb bins\\results_table.tsv", 
                                header = T)

#get mean scores for each region and label the 
epiAneufinder_res$states <- rowMeans(epiAneufinder_res[,-(1:3)])
epiAneufinder_res$cnv <- ifelse(epiAneufinder_res$states >= gain_thresh, "gain",
                                ifelse(epiAneufinder_res$states <= loss_thresh, 
                                       "loss", "base"))

#join neighbouring bins on each of the chromosomes or just transfer the data
if (to_merge == T){
  seg_table <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(seg_table) <- c("chr", "start", "end", "length", "cnv")
  
  for (bin in 1:nrow(epiAneufinder_res)){
    if (nrow(seg_table) == 0){                              #create the first line
      seg_table[1,] <- NA
      
      seg_table$chr[1] <- epiAneufinder_res$seq[bin]
      seg_table$start[1] <- epiAneufinder_res$start[bin]
      seg_table$cnv[1] <- epiAneufinder_res$cnv[bin]
    } else if (epiAneufinder_res$cnv[bin] != seg_table$cnv[nrow(seg_table)] ||
               epiAneufinder_res$seq[bin] != seg_table$chr[nrow(seg_table)]){ #end line if new bin has different cnv value OR chr
      seg_table$end[nrow(seg_table)] <- epiAneufinder_res$end[(bin - 1)]
      seg_table$length[nrow(seg_table)] <- seg_table$end[nrow(seg_table)]-seg_table$start[nrow(seg_table)]+1
      
      seg_table[nrow(seg_table)+1,] <- NA
      seg_table$chr[nrow(seg_table)] <- epiAneufinder_res$seq[bin]
      seg_table$start[nrow(seg_table)] <- epiAneufinder_res$start[bin]
      seg_table$cnv[nrow(seg_table)] <- epiAneufinder_res$cnv[bin]
    } 
    if (bin == nrow(epiAneufinder_res)){                    #end the last segment when the last bin is reached
      seg_table$end[nrow(seg_table)] <- epiAneufinder_res$end[bin]
      seg_table$length[nrow(seg_table)] <- seg_table$end[nrow(seg_table)]-seg_table$start[nrow(seg_table)]+1
    }
  }
} else {
  seg_table <- data.frame("chr" = epiAneufinder_res$seq,
                          "start" = epiAneufinder_res$start,
                          "end" = epiAneufinder_res$end,
                          "length" = (epiAneufinder_res$end - epiAneufinder_res$start +1),
                          "cnv" = epiAneufinder_res$cnv)
}

#change chromosome notations to numbers
seg_table$chr <- gsub("chr", "", seg_table$chr)

saveRDS(seg_table, "seg_table_500k_epiAneuFinder_SNU601.rds")