wd_path <- "//work//project//ladcol_014//thesis_cnvcalling//data//SNU601_scATACseq"

loss_thresh <- 0.5
gain_thresh <- 1.5
  
#set the path and load the input file
setwd(wd_path)

epiAneufinder_res <- read.table("SNU601_br7_epiAneufinder_results_table.tsv", 
                                header = T)

#get mean scores for each region and label the 
epiAneufinder_res$states <- rowMeans(epiAneufinder_res[,-(1:3)])
epiAneufinder_res$cnv <- ifelse(epiAneufinder_res$states >= gain_thresh, "gain",
                                ifelse(epiAneufinder_res$states <= loss_thresh, 
                                       "loss", "base"))

#join neighbouring bins on each of the chromosomes
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

write.table(seg_table, ".//seg_table_epiAneuFinder_SNU601.tsv")