setwd("//work//project//ladcol_014//thesis_cnvcalling//lmu-thesis-project")
source("data_prep_alleloscope//Alleloscope_fun.R")

seg_table_path <- "..//data//SNU601_scATACseq//seg_tables//seg_table_500k_epiAneuFinder_SNU601.rds"

input_table<-readRDS(seg_table_path)
input_table<-input_table[,c("chr","start","end","length")]

for (chr in unique(input_table$chr)){
  alleloscope_run(dir_path = paste("//work//project//ladcol_014//thesis_cnvcalling//output//",
                                   "alleloscope_batch//chr", chr),
                  seg_table = input_table[input_table$chr == chr,])
  file.rename(from = paste("//work//project//ladcol_014//thesis_cnvcalling//output//",
                    "alleloscope_batch//chr", chr, "//rds//theta_N_seg.rds"), 
              to = paste("//work//project//ladcol_014//thesis_cnvcalling//output//",
                         "alleloscope_batch//chr", chr, "//theta_N_seg.rds"))
  unlink(paste("//work//project//ladcol_014//thesis_cnvcalling//output//",
               "alleloscope_batch//chr", chr, "//rds//"), recursive = TRUE)
  unlink(paste("//work//project//ladcol_014//thesis_cnvcalling//output//",
               "alleloscope_batch//chr", chr, "//plots//"), recursive = TRUE)       
}