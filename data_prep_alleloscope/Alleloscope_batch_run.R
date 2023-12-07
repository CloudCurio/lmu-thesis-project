setwd("//work//project//ladcol_014//thesis_cnvcalling//lmu-thesis-project")
#setwd("C:\\Users\\liber\\Desktop\\Study\\LMU\\Thesis Project - MariaCT's Lab\\Data")

source("data_prep_alleloscope//Alleloscope_fun.R")
source("analysis_scripts//alleloscope_summary.R")

seg_table_path <- "..//data//SNU601_scATACseq//seg_tables//seg_table_500k_epiAneuFinder_SNU601.rds"
#seg_table_path <- "seg_table_500k_epiAneuFinder_SNU601.rds"

input_table<-readRDS(seg_table_path)
input_table<-input_table[,c("chr","start","end","length")]

#for (chr in unique(input_table$chr)){
for (chr in c(2)){
  alleloscope_run(dir_path = paste("//work//project//ladcol_014//thesis_cnvcalling//output//",
                                   "Alleloscope_batch//chr", chr, sep = ''),
                  seg_table = input_table[input_table$chr == chr,])
  summarize_alleloscope(chrom = paste("chr", chr, sep = ""),
                   seg_table = input_table[input_table$chr == chr,])
  print("summary files saved successfully")
  #file.rename(from = paste("//work//project//ladcol_014//thesis_cnvcalling//output//Alleloscope_batch//chr", 
  #                          chr, "//rds//theta_N_seg.rds", sep = ''), 
  #             to = paste("//work//project//ladcol_014//thesis_cnvcalling//output//Alleloscope_batch//chr", 
  #                        chr, "//theta_N_seg.rds", sep = ''))
  unlink(paste("//work//project//ladcol_014//thesis_cnvcalling//output//Alleloscope_batch//chr",
              chr, "//rds//", sep = ''), recursive = TRUE)
  unlink(paste("//work//project//ladcol_014//thesis_cnvcalling//output//Alleloscope_batch//chr",
              chr, "//plots//", sep = ''), recursive = TRUE)
}