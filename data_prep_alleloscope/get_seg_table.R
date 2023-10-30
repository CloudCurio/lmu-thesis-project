#This script extracts the segmentation table from an .rds file, provided by Alleloscope
#To use, change the working directory and path to the object to the ones you need
setwd("C:\\Users\\liber\\Desktop\\Study\\LMU\\Thesis Project - MariaCT's Lab\\Data")

Obj_scDNA=readRDS("SNU601_dna.rds")
seg_table <- Obj_scDNA$seg_table_filtered
saveRDS(seg_table, file = "seg_table_filtered_SNU601.rds")
