setwd("C:\\Users\\liber\\Desktop\\Study\\LMU\\Thesis Project - MariaCT's Lab\\Data")

Obj_scDNA=readRDS("SNU601_dna.rds")
seg_table <- Obj_scDNA$seg_table_filtered
write.table(seg_table, file = "seg_table_filtered_SNU601")
