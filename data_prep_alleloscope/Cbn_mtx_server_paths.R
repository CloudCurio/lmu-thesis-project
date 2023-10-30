setwd("//work//project//ladcol_014//thesis_cnvcalling//data//SNU601_scATACseq")
source("..//..//lmu-thesis-project//data_prep_alleloscope//Cbn_matrix.R")

Cbn_matrix(dir_path = "./", 
           vcf_path = "./vcf/SNU601/",
           barcodes_path = gzfile("./barcodes.tsv.gz"),
           samplename = "SNU601_scATAC",
           plot_stat = TRUE)