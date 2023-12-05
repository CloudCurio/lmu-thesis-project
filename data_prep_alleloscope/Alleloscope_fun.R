################################################################################
#This function runs Alleloscope and is intended to be used for batch runs
#'@param wd path to where the working directory should be set up. Standard parameter
#'values are built in relation to this directory.
#'@param dir_path path to where the output directory should be created.
#'For batch runs it is recommended that dir_path is different between iterations.
#'@param chr_size_path path to the table with chr sizes.
#'@param barcodes_path path to the barcodes file.
#'@param vartrix_output_path path to the folder where three vartrix output .rds files
#'(alt_all, ref_all and var_all) are stored.
#'@param raw_counts_path path to the raw counts file.
#'@param seg_table segmentation table to be used. NOT a path, but a dataframe,
#'containing four columns: chr, start, end and length.
#'Since all other files should be identical in the bounds of one experiment, this
#'parameter is the one to do iterations on.
################################################################################

alleloscope_run <- function(wd = "//work//project//ladcol_014//thesis_cnvcalling//lmu-thesis-project",
                            dir_path = "//work//project//ladcol_014//thesis_cnvcalling//output//",
                            chr_size_path = "..//data//SNU601_scATACseq//hg38//chrom.sizes.txt",
                            barcodes_path = "..//data//SNU601_scATACseq//barcodes.tsv.gz",
                            vartrix_output_path = "..//data//SNU601_scATACseq//vartrix_out//",
                            raw_counts_path = '..//data//SNU601_scATACseq//bin_by_cell_mtx//100k_fragments_sub.txt',
                            seg_table = NULL){
  library(Matrix)
  library(Alleloscope) # load the library
  setwd(wd) # set path to the github folder
  
  dir.create(dir_path) # set up output directory
  
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
  ################################################################################
  #Load the input files
  ################################################################################
  print("start input file loading")
  size=read.table(chr_size_path, stringsAsFactors = F) # read size file
  size=size[1:22,]
  
  # SNP by cell matrices for ref and alt alleles
  barcodes=read.table(gzfile(barcodes_path), 
                      sep='\t', 
                      stringsAsFactors = F, 
                      header=F)
  alt_all=readRDS(paste(vartrix_output_path, "alt_all.rds", sep = ''))
  ref_all=readRDS(paste(vartrix_output_path, "ref_all.rds", sep = ''))
  var_all=readRDS(paste(vartrix_output_path, "var_all.rds", sep = ''))
  
  # bin by cell matrices for tumor and normal for segmentation
  raw_counts=read.table(raw_counts_path, sep=' ', header=T, row.names = 1,stringsAsFactors = F)
  colnames(raw_counts)=gsub("[.]","-", colnames(raw_counts))
  
  #Optional: read known cell identity from peaks + visualization:
  #cell_type=readRDS('data-raw/SU008/scATAC/cell_type_from_peaks.rds')
  #clust_order=plot_scATAC_cnv(raw_mat = raw_counts , cell_type = cell_type, normal_lab=c("endo","fibro"), size = size, plot_path = paste0(dir_path,"/cov_cna_plot.pdf"))
  
  ################################################################################
  #Create an Alleloscope object for analysis
  ################################################################################
  print("creating alleloscope object")
  Obj=Createobj(alt_all =alt_all, ref_all = ref_all, var_all = var_all,
                samplename='Sample', 
                genome_assembly="GRCh38", 
                dir_path=dir_path, 
                barcodes=barcodes, 
                size=size, 
                assay='scATACseq')
  #filter out cells and SNPs with too few read counts
  Obj_filtered=Matrix_filter(Obj=Obj, cell_filter=5, SNP_filter=5, min_vaf = 0.1, 
                             max_vaf = 0.9) 
  
  # suggest setting min_vaf=0.1 and max_vaf=0.9 when SNPs are called in the tumor 
  # sample for higher confident SNPs
  
  ################################################################################
  #Unbiased segmentation based on matched WES/WGS data
  ################################################################################
  print("starting segmentation")
  Obj_filtered$seg_table<-seg_table
  
  Obj_filtered=Segments_filter(Obj_filtered=Obj_filtered, nSNP=100, len = 100000)
  
  ################################################################################
  #Estimate cell major haplotype proportion for each region
  ################################################################################
  print("estimating regions")
  #estimates theta_hat for each cell of each region in seg_table_filtered
  Obj_filtered=Est_regions(Obj_filtered = Obj_filtered, max_nSNP = 30000, plot_stat = T,cont = FALSE)
  
  # Recommend max_nSNP <50000
  # Regions without allelic imbalence do not coverge (Reach the max number of iterations.)
  print("Theta-hat obtained")
  ################################################################################
  #Identify/Assign normal cells and diploid regions
  ################################################################################
  
  #Obj_filtered$ref=Obj_filtered$seg_table_filtered$chrr[7] # choose one normal region
  
  ##Optional: assign "normal cells" from scATAC-seq genome-wide peak signals
  ##Obj_filtered$select_normal$barcode_normal=cell_type[which(cell_type[,2]!='tumor'),1]
  
  ################################################################################
  #Genotype each cell in each region
  ################################################################################
  
  #Select normal cells
  #Obj_filtered=Select_normal(Obj_filtered = Obj_filtered, 
  #                           raw_counts=raw_counts, plot_theta = TRUE)
  #Estimate cell-specific (rho-hat, theta-hat) values for each region
  #Obj_filtered=Genotype_value(Obj_filtered = Obj_filtered, type='tumor', 
  #                            raw_counts=raw_counts, cov_adj=1,
  #                            ref_gtv = NULL,mincell = NULL,
  #                            qt_filter = TRUE,cell_filter = TRUE,
  #                            refr = TRUE,cov_only = FALSE)  # for tumor
  
  #Genotype all cells and generate a genotype plot for each region
  #Obj_filtered=Genotype(Obj_filtered = Obj_filtered, #cell_type=cell_type, 
  #                      xmax=3)
  
  ################################################################################
  #Construct lineage structure using cell major haplotype proportions 
  #for each cell across all regions
  ################################################################################
  
  #Generate lineage tree based on cell-specific genotypes across the regions
  # tmp=Select_normal(Obj_filtered = Obj_filtered, 
  #                   raw_counts=raw_counts, 
  #                   plot_theta = TRUE, 
  #                   #cell_type = cell_type, 
  #                   mincell = 0)
  # rm(tmp)
  # print("Lineages constructed! Done!")
  ################################################################################
  #Potential downstream analysis
  ################################################################################
  
  #Integrate allele-specific CNAs and chromatin accessibility at the single-cell level
  #UMAP projection using genome-wide peak profile on the tumor cells
  #umap_peak=readRDS("./data-raw/SU008/scATAC/peak_umap_tumor.rds")
  
  #Integrate allele-specific CNAs and peak signals for each cell in the scATAC-seq data
  #example: CNA on chr 4
  #theta_hat_chr4=Obj_filtered$rds_list$`chr4:0`$theta_hat
  #theta_hat_chr4=theta_hat_chr4[match(rownames(umap_peak), names(theta_hat_chr4))]
  #umap_peak$theta_hat=theta_hat_chr4
  
  #visualize two signals simultaneously for each cell in the scATAC-seq data
  # library(ggplot2)
  # library(RColorBrewer)
  # # UMAP
  # pp=ggplot(umap_peak,aes(x = UMAP1, y=UMAP2)) +
  #   geom_point(size=1,alpha=0.5, aes(color=(theta_hat))) +
  #   scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100))+
  #   theme_bw()
  # print(pp)
  # 
  # # density plot
  # pd <-ggplot(umap_peak, aes(x=theta_hat, color=peak_group)) +
  #   geom_density()+
  #   scale_color_manual(values = c("peak2" = "#F8766D","peak1" = "#00BFC4")) +
  #   theme_bw()
  # print(pd)
}

