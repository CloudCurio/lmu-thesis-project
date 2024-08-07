#Step 0: Load the input files
library(Alleloscope) # load the library
setwd(".//Alleloscope//") # set path to the github folder

dir_path <- ".//samples//SU008//scATAC//output//"; dir.create(dir_path) # set up output directory

size=read.table("data-raw//sizes.cellranger-atac-hg19-1.2.0.txt", stringsAsFactors = F) # read size file
size=size[1:22,]

# SNP by cell matrices for ref and alt alleles
barcodes=read.table("data-raw//SU008//scATAC//barcodes.tsv", sep='\t', stringsAsFactors = F, header=F)
alt_all=readMM("data-raw//SU008//scATAC//alt_all.mtx")
ref_all=readMM("data-raw//SU008//scATAC//ref_all.mtx")
var_all=read.table("data-raw//SU008//scATAC//var_all.vcf", header = F, sep='\t', stringsAsFactors = F)

# bin by cell matrices for tumor and normal for segmentation
raw_counts=read.table('data-raw//SU008//scATAC//chr200k_fragments_sub.txt', sep='\t', header=T, row.names = 1,stringsAsFactors = F)
colnames(raw_counts)=gsub("[.]","-", colnames(raw_counts))

cell_type=readRDS('data-raw//SU008//scATAC//cell_type_from_peaks.rds')

clust_order=plot_scATAC_cnv(raw_mat = raw_counts , cell_type = cell_type, normal_lab=c("endo","fibro"), size = size, plot_path = paste0(dir_path,"/cov_cna_plot.pdf"))

#Step 1:
Obj=Createobj(alt_all =alt_all, ref_all = ref_all, var_all = var_all ,samplename='Sample', genome_assembly="GRCh37", dir_path=dir_path, barcodes=barcodes, size=size, assay='scATACseq')

Obj_filtered=Matrix_filter(Obj=Obj, cell_filter=5, SNP_filter=5, min_vaf = 0.1, max_vaf = 0.9) 

# suggest setting min_vaf=0.1 and max_vaf=0.9 when SNPs are called in the tumor sample for higher confident SNPs

#Step 2: Unbiased segmentation based on matched WES/WGS data
Obj_filtered$seg_table=readRDS(".//data-raw//SU008//scATAC//seg_table_WES.rds")

Obj_filtered=Segments_filter(Obj_filtered=Obj_filtered, nSNP=500)

#Step 3: Estimate cell major haplotype proportion for each region
Obj_filtered=Est_regions(Obj_filtered = Obj_filtered, max_nSNP = 30000, plot_stat = T,cont = FALSE)

# Recommend max_nSNP <50000
# Regions without allelic imbalence do not coverge (Reach the max number of iterations.)

#Step 4: Identify/Assign normal cells and diploid regions
Obj_filtered$ref=Obj_filtered$seg_table_filtered$chrr[7] # choose one normal region

Obj_filtered$select_normal$barcode_normal=cell_type[which(cell_type[,2]!='tumor'),1]

#Step 5: Genotype each cell in each regions
Obj_filtered=Genotype_value(Obj_filtered = Obj_filtered, type='tumor', raw_counts=raw_counts, cov_adj=1)  # for tumor

Obj_filtered=Genotype(Obj_filtered = Obj_filtered, cell_type=cell_type, xmax=3)

#Step 6: Construct lineage structure using cell major haplotype 
#proportions for each cell across regions
tmp=Select_normal(Obj_filtered = Obj_filtered, raw_counts=raw_counts, plot_theta = TRUE, cell_type = cell_type, mincell = 0)
rm(tmp)

#Potential downstream analysis
#Integrate allele-specific CNAs and chromatin accessibility at the single-cell
#level
umap_peak=readRDS("./data-raw/SU008/scATAC/peak_umap_tumor.rds")

#chr 4 as example
theta_hat_chr4=Obj_filtered$rds_list$`chr4:0`$theta_hat
theta_hat_chr4=theta_hat_chr4[match(rownames(umap_peak), names(theta_hat_chr4))]
umap_peak$theta_hat=theta_hat_chr4

library(ggplot2)
library(RColorBrewer)
# UMAP
pp=ggplot(umap_peak,aes(x = UMAP1, y=UMAP2)) +
  geom_point(size=1,alpha=0.5, aes(color=(theta_hat))) +
  scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100))+
  theme_bw()
print(pp)

# density plot
pd <-ggplot(umap_peak, aes(x=theta_hat, color=peak_group)) +
  geom_density()+
  scale_color_manual(values = c("peak2" = "#F8766D","peak1" = "#00BFC4")) +
  theme_bw()
print(pd)