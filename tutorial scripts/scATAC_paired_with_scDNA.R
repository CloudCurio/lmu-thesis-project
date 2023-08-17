#Step 0: load input files
library(Alleloscope) # load the library
setwd(".//Alleloscope//") # set path to the github folder

dir_path <- ".//samples//SNU601//scATAC//output//"; dir.create(dir_path) # set up output directory

data(centromere.GRCh38)
data(telomere.GRCh38)
size=read.table("data-raw//sizes.cellranger-GRCh38-1.0.0.txt", stringsAsFactors = F)

#read example files
# SNP by cell matrices for ref and alt alleles
barcodes=read.table("data-raw//SNU601//scATAC//barcodes.tsv", sep='\t', stringsAsFactors = F, header=F)
alt_all=readMM("data-raw//SNU601//scATAC//alt_all.mtx")
ref_all=readMM("data-raw//SNU601//scATAC//ref_all.mtx")
var_all=read.table("data-raw//SNU601//scATAC//var_all.vcf", header = F, sep='\t', stringsAsFactors = F)

# bin by cell matrices for tumor and normal for segmentation
raw_counts=read.table('data-raw//SNU601//scATAC//chr200k_fragments_sub.txt', sep='\t', header=T, row.names = 1,stringsAsFactors = F)
# Without paired normal sample, use matched scDNA-seq result to help normalize coverge for scATAC-seq data.

#load matched scDNA-seq
Obj_scDNA=readRDS("data-raw//SNU601//scATAC//SNU601_dna.rds")

#Step 1: Creating a Alleloscope object for the analysis
Obj=Createobj(alt_all =alt_all, ref_all = ref_all, var_all = var_all ,samplename='Sample', genome_assembly="GRCh38", dir_path=dir_path, barcodes=barcodes, size=size, assay='scATACseq')

#Filter out cells and SNPs with too few read counts
Obj_filtered=Matrix_filter(Obj=Obj, cell_filter=5, SNP_filter=5, centro=centromere.GRCh38, telo=telomere.GRCh38) 

# Since phasing information is estimated in the matched scDNA-seq dataset, 
# loose filter: cell_filter=5 and SNP_filter=5 can be used.  
# No further filter for extreme VAF values is needed.

#Step2: Segmentation results from matched scDNA-seq or WGS/WES
Obj_filtered$seg_table_filtered=Obj_scDNA$seg_table_filtered

#Step3: Estimate cell major haplotype proportion for each region
Obj_filtered = Est_regions(Obj_filtered = Obj_filtered, max_nSNP = 30000, min_cell = 20, phases = Obj_scDNA$rds_list, plot_stat = T, cont = TRUE)

# The phases for each SNP estimated from DNA sequencing data can help estimate the major haplotype proportion for each cell in scATAC-seq data. 
# Recommend max_nSNP <50000
# Regions without allelic imbalence do not coverge (Reach the max number of iterations.)

#Step4: Retrieve a diploid region from DNA-seq data
Obj_filtered$ref=Obj_scDNA$ref # choose one normal region

#Step5: Genotype each cell in each region
#Estimate cell-specific (rho_hat, theta_hat) values for each region.
#Set ref_gv = "genotype_values" (from scDNA-seq) to help with rho_hat estimation for scATAC-seq data.

Obj_filtered=Genotype_value(Obj_filtered = Obj_filtered, type='cellline', raw_counts=raw_counts, cov_adj =1 ,ref_gtv = Obj_scDNA$genotype_values) 

#Genotype all cells for each region and generate a genotype plot
#Set ref_gt = "genotypes" (from scDNA-seq) to help estimate haplotype profiles for each cell in scATAC-seq data.

Obj_filtered=Genotype(Obj_filtered = Obj_filtered, ref_gt = Obj_scDNA$genotypes,xmax=4)

#Step6: Infer clonal identity for each cell in the scATAC-seq data
clone.genotypes=readRDS(".//data-raw//SNU601//scATAC//clone.genotypes.rds")

Obj_filtered=AssignClones_ref(Obj_filtered=Obj_filtered, clone.genotypes=clone.genotypes)

#Potential downstream analysis
#Integrate DNA-level subclones and chromatin accessibility at the single-cell level
umap_peak=readRDS("./data-raw/SNU601/scATAC/peak_umap.rds")

Clone=Obj_filtered$cloneAssign$cloneAssign[match(rownames(umap_peak), names(Obj_filtered$cloneAssign$cloneAssign))]
umap_peak=cbind(umap_peak, Clone)

#identify differential accessible peaks (DAPs) between two clones with/without adjusting
#copy numbers: Basic CNV SNV:: Test_clonalDAP