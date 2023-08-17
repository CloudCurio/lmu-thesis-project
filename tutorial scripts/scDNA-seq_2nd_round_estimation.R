library(Alleloscope) # load the library
setwd("./Alleloscope/") # set path to the github folder

dir_path <- ".//samples//P5931//scDNA//output//"; dir.create(dir_path) # set up output directory

data(centromere.GRCh38)
data(telomere.GRCh38)
size=read.table("data-raw//sizes.cellranger-GRCh38-1.0.0.txt", stringsAsFactors = F)

# SNP by cell matrices for ref and alt alleles
barcodes=read.table("data-raw//P5931//scDNA//barcodes_sub.tsv", sep='\t', stringsAsFactors = F, header=F)
alt_all=readMM("data-raw//P5931//scDNA//alt_all_sub.mtx")
ref_all=readMM("data-raw//P5931//scDNA//ref_all_sub.mtx")
var_all=read.table("data-raw//P5931//scDNA//var_all_sub.vcf", header = F, sep='\t', stringsAsFactors = F)

# bin by cell matrices for tumor and normal for segmentation
raw_counts=read.table("data-raw//P5931//scDNA//tumor_sub.txt", sep='\t', header=T, row.names = 1,stringsAsFactors = F)
ref_counts=read.table("data-raw//P5931//scDNA//normal_sub.txt", sep='\t', header=T, row.names = 1,stringsAsFactors = F)

#Step 1: Creating an Alleloscope object for analysis
Obj=Createobj(alt_all =alt_all, ref_all = ref_all, var_all = var_all ,samplename='P5931', genome_assembly="GRCh38", dir_path=dir_path, barcodes=barcodes, size=size, assay='scDNAseq')

Obj_filtered=Matrix_filter(Obj=Obj, cell_filter=10, SNP_filter=10, min_vaf = 0, max_vaf = 1, centro=centromere.GRCh38, telo=telomere.GRCh38) 
#in the line above, 2nd-stage estimation uses much more lenient filters

#Step 2 is skipped in the 2nd-stage estimation, because chromosomes are
#used as segments

#Step 3: Estimate major haplotype proportion for each chromosome
set.seed(2021)
Obj_filtered=Est_regions(Obj_filtered = Obj_filtered, max_nSNP = 30000, plot_stat = T, cont=T)

# Recommend max_nSNP <50000
# Regions without allelic imbalence do not coverge (Reach the max number of iterations.)

#Step 4: Identify normal cells and diploid regions
Obj_filtered=Select_normal(Obj_filtered = Obj_filtered, raw_counts=raw_counts, plot_theta = TRUE)

# add "select_normal" list to the Obj_filtered object. 
# The list includes barcodes for the normal cells and some candidate "normal regions"

print(Obj_filtered$select_normal$region_normal)
Obj_filtered$ref=Obj_filtered$select_normal$region_normal[1] # choose one normal region

#Step 5: Genotype each cell in each region (tumor only)
Obj_filtered=Genotype_value(Obj_filtered = Obj_filtered, type='tumor', raw_counts=raw_counts)  # for tumor

Obj_filtered=Genotype(Obj_filtered = Obj_filtered)

#Optional: Improved estimation based on the abnormal cells
#(this is second-round estimation)
sub_cells=rownames(Obj_filtered$genotypes[which(Obj_filtered$genotypes[,2]!=4),]) #select abnormal cells (4 represents diploid)

Obj_filtered=Est_regions(Obj_filtered = Obj_filtered, max_nSNP = 50000, plot_stat = T, sub_cells =sub_cells, sub_region = '3', max_iter = 100 )

# theta_hat values are updated in the rds_list for the selected region.

#genotype each cell in each region again
# Estimate genotype values
Obj_filtered=Genotype_value(Obj_filtered = Obj_filtered, type='tumor', raw_counts=raw_counts)

# Genotype each cell
Obj_filtered=Genotype(Obj_filtered = Obj_filtered, plot_path="./samples/P5931/scDNA/output/plots/gtype_scatter_updated.pdf")

#Step 6: Construct lineage structure using genotypes for each cell across all regions
linplot=Lineage_plot(Obj_filtered = Obj_filtered, nSNP = 2000,  nclust = 3)