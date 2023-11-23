library(Alleloscope) # load the library

setwd(".//Alleloscope//") # set path to the github folder
dir_path <- ".//samples//SNU601//scDNA//output//"; dir.create(dir_path) # set up output directory

data(centromere.GRCh38)
data(telomere.GRCh38)
size=read.table("data-raw//sizes.cellranger-GRCh38-1.0.0.txt", stringsAsFactors = F)

# SNP by cell matrices for ref and alt alleles
barcodes=read.table("data-raw//SNU601//scDNA//barcodes_sub.tsv", sep='\t', stringsAsFactors = F, header=F)
alt_all=readMM("data-raw//SNU601//scDNA//alt_all_sub.mtx")
ref_all=readMM("data-raw//SNU601//scDNA//ref_all_sub.mtx")
var_all=read.table("data-raw//SNU601//scDNA//var_all_sub.vcf", header = F, sep='\t', stringsAsFactors = F)

# bin by cell matrices for tumor and normal for segmentation
raw_counts=read.table("data-raw//SNU601//scDNA//tumor_sub.txt", sep='\t', header=T, row.names = 1,stringsAsFactors = F)
ref_counts=read.table("data-raw//SNU601//scDNA//normal_sub.txt", sep='\t', header=T, row.names = 1,stringsAsFactors = F) # Normal sample from patient 6198 was used for the cell line.

#Please make sure the SNPs are located on chromosome 1-22.

#after this, all steps may be run with one line:
# Obj_filtered=Rundf_dna(alt_all =alt_all, ref_all = ref_all, var_all = var_all ,
#                        samplename='Sample', genome_assembly="GRCh38", dir_path=dir_path, 
#                        barcodes=barcodes, size=size, assay='scDNAseq',
#                        raw_counts=raw_counts, ref_counts=ref_counts, type='cellline',
#                        cell_filter = 1000, SNP_filter = 20, min_vaf = 0.1, max_vaf = 0.9)

#Step 1: Creating an Alleloscope object for the analysis
Obj=Createobj(alt_all =alt_all, ref_all = ref_all, var_all = var_all, samplename='Sample', genome_assembly="GRCh38", dir_path=dir_path, barcodes=barcodes, size=size, assay='scDNAseq')

Obj_filtered=Matrix_filter(Obj=Obj, cell_filter=1000, SNP_filter=20, min_vaf = 0.1, max_vaf = 0.9, centro=centromere.GRCh38, telo=telomere.GRCh38) 

# suggest setting min_vaf=0.1 and max_vaf=0.9 when SNPs are called in the tumor sample for higher confident SNPs.

#Step 2: Segmentation based on total coverage pooled across cells
Obj_filtered=Segmentation(Obj_filtered=Obj_filtered, 
                          raw_counts=raw_counts, # from matched DNA sequencing (bulk/single)
                          ref_counts=ref_counts, # from matched DNA sequencing (bulk/single)
                          plot_seg = TRUE)

Obj_filtered=Segments_filter(Obj_filtered=Obj_filtered, nSNP=2000)

#Step 3: Estimate cell major haplotype proportion for each region
Obj_filtered=Est_regions(Obj_filtered = Obj_filtered, max_nSNP = 30000, plot_stat = T,cont = TRUE)

# Recommend max_nSNP <50000
# Regions without allelic imbalence do not coverge (Reach the max number of iterations.)

#Step 4: Identify normal cells and diploid regions:
Obj_filtered=Select_normal(Obj_filtered = Obj_filtered, raw_counts=raw_counts, plot_theta = TRUE)

# add "select_normal" list to the Obj_filtered object. 
# The list includes barcodes for the normal cells and some candidate "normal regions"

print(Obj_filtered$select_normal$region_normal)
Obj_filtered$ref=Obj_filtered$select_normal$region_normal[1] # choose one normal region

#Step 5: Genotype each cell in each region
Obj_filtered=Genotype_value(Obj_filtered = Obj_filtered, type='tumor', raw_counts=raw_counts)  # for tumor
Obj_filtered=Genotype_value(Obj_filtered = Obj_filtered, type='cellline', raw_counts=raw_counts, ref_counts = ref_counts ) # for cell line without normal cells in the tumor sample.

Obj_filtered=Genotype(Obj_filtered = Obj_filtered)

#Step 6:
linplot=Lineage_plot(Obj_filtered = Obj_filtered, nSNP = 2000, nclust = 10)