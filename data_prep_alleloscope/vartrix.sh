#!/bin/bash
#SBATCH --error="vartrix"%J".err"
#SBATCH --output="vartrix"%J".out"
#SBATCH -J gatk
#SBATCH -p slim18 
#SBATCH --mem-per-cpu 64000 

# --------------------------------------------------
# Running VarTrix
# --------------------------------------------------
# Geneate SNP by cell matrix from the bam file
# Inputs: vcf file, bam file, fasta file, barcodes.tsv 
# VarTrix only allows parsing one chromosome in each run. 
# This (python) function can help to generate separate vcf files for each chromosome: 
# https://github.com/seasoncloud/Basic_CNV_SNV/blob/main/scripts/vcf_sep_gatk.py
# Install VarTrix (https://github.com/10XGenomics/vartrix)

mkdir chr

for chr in {1..22};
do
	mkdir chr/chr"$chr"_matrix
	../../tools/VarTrix/vartrix_linux -v "./vcf/SNU601/chr"$chr".vcf" -b ./bam/SNU601/possorted_rmdup_sort.bam -f ./hg38.fa -c ./barcodes.tsv.gz -s "coverage"	
	mv out_matrix.mtx chr/chr"$chr"_matrix/out_matrix.mtx
	mv ref_matrix.mtx chr/chr"$chr"_matrix/ref_matrix.mtx
done

# After the matrices are generated for each chromosome, separated matrices can be combined using the Cbn_matrix.R script.
