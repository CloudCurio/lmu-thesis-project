#!/bin/bash
#SBATCH --error=%J.err
#SBATCH --output=%J.out
#SBATCH -J filter_bam
#SBATCH -p slim16
#SBATCH -c 1 #number of CPUs needed 
#SBATCH --mem-per-cpu=16G 

#Filter samtools for the autosomal chromosomes
module load ngs/samtools/1.9
#samtools view -b /work/project/ladcol_014/thesis_cnvcalling/data/SNU601_scATACseq/bam/SNU601/possorted_rmdup_sort.bam chr{1..22} > /work/project/ladcol_010/possorted_rmdup_sort_filtered.bam

samtools index /work/project/ladcol_010/possorted_rmdup_sort_filtered.bam


