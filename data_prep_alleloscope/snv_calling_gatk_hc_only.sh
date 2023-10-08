#!/bin/bash
#SBATCH --error=%J.err
#SBATCH --output=%J.out
#SBATCH -J gatk
#SBATCH -p slim18 
#SBATCH --mem-per-cpu 64000 

source /home/dhlushchenko/miniconda3/etc/profile.d/conda.sh
conda activate thesis-env
name=SNU601

module load ngs/samtools/1.9

gatk --java-options "-Xmx4g" HaplotypeCaller \
-R "./hg38.fa" \
-I "./bam/"$name"/possorted_rmdup_sort.bam" \
-O "./vcf/"$name"/possorted_rmdup_sort.hc.vcf.gz"

## decompress file
gunzip "./vcf/"$name"/possorted_rmdup_sort.hc.vcf.gz"