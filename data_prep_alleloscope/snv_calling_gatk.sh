#!/bin/bash
#SBATCH --error="gatk"%J".err"
#SBATCH --output="gatk"%J".out"
#SBATCH -J gatk
#SBATCH -p fat 
#SBATCH --mem-per-cpu 80000

source /home/dhlushchenko/miniconda3/etc/profile.d/conda.sh
conda activate thesis-env
name=SNU601

mkdir "./bam/"
mkdir "./vcf/"
mkdir "./bam/"$name
mkdir "./vcf/"$name

java -jar ../../tools/picard/picard.jar SortSam \
      I=possorted_new_header.bam \
      O="./bam/"$name"/sorted_bam.bam" \
      SORT_ORDER=coordinate

#rm.dup
java -jar ../../tools/picard/picard.jar MarkDuplicates \
I= "./bam/"$name"/sorted_bam.bam" \
O= "./bam/"$name"/possorted_rmdup.bam" \
M= "./bam/"$name"/possorted_rmdup_marked_dup_metics.txt" \
REMOVE_DUPLICATES=true \
VALIDATION_STRINGENCY=LENIENT

#sort
module load ngs/samtools/1.9
samtools sort -o "./bam/"$name"/possorted_rmdup_sort.bam" "./bam/"$name"/possorted_rmdup.bam"

# build index
samtools index "./bam/"$name"/possorted_rmdup_sort.bam" -o "./bam/"$name"/possorted_rmdup_sort.bam.bai"

# mutation calling
#gatk Mutect2 \
#-R ~/genome.fa \
#-I "./bam/"$name"/possorted_rmdup_sort.bam" \
#-O "./vcf/"$name"/possorted_rmdup_sort.m2.vcf.gz" \
#-tumor sample_name \
#--max-population-af 1 \
#--dont-use-soft-clipped-bases 


gatk --java-options "-Xmx4g" HaplotypeCaller \
-R "./hg38.fa" \
-I "./bam/"$name"/possorted_rmdup_sort.bam" \
-O "./vcf/"$name"/possorted_rmdup_sort.hc.vcf.gz"

## decompress file
gunzip "./vcf/"$name"/possorted_rmdup_sort.hc.vcf.gz"