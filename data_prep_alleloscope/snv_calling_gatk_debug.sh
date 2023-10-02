#!/bin/bash
#SBATCH --error=%J.err
#SBATCH --output=%J.out
#SBATCH -J gatk
#SBATCH -p slim18 
#SBATCH --mem-per-cpu 64000 

source /home/dhlushchenko/miniconda3/etc/profile.d/conda.sh
conda activate thesis-env
name=SNU601

mkdir "./bam_debug/"$name
mkdir "./vcf_debug/"$name

java -jar ../../tools/picard/picard.jar SortSam \
      I=possorted_new_header.bam \
      O=bam_debug/sorted_bam.bam \
      SORT_ORDER=coordinate

#rm.dup
java -jar ../../tools/picard/picard.jar MarkDuplicates \
I= "bam_debug/sorted_bam.bam" \
O= "./bam_debug/"$name"/possorted_rmdup.bam" \
M= "./bam_debug/"$name"/possorted_rmdup_marked_dup_metics.txt" \
REMOVE_DUPLICATES=true \
VALIDATION_STRINGENCY=LENIENT

#sort
module load ngs/samtools/1.9
samtools sort -o "./bam_debug/"$name"/possorted_rmdup_sort.bam" "./bam_debug/"$name"/possorted_rmdup.bam"

# build index
java -jar ../../tools/picard/picard.jar BuildBamIndex \
I= "./bam_debug/"$name"/possorted_rmdup_sort.bam" \
VALIDATION_STRINGENCY=LENIENT


# mutation calling
#gatk Mutect2 \
#-R ~/genome.fa \
#-I "./bam/"$name"/possorted_rmdup_sort.bam" \
#-O "./vcf/"$name"/possorted_rmdup_sort.m2.vcf.gz" \
#-tumor sample_name \
#--max-population-af 1 \
#--dont-use-soft-clipped-bases 


gatk --java-options "-Xmx4g" HaplotypeCaller \
-R "/work/project/ladcol_004/epiAneufinder/ReferenceBWA/hg38_genome.fa" \
-I "./bam_debug/"$name"/possorted_rmdup_sort.bam" \
-O "./vcf_debug/"$name"/possorted_rmdup_sort.hc.vcf.gz"

## decompress file
gunzip "./vcf_debug/"$name"/possorted_rmdup_sort.hc.vcf.gz"