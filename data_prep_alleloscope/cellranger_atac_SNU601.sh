#!/bin/bash
#SBATCH --error=%J.err
#SBATCH --output=%J.out
#SBATCH -J cellranger
#SBATCH -p slim18 
#SBATCH --mem-per-cpu 64000 

/work/project/ladcol_005/tools/cellranger-atac-2.1.0/bin/cellranger-atac count --reference=/work/project/ladcol_005/genomes/hg38/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 --fastqs=/work/project/ladcol_014/thesis_cnvcalling/data/SNU601_scATACseq --sample=SRR12995620 --localcores 8 --localmem 64 --id=SNU601_atac



