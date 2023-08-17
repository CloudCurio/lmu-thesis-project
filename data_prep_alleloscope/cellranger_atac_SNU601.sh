#!/bin/bash
#SBATCH --error=%J.err
#SBATCH --output=%J.out
#SBATCH -J cellranger
#SBATCH -p slim18 
#SBATCH --mem-per-cpu 64000 

cellranger-atac count --reference=/work/project/ladcol_014/thesis_cnvcalling/data/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 --fastqs=/work/project/ladcol_014/thesis_cnvcalling/data/SNU601_scATACseq --sample=SRR12995620 --localcores 8 --localmem 64 --id=SNU601_atac



