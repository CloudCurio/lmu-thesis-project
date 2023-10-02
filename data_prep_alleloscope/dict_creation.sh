#!/bin/bash
#SBATCH --error=%J.err
#SBATCH --output=%J.out
#SBATCH -J gatk
#SBATCH -p slim18 
#SBATCH --mem-per-cpu 64000 

source /home/dhlushchenko/miniconda3/etc/profile.d/conda.sh
conda activate thesis-env

java -jar /work/project/ladcol_014/thesis_cnvcalling/tools/picard/picard.jar CreateSequenceDictionary \ 
      R=genome.fa \ 
      O=genome.dict