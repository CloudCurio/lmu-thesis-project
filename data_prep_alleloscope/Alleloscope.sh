#!/bin/bash
#SBATCH --error="alleloscope"%J".err"
#SBATCH --output="alleloscope"%J".out"
#SBATCH -J alleloscope
#SBATCH -p slim24 
#SBATCH --mem-per-cpu 64000 

source /home/dhlushchenko/miniconda3/etc/profile.d/conda.sh
conda activate alleloscope

Rscript ../../lmu-thesis-project/data_prep_alleloscope/Alleloscope.R