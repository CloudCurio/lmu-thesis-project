#!/bin/bash
#SBATCH --error="alleloscope"%J".err"
#SBATCH --output="alleloscope"%J".out"
#SBATCH -J alleloscope
#SBATCH -p slim16 
#SBATCH --mem-per-cpu 64000 

source /home/kschmid/miniconda3/etc/profile.d/conda.sh
conda activate alleloscopev2-r41

Rscript ../../lmu-thesis-project/data_prep_alleloscope/Alleloscope_batch_run.R