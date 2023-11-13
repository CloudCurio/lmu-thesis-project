#!/bin/bash
#SBATCH --error="epianeufinder"%J".err"
#SBATCH --output="epianeufinder"%J".out"
#SBATCH -J epianeufinder
#SBATCH -p slim16 
#SBATCH --mem-per-cpu 64000 

source /home/kschmid/miniconda3/etc/profile.d/conda.sh
conda activate alleloscopev2-r41

Rscript ../../lmu-thesis-project/data_prep_alleloscope/epiAneufinder_run.R