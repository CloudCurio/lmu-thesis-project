# Identifying copy number alterations from scATAC-seq data using allelic information (Dmytro Hlushchenko, LMU, 2024)
## Preface
This repository contains the code written for my Master's thesis and describes a computational tool for Copy Number Variation (CNV) prediction from scATAC-seq data, as well as all the necessary preparatory steps. This README describes the workflow of the tool, as well as the purpose of individual scripts. For a complete description of the tool's performance and creation, please refer to the full text of the thesis:
[Hlushchenko, 2024.pdf](https://github.com/CloudCurio/lmu-thesis-project/files/14616168/Hlushchenko.2024.pdf)
## Graphical abstract
<p align="center">
 <img src="https://github.com/CloudCurio/lmu-thesis-project/assets/66057046/1de4c047-079d-492d-a6ce-c341726da385" width="720">
</p>

## Abstract
With the increasing number and improving quality of sequencing technologies in the recent years, many spheres of biology were revolutionized, enabling new discoveries at unprecedented definition, and changing our understanding of many biological subjects, including cancer research. This became especially clear with the advent of single-cell sequencing technologies, allowing to observe processes and their changes on the scale of an individual cell, removing potential heterogeneity of a bulk sample. As a result, the mechanism of oncogenesis and its effect on cell’s biology were greatly illuminated, and new contributing factors were identified. One of such factors are copy number variations (CNVs) – changes in copy number of certain genomic fragments. CNVs were found to be prolific in cancerous cells and recognized as major cancer drivers. Consequently, great effort is directed at developing methods for their detection and characterization.

However, due to the difficulty of performing single cell whole genome sequencing (scWGS) experiments and lack of available datasets, other sequencing technologies are also explored for this task. Furthermore, other technologies, like scRNA-seq or scATAC-seq, provide new facets of information, unavailable in scWGS, e.g., the state of the cell, making them valuable objects to analyze. At the moment, many methods were developed for CNV prediction from single cell RNA sequencing (scRNA-seq) assays, but single cell Assay for Transposase-Accessible Chromatin (scATAC-seq) data is less studied, with only a couple of CNV callers using it as input. However, scATAC-seq has a major advantage over scRNA-seq, providing better coverage of the entire genome instead of only genes. Current tools for scATAC-seq analysis use read count data as input and show good performance, but previous experience with scRNA-seq shows merit in supplementing read counts with allelic frequency (AF) information too. Unfortunately, calling AF for scATAC-seq is complicated due to the sparsity of data.

In this thesis project, we explore whether CNV prediction from scATAC-seq data could be improved with the introduction of AF information. To this end, we predict AF values for the gastric cancer cell line SNU601 using the Alleloscope tool and design multiple dependent mixture Hidden Markov Models (HMM) using the depmixS4 R package for per cell CNV identification on 500kbp-1mbp definition. As a results, 2-state and 3-state models were built with the use of counts only and counts+AF inputs at an individual cell definition. Due to the sparsity of estimated AF values, data imputation was also attempted for the regions with no prediction available. After analyzing the model’s performance in comparison with scWGS ground truth, we report a mean across cells accuracy of 0.531 with the 3-state HMM using counts and imputed AF values. It was also compared with the previously developed scATAC-seq CNV-caller epiAneufinder (which uses read count inputs and is based on the segmentation algorithm), showing improved predictions for gains, but worse predictions for loss. Although the current model is not suitable for practical applications, this attempt still shows potential of using AF information to improve CNV prediction in future research.
## Project structure
Overall workflow of the tool can be described as follows:
1. Get a scATAC-seq dataset to analyze and prepare the read counts data as usual
2. Prepare data for Alleloscope analysis
3. Predict Allele Frequency (AF) information using Alleloscope
4. Use read counts and AF information to train an HMM and predict CNVs

Scripts for the project are split into 3 folders: "epiAneufinder_scripts", "data_prep_alleloscope", and "analysis_scripts". Folders "debug_function_versions" and "tutorial_scripts" are technical and should be ignored.

## "epiAneufinder_scripts" folder
This folder contains scripts needed to run epiAneufinder, outputs of which are used throughout the analysis
- epiAneufinder_run.R - runs epiAneufinder with selected parameters;
- epiAneufinder_run.sh - bash script to run the previous script on the server with the benefits of the process manager.

## "data_prep_alleloscope" folder
This folder covers steps 1-3 of the workflow above. The inputs required are:
- hg38.fa - reference genome sequence in the FASTA format;
- possorted_new_header.bam - position-sorted BAM file of the scATAC assay. Make sure that the header of the BAM is correct, since it creates errors with the tools used;
- barcodes.tsv - a .tsv file containing all cellular barcodes in the scATAC assay;
- chrom_sizes.txt - a text file containing sizes of all chromosomes in the data;
- fragments.tsv.gz - a compressed TSV file listing scATAC fragments in the data.

### Initial data prep:
- dict_creation.sh - constructs a .dict file for the provided reference genome .fa data;
- download_SNU601.sh - downloads the SNU601 WGS data in the .fastq format;
- cellranger_atac_SNU601.sh - downloads the scATAC-seq data in the .fastq format;
- filter_bam.sh - filters the bam file for autosomal chromosomes (taken from Katharina Schmid).

### Alleloscope data preparation
- snv_calling_gatk.sh - gets the .vcf file with SNP info, using the reference genome .fa file and the .bam file of the scATAC-data;
- snv_calling_gatk_hc_only.sh - a shorter version of the above script, only including the HaplotypeCaller command. Useful for saving time rerunning the script in case of HaplotypeCaller errors;
- vcf_sep_gatk.py - splits the .vcf output of the snv_calling_gatk.sh by chromosome, producing a number of smaller .vcf files;
- vartrix.sh - gets SNP-by-cell matrices (alt_all.rds, ref_all.rds, var_all.rds) for each chromosome using per-chromosome .vcf files;
- Cbn_matrix.R - an R function to merge per-chromosome SNP-by-cell matrices into whole-genome forms. **Technical file, sourced in the next script**;
- Cbn_matrix_server_paths.R - uses the function above and allows to set the paths to inputs/outputs for the server;
- Signac_gen_bin_cell.R - creates the bin-by-cell matric using Signac, using chrom_sizes.txt, barcodes.tsv.gz, and fragments.tsv.gz files;
- Gen_bin_cell_atac.R - old script for the same task as Signac_gen_bin_cell.R. **Not used anymore, ignore this file**;
- generate_seg_table.R - an R function that calculates pseudobulk CNV predictions using the epiAneufinder's segmentation table (results_table.tsv) and writes chromosome labels in character form instead of integer. **A function, used in the next step, don't run individually**;
- merge_segments.R - a standalone alternative to the previous function. Does the same thing, but can also merge neighbouring segments with the same classes. **Was used earlier, but not anymore**.

### Running Alleloscope
- Alleloscope_fun.R - runs Alleloscope using previously created files (see documentation at the start of the script). Intended for per-chromosome batch runs, sourced by the next script;
- Alleloscope_batch_run.R - runs the Alleloscope_fun.R script on each chromosome individually. Uses the generate_seg_table.R function script from before to get the segmentation results. After running Alleloscope on a chromosome, immediately uses the function from alleloscope_summary.R (see below, in the "analysis_scripts":"Analyzing Alleloscope results" section) function on the outputs to avoid flooding the output folder with hundreds of files. Called by the bash script below;
- Alleloscope.sh - performs Alleloscope analysis using the above scripts;
- get_seg_table.R - extracts the filtered segmentation table from the Alleloscope output .rds object, stores it in a separate .rds file.

## "analysis_scripts" folder
This folder covers the post-processing of Alleloscope results in the step 3, as well as the construction and testing of HMMs for the step 4 of the workflow. Since the HMM part was coded to run on my personal laptop, make sure to check that all paths are adapted to your preferred folder structure.

### Analyzing Alleloscope results
- alleloscope_summary.R - generates fragment_summary.csv and all_frags.rds files by merging averaged-by-region AF data from Alleloscope and the segmentation table from epiAneufinder. In addition, creates a number of plots describing Alleloscope output. **Intended for use with Alleloscope_batch_run.R**;
- bin_by_cell_theta_summary.R - summarizes Alleloscope outputs and produces the bin-by-cell AF matrix (bin_by_cell_theta.csv) and the SNP_by_region matrix (SNP_counts.csv);
- alleloscope_visualization.R - Prepares input files for the HMM: creates a pseudobulk and per-cell sets of epiAneufinder CNV calls, WGS CNV calls, AF information from Alleloscope, checks missing values percentage in AF data for each cell, removing the cells with less than 50% of the values from the outputs with the "filtered" suffix (as opposed to the "full" suffix, including all cells). AF imputation is also done in this script. Also generates a number of plots for data evaluation.

### Building and running the HMMs
- distribution_fit_test.R - 
- HMM_2.R, HMM_3.R - build a set of HMM models - one for each cell - with 2 or 3 hidden states respectively. Include a number of settings to set: seed value for replicability, stat_flag for distribution testing, impute_theta_flag for enabling analysis with imputation, and models_to_test, which allows to set the types of input data to test;
- HMM_2_posterior.R, HMM_3_posterior.R - same as above, but with posterior probability thresholds for the loss/gain class;
- plot_emissions_3_state.R - plots emission distributions of a 3-state HMM;
- plot_posteriors_3_state.R - plots posterior probability distribution for the 3-state HMM;
- emission_dist_group_metric_analysis.R - convenience script to compare emission distributions for different cells. **Not needed for analysis**
