#Data downloaded from SRA
module load ngs/sratoolkit/2.10.0
fastq-dump --split-files --origfmt --gzip SRR12995620

#Renaming files to run cellranger
mv SRR12995620_1.fastq.gz SRR12995620_S1_L001_R1_001.fastq.gz
mv SRR12995620_2.fastq.gz SRR12995620_S1_L001_R2_001.fastq.gz