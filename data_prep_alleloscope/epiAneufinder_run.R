library(epiAneufinder)
epiAneufinder(input="fragments.tsv.gz", #Enter path to your fragments.tsv file or the folder containing bam files
              outdir="..//..//output//epiAneufinder_results", #Path to the directory where results should be written 
              blacklist="hg38-blacklist.v2.bed", #Path to bed file that contains the blacklisted regions of your genome
              windowSize=5e5, 
              genome="BSgenome.Hsapiens.UCSC.hg38", #Substitute with relevant BSgenome
              exclude=c('chrX','chrY','chrM'), 
              reuse.existing=TRUE,
              title_karyo="Karyogram of sample data", 
              ncores=4,
              minFrags=20000,
              minsizeCNV=0,
              k=4,
              plotKaryo=TRUE)