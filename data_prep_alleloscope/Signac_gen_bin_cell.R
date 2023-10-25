library(Signac)

sizes_file <- read.table("./hg38/chrom.sizes.txt")
chrom_sizes <- sizes_file$V2
names(chrom_sizes) <- sizes_file$V1

GenomeBinMatrix(fragments = "./fragments.tsv",
                genome = chrom_sizes,
                binsize = 100000)