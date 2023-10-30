library(Signac)

sizes_file <- read.table("./hg38/chrom.sizes.txt")
chrom_sizes <- sizes_file$V2
names(chrom_sizes) <- sizes_file$V1

barcodes<-read.table("barcodes.tsv.gz")
print("line 8 ok")
tmp <- CreateFragmentObject(path="fragments.tsv.gz",
                            cells=barcodes$V1)
print("line 11 ok")
mtx <- GenomeBinMatrix(fragments = list(tmp),
                       genome = chrom_sizes,
                       binsize = 100000)
print("line 15 ok")
mtx <- as.matrix(mtx)
write.table(mtx, file = "100k_fragments_sub.txt")