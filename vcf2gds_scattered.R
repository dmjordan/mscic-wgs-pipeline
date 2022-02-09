library(SeqArray, quietly=T)
library(stringr, quietly=T)

dir.create(dirname(snakemake@output[[1]]), showWarnings=FALSE)
seqVCF2GDS(snakemake@input[[1]], snakemake@output[[1]], ignore.chr.prefix="")

# create a nice variant id
f <- seqOpen(snakemake@output[[1]])
chr <- seqGetData(f, "chromosome")
pos <- seqGetData(f, "position")
alleles <- seqGetData(f, "allele")
seqClose(f)
alleles <- str_replace(alleles, ",", ":")
variant_id <- str_glue("{chr}:{pos}:{alleles}")
setVariantID(snakemake@output[[1]], variant_id)

NULL
