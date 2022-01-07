library("GENESIS")
library("GWASTools")

run_pcrelate <- function(prefix) {
    gdsfile <- paste(prefix, "GWAS_filtered", "LD_pruned", "snp", "gds", sep=".")
    in_rdsfile <- paste(prefix, "PCAir", "RDS", sep=".")
    out_rdsfile <- paste(prefix, "PCRelate", "RDS", sep=".")

    pcair_result <- readRDS(in_rdsfile)
    geno <- GdsGenotypeReader(filename=gdsfile)
    genoData <- GenotypeData(geno)
    genoIterator <- GenotypeBlockIterator(genoData)
    pcrel_result <- pcrelate(genoIterator, pcs=pcair_result$vectors[,1:4], training.set=pcair_result$unrels)

    saveRDS(pcrel_result, file=out_rdsfile)
    pcrel_result
}

run_pcrelate(snakemake@wildcards[["prefix"]])

