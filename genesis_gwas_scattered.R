library(SeqArray, quietly=TRUE, warn.conflicts=FALSE)
library(SeqVarTools, quietly=TRUE, warn.conflicts=FALSE)
library(GENESIS, quietly=TRUE, warn.conflicts=FALSE)
library(tibble, quietly=TRUE, warn.conflicts=FALSE)
library(dplyr, quietly=TRUE, warn.conflicts=FALSE)
library(fs, quietly=TRUE, warn.conflicts=FALSE)

doGwasShard <- function (seqFile, nullmod) {
    seqData <- SeqVarData(seqFile)
    seqIter <- SeqVarBlockIterator(seqData)
    assoc_df <- assocTestSingle(seqIter,
                                null.model=nullmod,
                                sparse=FALSE, genome.build="hg38")
    seqResetFilter(seqData)

    alleles <- tibble(variant.id=seqGetData(seqData, "variant.id"),
                      effect.allele=altChar(seqData),
                      other.allele=refChar(seqData))
    assoc_df <- left_join(assoc_df, alleles)
    if (nullmod$family$family == "binomial") {
        # calculate cases and controls
        seqSetFilter(seqData, variant.id=assoc_df$variant.id,
                     sample.id=row.names(nullmod$outcome)[nullmod$outcome == 0])
        n.controls <- colSums(!is.na(getGenotype(seqData)))
        seqSetFilter(seqData, variant.id=assoc_df$variant.id,
                     sample.id=row.names(nullmod$outcome)[nullmod$outcome == 1])
        n.cases <- colSums(!is.na(getGenotype(seqData)))
        assoc_df <- add_column(assoc_df, n.cases=n.cases, n.controls=n.controls)
    }
    assoc_df
}

seqFile <- seqOpen(snakemake@input[["gds_shard"]])

if (length(seqGetData(seqFile, "variant.id")) > 0) {
    file_touch(snakemake@output[[1]])
} else {
    nullmod <- readRDS(snakemake@input[["null_model"]])
    assoc_shard <- doGwasShard(seqFile, nullmod)
    write_delim(assoc_shard, snakemake@output[[1]])
}
seqClose(seqFile)
NULL
