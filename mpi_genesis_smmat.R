suppressPackageStartupMessages({
  library(doMPI)
  library(SeqArray)
  library(SeqVarTools)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(tidyverse)
})

cl <- startMPIcluster()
registerDoMPI(cl)
cluster_size <- clusterSize(cl)

geneRanges <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
splitGeneRanges <- split(geneRanges, 1:length(geneRanges) %% cluster_size)

args <- commandArgs(trailingOnly=TRUE)
gds_filename <- args[[1]]
null_model_filename <- args[[2]]
model_name <- args[[3]]

nullmod <- readRDS(null_model_filename)
assoc <- foreach(gRanges=splitGeneRanges,
                 .combine=rbind,
                 .packages=c("GENESIS", "SeqArray", "SeqVarTools", "tibble", "magrittr")) %dopar% {
        gdsFile <- seqOpen(gds_filename, allow.duplicate=TRUE)
        tryCatch({
                 seqData <- SeqVarData(gdsFile)
                 seqIter <- SeqVarRangeIterator(seqData, gRanges)

                 assoc <- assocTestAggregate(seqIter, nullmod, sparse=FALSE, test="SMMAT")
                 as_tibble(assoc$results, rownames="geneid") %>% add_column(chr=as.factor(seqnames(gRanges)),
                                                                            start=start(gRanges),
                                                                            end=end(gRanges), .after="geneid")
            }, finally=seqClose(gdsFile))
}

write_delim(assoc, paste(model_name, 'GENESIS', 'SMMAT', 'assoc', 'txt', sep='.'))

closeCluster(cl)
mpi.quit()
NULL
