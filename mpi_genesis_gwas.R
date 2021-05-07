library(doMPI, quietly=TRUE, warn.conflicts=FALSE)
library(SeqArray, quietly=TRUE, warn.conflicts=FALSE)
library(SeqVarTools, quietly=TRUE, warn.conflicts=FALSE)
library(GENESIS, quietly=TRUE, warn.conflicts=FALSE)
library(tidyverse, quietly=TRUE, warn.conflicts=FALSE)
cl <- startMPIcluster()
registerDoMPI(cl)

args <- commandArgs(trailingOnly=TRUE)
file_prefix <- args[[1]]
endpoint <- args[[2]]

gdsShardsDir <- paste(file_prefix, "GWAS_filtered", "shards", "seq", "gds", sep=".")

gdsFiles <- list.files(gdsShardsDir, "part-[0-9]{5}\\.seq\\.gds$")
gdsChunkSize <- length(gdsFiles) %/% clusterSize(cl)

nullmod <- readRDS(paste(file_prefix, endpoint, "null", "RDS", sep="."))

# sinkWorkerOutput("/dev/null")
assoc <- foreach(gdsPath=gdsFiles,
                 .combine=rbind,
                 .packages=c("GENESIS", "SeqArray", "SeqVarTools", "tibble", "dplyr"),
                 .options.mpi=list(chunkSize=gdsChunkSize)) %dopar% {
    seqFile <- seqOpen(file.path(gdsShardsDir, gdsPath))
    if (length(seqGetData(seqFile, "variant.id")) == 0) {
        seqClose(seqFile)
        return(NULL)
    } else{
        seqData <- SeqVarData(seqFile)
        seqIter <- SeqVarBlockIterator(seqData)
        assoc_subdf <- assocTestSingle(seqIter,
                        null.model=nullmod,
                        sparse=FALSE, genome.build="hg38")
        seqResetFilter(seqData)
        alleles <- tibble(variant.id=seqGetData(seqData, "variant.id"),
                          effect.allele=altChar(seqData),
                          other.allele=refChar(seqData))
        assoc_subdf <- left_join(assoc_subdf, alleles)
        seqClose(seqFile)
        return(assoc_subdf)
    }
}

write_delim(assoc, paste(endpoint, "GENESIS", "assoc", "txt", sep="."))

# library(qqman)
# png(paste(endpoint, "GENESIS", "qq", "png", sep="."))
# qq(assoc[,Score.pval])
# dev.off()
#
# assoc[,chr_numeric := as.numeric(factor(chr, levels=c(1:22, "X", "Y")))]
# png(paste(endpoint, "GENESIS", "manhattan", "png", sep="."))
# manhattan(assoc, chr="chr_numeric", bp="pos", p="Score.pval", snp="variant.id", chrlabs=c(1:22,"X","Y"), col=c("blue4","orange3"),ylim=c(0,10))
# dev.off()

closeCluster(cl)
mpi.quit()
NULL
