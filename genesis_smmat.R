library(doMPI, quietly=TRUE)
library(SeqArray, quietly=TRUE)
library(SeqVarTools, quietly=TRUE)
library(GENESIS, quietly=TRUE)
cl <- startMPIcluster()
registerDoMPI(cl)

args <- commandArgs(trailingOnly=TRUE)
endpoint <- args[[1]]
file_prefix <- "625_Samples.cohort.QC_filtered.sample_matched"

gdsShardsDir <- paste(file_prefix, "GWAS_filtered", "shards", "seq", "gds", sep=".")

gdsFiles <- list.files(gdsShardsDir, "part-[0-9]{5}\\.seq\\.gds$")
seqMerge(gdsFiles, paste(file_prefix, "GWAS_filtered", "seq", "gds", sep="."))

nullmod <- readRDS(paste(file_prefix, endpoint, "null", "RDS", sep="."))    

seqFile <- seqOpen(file.path(gdsShardsDir, gdsPath))
        seqIter <- SeqVarBlockIterator(SeqVarData(seqFile))
        assoc_subdf <- assocTestSingle(seqIter,
                        null.model=nullmod,
                        sparse=FALSE, genome.build="hg38")
        seqClose(seqFile)
        return(assoc_subdf)
    }
}

library(data.table)
assoc <- as.data.table(assoc)
write.table(assoc, paste(endpoint, "GENESIS", "assoc", "txt", sep="."))

library(qqman)
png(paste(endpoint, "GENESIS", "qq", "png", sep="."))
qq(assoc[,Score.pval])
dev.off()

assoc[,chr_numeric := as.numeric(factor(chr, levels=c(1:22, "X", "Y")))]
png(paste(endpoint, "GENESIS", "manhattan", "png", sep="."))
manhattan(assoc, chr="chr_numeric", bp="pos", p="Score.pval", snp="variant.id", chrlabs=c(1:22,"X","Y"), col=c("blue4","orange3"),ylim=c(0,10))
dev.off()

closeCluster(cl)
mpi.quit()
NULL
