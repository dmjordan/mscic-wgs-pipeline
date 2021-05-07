library(doMPI, quietly=T)
library(SeqArray, quietly=T)
library(stringr, quietly=T)
cl <- startMPIcluster()
registerDoMPI(cl)

vcfShardsDir <- commandArgs(trailingOnly=TRUE)
gdsShardsDir <- sub("\\.vcf\\.bgz", ".seq.gds", vcfShardsDir)
dir.create(gdsShardsDir, showWarnings=FALSE)
vcfFiles <- list.files(vcfShardsDir, "part-[0-9]{5}\\.bgz$")
vcfChunkSize <- length(vcfFiles) %/% clusterSize(cl)

foreach(vcfFile=vcfFiles,
        .packages=c("SeqArray", "SeqVarTools", "stringr"),
        .options.mpi=list(chunkSize=vcfChunkSize),
        .errorhandling="stop") %dopar% {
    gdsFile <- sub("\\.bgz$", ".seq.gds", vcfFile)
    seqVCF2GDS(file.path(vcfShardsDir, vcfFile), file.path(gdsShardsDir, gdsFile))
    # create a nice variant id
    f <- seqOpen(file.path(gdsShardsDir, gdsFile))
    #rsid <- seqGetData(f, "annotation/id")
    chr <- seqGetData(f, "chromosome")
    pos <- seqGetData(f, "position")
    alleles <- seqGetData(f, "allele")
    seqClose(f)
    alleles <- str_replace(alleles, ",", ":")
    variant_id <- str_glue("chr{chr}:{pos}:{alleles}")
    setVariantID(file.path(gdsShardsDir, gdsFile), variant_id)  # ifelse(rsid != "", rsid, variant_id))
}

closeCluster(cl)
mpi.quit()
NULL
