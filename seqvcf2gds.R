library(parallel)
library(SeqArray)
library(SeqVarTools)
library(stringr)

# start cluster
hosts <- strsplit(Sys.getenv("LSB_HOSTS"), " ")[[1]]
if (length(hosts) > 1) {
    if (length(hosts) > 65) {
      warning("parallel can't use more than 64 hosts; found ", length(hosts))
    }
    cluster_hosts <- hosts[2:min(length(hosts), 65)]
    cl <<- makeCluster(cluster_hosts, rshcmd="blaunch")
    message("running on a psock cluster with ", length(cl), " hosts")

} else {
    num_cores <- length(mcaffinity())
    cl <<- makeCluster(num_cores - 1)
    message("running on a fork cluster with ", length(cl), " cores")
}

# read in arguments
args <- commandArgs(trailingOnly=TRUE)
vcfShardsDir <- args[[1]]
gdsFile <- args[[2]]

# do the conversion
vcfFiles <- list.files(vcfShardsDir, "part-[0-9]{5}\\.bgz$", full.names=TRUE)
seqVCF2GDS(vcfFiles, gdsFile, parallel=cl)

# create a nice variant id
f <- seqOpen(gdsFile)
chr <- seqGetData(f, "chromosome")
pos <- seqGetData(f, "position")
alleles <- seqGetData(f, "allele")
seqClose(f)
alleles <- str_replace(alleles, ",", ":")
variant_id <- str_glue("chr{chr}:{pos}:{alleles}")
setVariantID(file.path(gdsShardsDir, gdsFile), variant_id)

stopCluster(cl)

NULL