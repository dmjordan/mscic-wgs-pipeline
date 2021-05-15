source("seqarray_genesis.R")

start_cluster()

args <- commandArgs(trailingOnly=TRUE)
gds_filename <- args[[1]]
null_model_filename <- args[[2]]
model_name <- args[[3]]

run_smmat(gds_filename, null_model_filename, model_name)