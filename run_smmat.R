args <- commandArgs(trailingOnly = F)  
scriptPath <- normalizePath(sub("^--file=", "", args[grep("^--file=", args)]))
source(file.path(dirname(scriptPath), "seqarray_genesis.R"))

args <- commandArgs(trailingOnly=TRUE)
gds_filename <- args[[1]]
null_model_filename <- args[[2]]
model_name <- args[[3]]

run_smmat(gds_filename, null_model_filename, model_name)
