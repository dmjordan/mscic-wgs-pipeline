library(SeqArray)

args <- commandArgs(trailingOnly=TRUE)

sample_file_basename <- args[1]
infile <- paste(sample_file_basename, "snp", "gds", sep=".")
outfile <- paste(sample_file_basename, "seq", "gds", sep=".")


seqSNP2GDS(infile, outfile)
