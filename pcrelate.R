library(SNPRelate)
library(GWASTools)
library(GENESIS)

sample_file_basename <- "625_Samples.cohort.QC_filtered.sample_matched"
gdsfile <- paste(sample_file_basename, "maf_and_ld", "snp", "gds", sep=".")
in_rdsfile <- paste(sample_file_basename, "PCAir", "RDS", sep=".")
out_rdsfile <- paste(sample_file_basename, "PCRelate", "RDS", sep=".")

pcair_result <- readRDS(in_rdsfile)
geno <- GdsGenotypeReader(filename=gdsfile)
genoData <- GenotypeData(geno)
genoIterator <- GenotypeBlockIterator(genoData)
pcrel_result <- pcrelate(genoIterator, pcs=pcair_result$vectors[,1:4], training.set=pcair_result$unrels)

saveRDS(pcrel_result, file=out_rdsfile)
