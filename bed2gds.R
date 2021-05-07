library(SNPRelate)

args <- commandArgs(trailingOnly=TRUE)

sample_file_basename <- args[1]
bedfile <- paste(sample_file_basename, "bed", sep=".")
bimfile <- paste(sample_file_basename, "bim", sep=".")
famfile <- paste(sample_file_basename, "fam", sep=".")
gdsfile <- paste(sample_file_basename, "snp", "gds", sep=".")

# convert PLINK to GDS
snpgdsBED2GDS(bed.fn=bedfile,
              bim.fn=bimfile,
              fam.fn=famfile,
              out.gdsfn=gdsfile)

