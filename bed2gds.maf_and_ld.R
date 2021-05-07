library(SNPRelate)

sample_file_basename <- "625_Samples.cohort.QC_filtered.sample_matched"
bedfile <- paste(sample_file_basename, "maf_and_ld", "plink_recoded", "bed", sep=".")
bimfile <- paste(sample_file_basename, "maf_and_ld", "plink_recoded", "bim", sep=".")
famfile <- paste(sample_file_basename, "maf_and_ld", "plink_recoded", "fam", sep=".")
gdsfile <- paste(sample_file_basename, "maf_and_ld", "gds", sep=".")

# convert PLINK to GDS
snpgdsBED2GDS(bed.fn=bedfile,
              bim.fn=bimfile,
              fam.fn=famfile,
              out.gdsfn=gdsfile)

