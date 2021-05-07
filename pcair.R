library(SNPRelate)
library(GWASTools)
library(GENESIS)
library(data.table)
library(ggplot2)

sample_file_basename <- "625_Samples.cohort.QC_filtered.sample_matched"
gdsfile <- paste(sample_file_basename, "maf_and_ld", "snp", "gds", sep=".")
rdsfile <- paste(sample_file_basename, "PCAir", "RDS", sep=".")
txtfile <- paste(sample_file_basename, "PCAir", "txt", sep=".")

# load KING result
KINGmat <- fread("625_Samples.cohort.KING.txt", header=TRUE)
# PICR7147T1a and PICR7147T1b are duplicated
KINGmat <- KINGmat[(DNA_sample_1 != "PICR7147T1a") & (DNA_sample_2 != "PICR7147T1a")]
KINGmat[,DNA_sample_1 := tstrsplit(DNA_sample_1,'T')[1]]
KINGmat[,DNA_sample_2 := tstrsplit(DNA_sample_2,'T')[1]]
KINGmat <- dcast(KINGmat, DNA_sample_1~DNA_sample_2)
KINGmat <- as.matrix(KINGmat, rownames="DNA_sample_1")

# run PC-AIR
geno <- GdsGenotypeReader(filename=gdsfile)
genoData <- GenotypeData(geno)

pcair_result <- pcair(genoData, kinobj=KINGmat, divobj=KINGmat)

# Plot first 10 PCs
pc_table <- data.table(pcair_result$vectors)
names(pc_table) <- paste("PC", 1:length(pc_table), sep="")
pc_table[,Subject_ID := pcair_result$sample.id]
setkey(pc_table, Subject_ID)
clinical_covariates <- as.data.table(readRDS("../../data/covariates/clinical_data_deidentified_allsamples/Biobank_clinical_data_table_by_blood_sample_deidentified_UNCONSENTED.RDS"))
race_table <- clinical_covariates[,.(Race_From_Consent=unique(Race_From_Consent), 
                                     Ethnicity_From_Consent=unique(Ethnicity_From_Consent)),
                                    keyby=Subject_ID]
pc_table <- merge(pc_table, race_table)
for (pc1 in 1:9) {
    for (pc2 in pc1:10) {
        pdf(paste(sample_file_basename, ".PC", pc1, "v", pc2, ".pdf", sep=""))
        ggplot(pc_table, aes_string(x=paste("PC", pc1, sep=""), 
                               y=paste("PC", pc2, sep=""),
                               color="Race_From_Consent"),
                               shape="Ethnicity_From_Consent") + geom_point()
        dev.off()
    }
}

saveRDS(pcair_result, file=rdsfile)
write.table(pcair_result$vectors, file=txtfile)
