library(SNPRelate)
library(GWASTools)
library(GENESIS)
library(data.table)
library(ggplot2)

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


# load KING result
KINGmat <- read.table("625_Samples.cohort.KING.txt", header=TRUE)
KINGmat <- reshape(KINGmat, idvar="DNA_sample_1", timevar="DNA_sample_2", direction="wide")
rownames(KINGmat) <- KINGmat$DNA_sample_1
KINGmat <- KINGmat[-1]
colnames(KINGmat) <- substring(colnames(KINGmat), 5)
KINGmat <- as.matrix(KINGmat)

# run PC-AIR
geno <- GdsGenotypeReader(filename=gdsfile)
genoData <- GenotypeData(geno)

pcair_result <- pcair(genoData, kinobj=KINGmat, divobj=KINGmat)

# Plot first 10 PCs
pc_table <- data.table(pcair_result$vectors)
names(pc_table) <- paste("PC", 1:length(pc_table), sep="")
pc_table[,Subject_ID := sapply(strsplit(pcair_result$sample.id, "T"), function(x) { x[1] }) ]
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

# run PCRelate

genoIterator <- GenotypeBlockIterator(genoData)
pcrel_result <- pcrelate(genoIterator, pcs=pcair_result$vectors[,1:4], training.set=pcair_result$unrels)

# plot PCRelate
pdf(paste(sample_file_basename, ".PCRelate.pdf", sep=""))
plot(pcrel_result$kinBtwn$k0, pcrel_result$kinBtwn$kin, xlab="k0", ylab="kinship")
dev.off()

save(c("pcair_result", "pcrel_result"), file=paste(sample_file_basename, "GENESIS_PCA.Rdata", sep="."))


