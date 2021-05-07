library(data.table)
library(GENESIS)
library(GWASTools)
library(SNPRelate)
library(parallel)
library(qqman)

args <- commandArgs(trailingOnly=TRUE)
endpoint <- args[[1]]

sample_file_basename <- "625_Samples.cohort.QC_filtered.sample_matched"
gdsfile <- paste(sample_file_basename, "GWAS_filtered", "snp", "gds", sep=".")
rdsfile <- paste(sample_file_basename, "PCRelate", "RDS", sep=".")

cat("loading clinical data\n", file=stderr())
clinical_table <- read.csv("../../data/covariates/clinical_data_deidentified_allsamples/jordad05/625_Samples.cohort.QC_filtered.sample_matched.age_flowcell_PCAir_dmatrix.csv", header=TRUE)
names(clinical_table)[1] <- "scanID"
covars <- c("age", "sex", "age_sex", #"is_hispanic",
               names(clinical_table)[startsWith(names(clinical_table), "pc")],
               #names(clinical_table)[startsWith(names(clinical_table), "race")],
               names(clinical_table)[startsWith(names(clinical_table), "flowcell")])

scan_annot <- ScanAnnotationDataFrame(clinical_table)
clinical_table <- as.data.table(clinical_table)

cat("loading PCRelate\n", file=stderr())
pcrel_result <- readRDS(rdsfile)
grm <- pcrelateToMatrix(pcrel_result)

cat("loading genotypes\n", file=stderr())
genoFile <- openfn.gds(filename=gdsfile, allow.fork=TRUE)
genoReader <- GdsGenotypeReader(filename=genoFile)
genoData <- GenotypeData(genoReader)
snpSplit <- split(getSnpID(genoData), 1:nsnp(genoData) %% 32)

run_assoc <- function (nullmod, endpoint) {

    cat("running SNP assocation\n", file=stderr())
    assoc <- do.call(rbind, mclapply(snpSplit,                
        function (snps) { 
            assocTestSingle(
                GenotypeBlockIterator(genoData, snpInclude=snps),
                null.model=nullmod, verbose=FALSE)
            }, mc.cores=32))
    
    assoc <- as.data.table(assoc)
    cat("found", assoc[Score.pval < 5e-8,.N],
        "GW significant loci", file=stderr())

    min_pval <- assoc[,min(Score.pval)]

    png(paste(endpoint, "GENESIS", "qq", "png", sep="."))
    qq(assoc[,Score.pval])
    dev.off()


    assoc[,chr_numeric := as.numeric(factor(chr, levels=c(1:22, "X", "XY")))]
    png(paste(endpoint, "GENESIS", "manhattan", "png", sep="."))
    manhattan(assoc, chr="chr_numeric", bp="pos", p="Score.pval", snp="variant.id", chrlabs=c(1:22,"X","XY"), col=c("blue4","orange3"),ylim=c(0,10), annotatePval=5e-8)
    dev.off()

    write.table(assoc, paste(endpoint, "GENESIS", "assoc", "txt", sep="."))
}

if (clinical_table[get(endpoint) != 0 & get(endpoint) != 1, .N] > 0) {
    cat("fitting model for", 
        endpoint,
        "with", clinical_table[!is.na(get(endpoint)),.N],
        "samples, ", clinical_table[get(endpoint) != 0,.N],
        "nonzero",
        fill=TRUE, file=stderr())
    nullmod <- fitNullModel(scan_annot, outcome=endpoint,
                                covars=covars,
                                cov.mat=grm, family=gaussian,
                                verbose=FALSE)
} else {
    cat("fitting model for", 
        endpoint,
        "with", clinical_table[get(endpoint),.N],
        "positives and ", clinical_table[(!get(endpoint)),.N],
        "negatives",
        fill=TRUE, file=stderr())
    nullmod <- fitNullModel(scan_annot, outcome=endpoint,
                            covars=covars,
                            cov.mat=grm, family=binomial,
                            verbose=FALSE)
}

run_assoc(nullmod, endpoint)

