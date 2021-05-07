library(data.table)
library(GENESIS)
library(GWASTools)
library(SNPRelate)
library(parallel)
library(qqman)

sample_file_basename <- "625_Samples.cohort.QC_filtered.sample_matched"
gdsfile <- paste(sample_file_basename, "maf_and_ld", "gds", sep=".")
rdsfile <- paste(sample_file_basename, "PCRelate", "RDS", sep=".")

clinical_table <- read.csv("../../data/covariates/clinical_data_deidentified_allsamples/jordad05/625_Samples.cohort.QC_filtered.sample_matched.age_flowcell_PCAir_dmatrix.csv", header=TRUE)
names(clinical_table)[1] <- "scanID"
covars <- c("age", "sex", "age_sex", #"is_hispanic",
               names(clinical_table)[startsWith(names(clinical_table), "pc")],
               #names(clinical_table)[startsWith(names(clinical_table), "race")],
               names(clinical_table)[startsWith(names(clinical_table), "flowcell")])

scan_annot <- ScanAnnotationDataFrame(clinical_table)
clinical_table <- as.data.table(clinical_table)

pcrel_result <- readRDS(rdsfile)
grm <- pcrelateToMatrix(pcrel_result)

binary_endpoints <- c("ever_covid", "ever_icu_covid_only", "deceased", "survived")
continuous_endpoints <- c("max_sofa", "max_sofa_icu_only", "sofa_range", "sofa_range_icu_only", "icu_hours", "icu_hours_with_zero", "max_viral_load_n1", "max_viral_load_n2", "max_viral_load_n3", "max_viral_load_cn1", "max_viral_load_cn2", "max_viral_load_cn3", "max_antibody_titer")
continuous_endpoints <- c(continuous_endpoints, paste(continuous_endpoints, "_irnt", sep=""))
bonferroni_threshold <- 5e-8 / (length(binary_endpoints) + length(continuous_endpoints))

genoFile <- openfn.gds(filename=gdsfile, allow.fork=TRUE)
genoReader <- GdsGenotypeReader(filename=genoFile)
genoData <- GenotypeData(genoReader)
snpSplit <- split(getSnpID(genoData), 1:nsnp(genoData) %% 32)


run_assoc <- function (nullmod, endpoint) {

    assoc <- do.call(rbind, mclapply(snpSplit,                
        function (snps) { 
            assocTestSingle(
                GenotypeBlockIterator(genoData, snpInclude=snps),
                null.model=nullmod, verbose=FALSE)
            }, mc.cores=32))
    
    assoc <- as.data.table(assoc)
    cat("found", assoc[Score.pval < 5e-8,.N],
        "GW significant loci,",
        assoc[Score.pval < bonferroni_threshold,.N],
        "passing Bonferroni correction", fill=TRUE)

    min_pval <- assoc[,min(Score.pval)]

    png(paste(endpoint, "GENESIS", "qq", "png", sep="."))
    qq(assoc[,Score.pval])
    dev.off()


    assoc[,chr_numeric := as.numeric(factor(chr, levels=c(1:22, "X", "XY")))]
    png(paste(endpoint, "GENESIS", "manhattan", "png", sep="."))
    manhattan(assoc, chr="chr_numeric", bp="pos", p="Score.pval", snp="variant.id", suggestiveline=-log10(5e-8), genomewideline=-log10(bonferroni_threshold), chrlabs=c(1:22,"X","XY"), col=c("blue4","orange3"),ylim=c(0,10), annotatePval=ifelse(min_pval < bonferroni_threshold, bonferroni_threshold, 5e-8))
    dev.off()

    write.table(assoc, paste(endpoint, "GENESIS", "assoc", "txt", sep="."))
}

for (endpoint in continuous_endpoints) {
    cat("fitting model for", 
        endpoint,
        "with", clinical_table[!is.na(get(endpoint)),.N],
        "samples, ", clinical_table[get(endpoint) != 0,.N],
        "nonzero",
        fill=TRUE)
    try({
        nullmod <- fitNullModel(scan_annot, outcome=endpoint,
                                covars=covars,
                                cov.mat=grm, family=gaussian,
                                verbose=FALSE)
        try(run_assoc(nullmod, endpoint))
    })
}

for (endpoint in binary_endpoints) {
    cat("fitting model for", 
        endpoint,
        "with", clinical_table[get(endpoint),.N],
        "positives and ", clinical_table[(!get(endpoint)),.N],
        "negatives",
        fill=TRUE)
    try({
        nullmod <- fitNullModel(scan_annot, outcome=endpoint,
                            covars=covars,
                            cov.mat=grm, family=binomial,
                            verbose=FALSE)
        run_assoc(nullmod, endpoint)
    })
}

