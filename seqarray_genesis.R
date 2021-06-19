suppressPackageStartupMessages({
    library(foreach)
    library(parallel)
    library(doParallel)
    library(GWASTools)
    library(SNPRelate)
    library(SeqArray)
    library(SeqVarTools)
    library(GENESIS)
    library(ggplot2)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(tidyverse)
    library(qqman)
    library(SummarizedExperiment)
})

cl <- NULL

start_cluster <- function() {
    hosts <- strsplit(Sys.getenv("LSB_HOSTS"), " ")[[1]]
    if (length(hosts) > 1) {
        cluster_hosts <- hosts[2:min(length(hosts), 65)]
        cl <<- makeCluster(cluster_hosts, rshcmd="blaunch")
    } else {
        num_cores <- length(mcaffinity())
        cl <<- makeCluster(num_cores - 1)
    }
    setDefaultCluster(cl)
    registerDoParallel(cl)
    seqParallelSetup(cl)
    1
}


stop_cluster <- function() {
    stopCluster(cl)
    1
}


convert_vcf2seqgds <- function(dependencies, targets) { 
    seqVCF2GDS(unlist(dependencies), unlist(targets), ignore.chr.prefix="")
    1
}


build_seq_gds <- function(dependencies, targets) {
    gdsFiles <- list.files(unlist(dependencies), "part-[0-9]{5}\\.seq\\.gds", full.names=TRUE)
    seqMerge(gdsFiles, unlist(targets))
    1
}

build_snp_gds <- function(prefix) { 
    bedfile <- paste(prefix, "bed", sep=".")
    bimfile <- paste(prefix, "bim", sep=".")
    famfile <- paste(prefix, "fam", sep=".")
    gdsfile <- paste(prefix, "snp", "gds", sep=".")
    snpgdsBED2GDS(bedfile, famfile, bimfile, gdsfile)
    1
}

run_pcair <- function(gdsfile, out_prefix) { 
    kingfile <- "625_Samples.cohort.QC_filtered.sample_matched.kin0"
    rdsfile <- paste(out_prefix, "PCAir", "RDS", sep=".")
    txtfile <- paste(out_prefix, "PCAir", "txt", sep=".")

    # load genotype data
    geno <- GdsGenotypeReader(filename=gdsfile)
    genoData <- GenotypeData(geno)

    # load KING result
    KINGmat <- kingToMatrix(kingfile, estimator="Kinship", sample.include=getVariable(genoData, "sample.id"))

    # run PC-AIR
    pcair_result <- pcair(genoData, kinobj=KINGmat, divobj=KINGmat)

    # Plot first 10 PCs
    pc_table <- as_tibble(pcair_result$vectors)
    pc_table %>% rename_all(~ str_replace(.x, 'V', 'PC')) %>% 
                 add_column(Subject_ID=pcair_result$sample.id, .before=1) -> pc_table
    clinical_covariates <- readRDS("../../data/covariates/clinical_data_deidentified_allsamples/Biobank_clinical_data_table_by_blood_sample_deidentified_UNCONSENTED.RDS")
    clinical_covariates %>% group_by(Subject_ID) %>% 
                            summarise(Race_From_Consent=unique(Race_From_Consent),
                                      Ethnicity_From_Consent=unique(Ethnicity_From_Consent)) -> race_table
    pc_table <- left_join(x=pc_table, y=race_table, by="Subject_ID")
    for (pc1 in 1:9) {
        for (pc2 in pc1:10) {
            pdf(paste(out_prefix, ".PC", pc1, "v", pc2, ".pdf", sep=""))
            ggplot(pc_table, aes_string(x=paste("PC", pc1, sep=""), 
                                   y=paste("PC", pc2, sep=""),
                                   color="Race_From_Consent"),
                                   shape="Ethnicity_From_Consent") + geom_point()
            dev.off()
        }
    }

    saveRDS(pcair_result, file=rdsfile)
    pc_table %>% select(c(Subject_ID, starts_with("PC"))) %>% write_delim(path=txtfile)
    1
}


run_pcrelate <- function(prefix) {
    gdsfile <- paste(prefix, "GWAS_filtered", "LD_pruned", "snp", "gds", sep=".")
    in_rdsfile <- paste(prefix, "PCAir", "RDS", sep=".")
    out_rdsfile <- paste(prefix, "PCRelate", "RDS", sep=".")

    pcair_result <- readRDS(in_rdsfile)
    geno <- GdsGenotypeReader(filename=gdsfile)
    genoData <- GenotypeData(geno)
    genoIterator <- GenotypeBlockIterator(genoData)
    pcrel_result <- pcrelate(genoIterator, pcs=pcair_result$vectors[,1:4], training.set=pcair_result$unrels)

    saveRDS(pcrel_result, file=out_rdsfile)
    1
}

genesis_null_model <- function(file_prefix, endpoint, covars) {
    clinical_table <- read.csv("/sc/arion/projects/mscic1/data/covariates/clinical_data_deidentified_allsamples/jordad05/625_Samples.cohort.QC_filtered.sample_matched.age_flowcell_PCAir_dmatrix.csv", header=TRUE)
    names(clinical_table)[1] <- "scanID"
    # covars <- c("age", "age_squared", "sex", "age_sex", "recruitment_date", #"is_hispanic",
    #                names(clinical_table)[startsWith(names(clinical_table), "pc")],
    #                names(clinical_table)[startsWith(names(clinical_table), "flowcell")])

    scan_annot <- ScanAnnotationDataFrame(clinical_table)
    clinical_table <- as_tibble(clinical_table)

    pcrel_result <- readRDS(paste(file_prefix, "PCRelate", "RDS", sep="."))
    grm <- pcrelateToMatrix(pcrel_result)
    nullmod <- tryCatch({
        if (clinical_table %>% filter(.[[endpoint]] != 0 & .[[endpoint]] != 1) %>% tally > 0) {
            fitNullModel(scan_annot, outcome=endpoint,
                         covars=covars,
                         cov.mat=grm, family="gaussian")
        } else {
            fitNullModel(scan_annot, outcome=endpoint,
                         covars=covars,
                         cov.mat=grm, family="binomial")
            }
        }, error = function(c) list(converged=FALSE) )

    saveRDS(nullmod, paste(file_prefix, endpoint, "null", "RDS", sep="."))
    if (nullmod$converged) {
        1
    } else {
        0
    }
}

genesis_null_model_ffs <- function(file_prefix, endpoint) {
    read_csv("/sc/arion/projects/mscic1/data/covariates/clinical_data_deidentified_allsamples/jordad05/625_Samples.cohort.QC_filtered.sample_matched.age_flowcell_PCAir_dmatrix.csv") %>% 
        rename(scanID=X1) -> clinical_table
    scan_annot <- ScanAnnotationDataFrame(as.data.frame(clinical_table))  # somehow GWASTools doesn't recognize tibble columns?

    clinical_table %>% names %>% str_subset("^flowcell") %>% list %>%
                        append(list("age",
                                    "age_squared",
                                    "sex",
                                    "age_sex",
                                    "recruitment_date",
                                    "is_hispanic",
                                    "race_factor",
                                    "multi_batch",
                                    "any_comorbidity")) %>%
                         append(map(1:10, ~paste0("pc", 1:.x))) -> all_covars

    family <- if (clinical_table %>% filter(.[[endpoint]] != 0 & .[[endpoint]] != 1) %>% tally > 0) "gaussian" else "binomial"

    pcrel_result <- readRDS(paste(file_prefix, "PCRelate", "RDS", sep="."))
    grm <- pcrelateToMatrix(pcrel_result)

    # do forward feature selection
    last_model <- list(AIC=Inf, covars=NULL)
    current_model <- tryCatch(
      fitNullModel(scan_annot, outcome=endpoint,
                               covars=NULL,
                               cov.mat=grm, family=family),
      error=function(e) list(AIC=.Machine$double.xmax,  # smaller than Inf, larger than any other number
                              covars=NULL)
    )
    while (current_model$AIC < last_model$AIC) {
        cat("feature selection round", length(current_model$covars) + 1, "\n")
        last_model <- current_model
        models <- foreach(covar=setdiff(all_covars, current_model$covars),
                         .packages=c("GENESIS", "magrittr", "purrr"), .errorhandling="remove") %dopar% {
            cat(covar, "\n")
            this_loop_covars <- append(current_model$covars, list(covar))
            model <- fitNullModel(scan_annot, outcome=endpoint,
                         covars=this_loop_covars %>% flatten_chr %>% unique,
                         cov.mat=grm, family=family)
            model$covars <- this_loop_covars
            model
        }
        if (length(models) == 0) break
        models %>% map("AIC") %>% compact %>% which.min %>% pluck(models, .) -> current_model
        cat("Selected:", last(current_model$covars), "\nAIC:", current_model$AIC, "\n")
    }

    saveRDS(last_model, paste(file_prefix, endpoint, "null", "RDS", sep="."))

    if (!is.null(last_model$converged) && last_model$converged) {
        1
    } else {
        0
    }
}


genesis_null_model_exhaustive <- function(file_prefix, endpoint) {
    read_csv("/sc/arion/projects/mscic1/data/covariates/clinical_data_deidentified_allsamples/jordad05/625_Samples.cohort.QC_filtered.sample_matched.age_flowcell_PCAir_dmatrix.csv") %>%
      rename(scanID=X1) %>% mutate(race_factor = factor(race_factor)) -> clinical_table
    scan_annot <- ScanAnnotationDataFrame(as.data.frame(clinical_table))  # somehow GWASTools doesn't recognize tibble columns?

    clinical_table %>% names %>% str_subset("^flowcell") %>% list %>%
      append(list("age",
                  "age_squared",
                  "sex",
                  "age_sex",
                  "recruitment_date",
                  "is_hispanic",
                  "race_factor",
                  "multi_batch",
                  "any_comorbidity")) -> non_pc_covars

    map(1:10, ~paste0("pc", 1:.x)) %>% prepend(list(character())) -> pcs
    map(1:length(non_pc_covars), ~ combn(non_pc_covars, .x, simplify=FALSE)) %>%
      flatten %>% cross2(pcs) %>% map(unlist) -> all_covar_combinations

    family <- if (clinical_table %>% filter(.[[endpoint]] != 0 & .[[endpoint]] != 1) %>% tally > 0) "gaussian" else "binomial"

    pcrel_result <- readRDS(paste(file_prefix, "PCRelate", "RDS", sep="."))
    grm <- pcrelateToMatrix(pcrel_result)

    ptm <- proc.time()
    models <- foreach (covars=all_covar_combinations,
             .packages="GENESIS",
             .errorhandling="remove") %dopar% {
            model <- fitNullModel(scan_annot, outcome=endpoint,
                                  covars=covars,
                                  cov.mat=grm, family=family)
            model$covars <- covars
            model
        }
    exec_time <- proc.time() - ptm
    cat(length(models), "of", length(all_covar_combinations), "covar combinations tested in", exec_time[["elapsed"]], "seconds\n")
    map_dbl(models, "AIC") %>% compact %>% which.min %>% pluck(models, .) -> model
    cat("Best combination:", model$covars, "\n")

    saveRDS(model, paste(file_prefix, endpoint, "null", "RDS", sep="."))
    model
    #if (!is.null(model$converged) && model$converged) {
    #    1
    #} else {
    #    0
    #}
}

gwas_analytics <- function(endpoints) {
    covar_column_names <- c(race="race_factor", ethnicity="is_hispanic", flowcell="flowcell_HFHFYDSXY",
                            multi_batch="multi_batch", recruitment_date="recruitment_date", age="age",
                            sex="sex", age_sex="age_sex", age_squared="age_squared", any_comorbidity="any_comorbidity")
    output_table <- foreach (trait=endpoints,
             .errorhandling="remove",
             .packages=c("tibble", "readr", "purrr"),
             .combine=bind_rows,
             .multicombine=TRUE) %dopar% {
        assoc <- read_delim(paste(trait, "GENESIS", "assoc", "txt", sep="."), delim=" ")
        lambda <- median(assoc$Score.Stat^2) / qchisq(0.5, 1)
        hits <- sum(assoc$Score.pval < 5e-8)

        nullmod <- readRDS(paste("625_Samples.cohort.QC_filtered.sample_matched", trait, "null.RDS", sep="."))
        n_pcs <- sum(paste0("pc", 1:10) %in% nullmod$covars)
        covars_present <- map(covar_column_names, ~ .x %in% nullmod$covars)  # map preserves names
        tibble(trait=trait,
               lambda=lambda,
               gwas_hits=hits,
               pcs=n_pcs,
               !!!covars_present)
    }
    write_tsv(output_table, "GENESIS.GWAS_report.tsv")
    1
}

make_gwas_plots <- function(endpoint) {
    assoc <- read_delim(paste(endpoint, "GENESIS", "assoc", "txt", sep="."), delim=" ")

    png(paste(endpoint, "GENESIS", "qq", "png", sep="."))
    qq(assoc$Score.pval)
    dev.off()

    png(paste(endpoint, "GENESIS", "manhattan", "png", sep="."))
    mutate(assoc, chr_numeric=as.numeric(as.factor(chr))) %>% filter(chr_numeric < 25) %>%
        manhattan(chr="chr_numeric", bp="pos", p="Score.pval", snp="variant.id",
                    chrlabs=c(1:22,"X","Y"), col=c("blue4","orange3"))#,ylim=c(0,10))
    dev.off()
    1
}

run_smmat <- function(gds_filename, null_model_filename, phenotype_name) { 
    geneRanges <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
    nullmod <- readRDS(null_model_filename)
    gdsFile <- seqOpen(gds_filename)
    tryCatch({
             seqData <- SeqVarData(gdsFile)
             seqIter <- SeqVarRangeIterator(seqData, geneRanges)

             assoc <- assocTestAggregate(seqIter, nullmod, sparse=FALSE, test="SMMAT")
             as_tibble(assoc$results, rownames="geneid") %>% add_column(chr=as.factor(seqnames(geneRanges)),
                                                                        start=start(geneRanges),
                                                                        end=end(geneRanges), .after="geneid")
    }, finally=seqClose(gdsFile)) -> assoc
    write_delim(assoc, paste(phenotype_name, 'GENESIS', 'SMMAT', 'assoc', 'txt', sep='.'))
    
    library(qqman)
    png(paste(phenotype_name, 'GENESIS', 'SMMAT', 'qq', 'png', sep="."))
    qq(assoc$pval_SMMAT)
    dev.off()

    mutate(assoc, chr_numeric=as.numeric(chr)) %>% filter(chr_numeric < 24) -> assoc_manhattan
    png(paste(phenotype_name, "GENESIS", "SMMAT", "manhattan", "png", sep="."))
    manhattan(assoc_manhattan, chr="chr_numeric", bp="start", p="pval_SMMAT", chrlabs=c(1:22,"X"), 
              col=c("blue4","orange3"), genomewideline=-log10(0.05 / length(geneRanges)), ylim=c(0,10))
    dev.off()
    1
}

make_isoform_tpms_table <- function (infile, outfile) {
    transcript_abundance_experiment <- readRDS(infile)
    assay(transcript_abundance_experiment, "abundance") %>%
      as_tibble %>% rename_all(~ paste0("WholeBlood.", 1:length(.x))) %>%
      add_column(transcript_id=rowData(transcript_abundance_experiment)$tx_id,
                 gene_id=rowData(transcript_abundance_experiment)$gene_id,
                 .before=1) %>% write_tsv(outfile)
    1
}

if (exists("snakemake")) {
    switch(
        snakemake@params[["genesis_cmd"]],
        plink2snpgds=build_snp_gds(snakemake@wildcards[["prefix"]]),
        pcair=run_pcair(snakemake@input[["gds"]], snakemake@params[["output_stem"]]),
        pcrelate=run_pcrelate(snakemake@params[["prefix"]]),
        gwas_plots=make_gwas_plots(snakemake@wildcards[["phenotype"]]),
        isoform_table=make_isoform_tpms_table(snakemake@input[[1]], snakemake@output[[1]])
    )
}